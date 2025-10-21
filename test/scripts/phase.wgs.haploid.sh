#!/bin/bash
set -euo pipefail

SCRIPT_PATH=$(realpath "${BASH_SOURCE[0]}")
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
TEST_DIR=$(realpath "${SCRIPT_DIR}/..")
cd "$TEST_DIR"

source "$SCRIPT_DIR/lib/test_utils.sh"

tmp_dir=$(mktemp -d)
trap 'rm -rf "$tmp_dir"' EXIT

region="${TEST_REGION:-1:5000000-6000000}"
chunk_file="${TEST_CHUNK_FILE:-info/chunks.coordinates.small.txt}"
if [[ ! -f "$chunk_file" ]]; then
  chunk_file="info/chunks.coordinates.txt"
fi

HAP=info/target.haploid.txt
list_file="$tmp_dir/chunks.files.txt"
scaffold_bcf="$tmp_dir/target.scaffold.bcf"
: >"$list_file"

../phase_common/bin/phase_common \
  --input wgs/target.haploid.bcf \
  --filter-maf 0.001 \
  --haploids "$HAP" \
  --region "$region" \
  --map info/chr1.gmap.gz \
  --output "$scaffold_bcf" \
  --thread 8

while read -r CHK CHR SRG IRG; do
  [[ -z "$CHK" ]] && continue
  OUT="$tmp_dir/target.phased.chunk${CHK}.bcf"
  log_file="$tmp_dir/phase_rare_${CHK}.log"
  if ../phase_rare/bin/phase_rare \
      --input wgs/target.haploid.bcf \
      --scaffold "$scaffold_bcf" \
      --haploids "$HAP" \
      --map info/chr1.gmap.gz \
      --input-region "$IRG" \
      --scaffold-region "$SRG" \
      --output "$OUT" \
      --thread 8 >"$log_file" 2>&1; then
    if [[ -f "$OUT" ]]; then
      echo "$OUT" >>"$list_file"
    fi
  else
    if grep -q "No variants to be phased" "$log_file"; then
      continue
    fi
    cat "$log_file" >&2
    exit 1
  fi
done <"$chunk_file"

output_bcf="$tmp_dir/target.phased.bcf"
if [[ -s "$list_file" ]]; then
  SSH_AUTH_SOCK= bcftools concat -n -Ob -o "$output_bcf" -f "$list_file"
  SSH_AUTH_SOCK= bcftools index "$output_bcf"
else
  cp "$scaffold_bcf" "$output_bcf"
  cp "${scaffold_bcf}.csi" "${output_bcf}.csi"
fi

assert_same_variants "$output_bcf" "$SCRIPT_DIR/expected/phase.wgs.haploid.vcf"
