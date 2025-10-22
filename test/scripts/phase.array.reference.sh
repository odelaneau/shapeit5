#!/bin/bash
set -euo pipefail

SCRIPT_PATH=$(realpath "${BASH_SOURCE[0]}")
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
TEST_DIR=$(realpath "${SCRIPT_DIR}/..")
cd "$TEST_DIR"

source "$SCRIPT_DIR/lib/test_utils.sh"

tmp_dir=$(mktemp -d)
trap 'rm -rf "$tmp_dir"' EXIT

output_bcf="$tmp_dir/target.phased.bcf"
region="${TEST_REGION:-1:5000000-6000000}"

../phase_common/bin/phase_common \
  --input array/target.unrelated.bcf \
  --reference array/reference.bcf \
  --region "$region" \
  --map info/chr1.gmap.gz \
  --output "$output_bcf" \
  --thread 8

assert_same_variants "$output_bcf" "$SCRIPT_DIR/expected/phase.array.reference.vcf"
