#!/bin/bash

set -euo pipefail

_prepend_unique_path() {
  local var_name=$1
  local dir_path=$2
  if [[ -d "$dir_path" ]]; then
    local current=""
    if [[ -n ${!var_name+x} ]]; then
      current="${!var_name}"
    fi
    case ":${current}:" in
      *":${dir_path}:"*) ;;
      *)
        if [[ -n "$current" ]]; then
          local new_value="${dir_path}:${current}"
          printf -v "$var_name" '%s' "$new_value"
        else
          printf -v "$var_name" '%s' "$dir_path"
        fi
        export "$var_name"
        ;;
    esac
  fi
}

# Ensure shared libraries from user and linuxbrew prefixes are visible
_prepend_unique_path LD_LIBRARY_PATH "$HOME/usr/local/lib"
_prepend_unique_path LD_LIBRARY_PATH "$HOME/.linuxbrew/lib"

normalize_bcf() {
  local bcf_path=$1
  if [[ ! -f $bcf_path ]]; then
    echo "normalize_bcf: missing file $bcf_path" >&2
    return 1
  fi
  SSH_AUTH_SOCK= bcftools view -Ov "$bcf_path" |
    grep -Ev '^##(fileDate|bcftools_view(Command|Version))='
}

assert_same_variants() {
  local actual_bcf=$1
  local expected_vcf=$2
  if [[ ! -f $expected_vcf ]]; then
    echo "assert_same_variants: missing expected VCF $expected_vcf" >&2
    return 1
  fi
  local tmp_dir
  tmp_dir=$(mktemp -d)
  trap 'rm -rf "$tmp_dir"' RETURN
  local normalized_actual=$tmp_dir/actual.vcf
  normalize_bcf "$actual_bcf" >"$normalized_actual"
  diff -u "$expected_vcf" "$normalized_actual"
}
