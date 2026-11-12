#!/usr/bin/env bash
set -euo pipefail

fail() {
  printf 'Bioconductor matrix: %s\n' "$*" >&2
  exit 1
}

[[ -n ${GITHUB_OUTPUT:-} ]] || fail 'GITHUB_OUTPUT is required.'
command -v jq >/dev/null || fail 'jq is required.'

bioc_config=$(curl --fail --silent --show-error --location --retry 2 \
  https://bioconductor.org/config.yaml)
r_release_json=$(curl --fail --silent --show-error --location --retry 2 \
  https://api.r-hub.io/rversions/r-release)

# These are simple top-level scalars in Bioconductor's published configuration.
version_scalar() {
  local value
  value=$(printf '%s\n' "$bioc_config" | awk -v key="$1" '
    $0 ~ ("^" key ":[[:space:]]") {
      count++
      value = substr($0, length(key) + 2)
      sub(/[[:space:]]+#.*$/, "", value)
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
    }
    END { if (count != 1) exit 1; print value }
  ') || fail "Expected one $1 scalar."
  case "$value" in
    \"*\") value=${value#\"}; value=${value%\"} ;;
    \'*\') value=${value#\'}; value=${value%\'} ;;
  esac
  [[ $value =~ ^[0-9]+\.[0-9]+(\.[0-9]+)?$ ]] ||
    fail "$1 must be a numeric version."
  printf '%s\n' "$value"
}

bioc=$(version_scalar devel_version)
[[ $bioc =~ ^[0-9]+\.[0-9]+$ ]] || fail 'devel_version must have two components.'
required_r=$(version_scalar r_version_associated_with_devel)
released_r=$(printf '%s\n' "$r_release_json" | jq -er '
  .version | select(type == "string") |
  select(test("^[0-9]+\\.[0-9]+\\.[0-9]+$"))
') || fail 'The published R release must contain a numeric version.'

IFS=. read -r required_major required_minor required_patch <<< "$required_r"
IFS=. read -r released_major released_minor released_patch <<< "$released_r"
r_minor="$required_major.$required_minor"
if (( 10#$required_major < 10#$released_major ||
      (10#$required_major == 10#$released_major &&
       10#$required_minor <= 10#$released_minor) )); then
  r_selector=$r_minor
else
  r_selector=devel
fi

matrix=$(jq -cn --arg r "$r_selector" --arg r_minor "$r_minor" --arg bioc "$bioc" '
  {config: [
    {os: "linux", runner: "ubuntu-latest",
     cont: "bioconductor/bioconductor_docker:devel"},
    {os: "macos", runner: "macOS-latest"},
    {os: "windows", runner: "windows-latest"}
  ] | map(. + {r: $r, r_minor: $r_minor, bioc: $bioc})}
')
printf 'matrix=%s\n' "$matrix" >> "$GITHUB_OUTPUT"
printf 'Bioconductor %s; required R %s; published R %s; setup-r %s.\n' \
  "$bioc" "$r_minor" "$released_r" "$r_selector"
