#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "$0")" && pwd)
resolver="$script_dir/resolve-bioc.sh"
fixture_dir=$(mktemp -d)
trap 'rm -rf "$fixture_dir"' EXIT
mkdir "$fixture_dir/bin"
cat > "$fixture_dir/bin/curl" <<'CURL'
#!/usr/bin/env bash
set -euo pipefail
for url; do :; done
case "$url" in
  https://bioconductor.org/config.yaml)
    [[ ${FAIL_CONFIG:-0} == 0 ]] || exit 22
    cat "$FIXTURE_DIR/config.yaml" ;;
  https://api.r-hub.io/rversions/r-release)
    [[ ${FAIL_RELEASE:-0} == 0 ]] || exit 22
    cat "$FIXTURE_DIR/release.json" ;;
  *) exit 22 ;;
esac
CURL
chmod +x "$fixture_dir/bin/curl"

write_fixture() {
  printf 'devel_version: "3.24"\nr_version_associated_with_release: "4.6.0"\nr_version_associated_with_devel: "%s"\n' \
    "$1" > "$fixture_dir/config.yaml"
  printf '{"version":"%s"}\n' "$2" > "$fixture_dir/release.json"
}

run_resolver() {
  : > "$fixture_dir/output"
  env PATH="$fixture_dir/bin:$PATH" FIXTURE_DIR="$fixture_dir" \
    GITHUB_OUTPUT="$fixture_dir/output" \
    FAIL_CONFIG="${FAIL_CONFIG:-0}" FAIL_RELEASE="${FAIL_RELEASE:-0}" \
    bash "$resolver" > "$fixture_dir/log" 2>&1
}

check_selection() {
  local name=$1 required=$2 released=$3 expected=$4 matrix minor
  printf 'Checking %s\n' "$name"
  write_fixture "$required" "$released"
  run_resolver
  matrix=$(cat "$fixture_dir/output")
  [[ $matrix == matrix=* ]]
  minor=$(printf '%s\n' "$required" | cut -d. -f1,2)
  printf '%s\n' "${matrix#matrix=}" | jq -e \
    --arg selector "$expected" --arg minor "$minor" '
      (.config | length) == 3 and
      [.config[].os] == ["linux", "macos", "windows"] and
      [.config[].runner] == ["ubuntu-latest", "macOS-latest", "windows-latest"] and
      all(.config[]; .r == $selector and .r_minor == $minor and .bioc == "3.24") and
      .config[0].cont == "bioconductor/bioconductor_docker:devel" and
      (.config[1] | has("cont") | not) and (.config[2] | has("cont") | not)
    ' >/dev/null
}

check_failure() {
  printf 'Checking %s fails without emitting a matrix\n' "$1"
  if run_resolver; then
    printf 'Unexpected success: %s\n' "$1" >&2
    exit 1
  fi
  [[ ! -s "$fixture_dir/output" ]]
}

check_selection released 4.6.0 4.6.2 4.6
check_selection older_released_minor 4.5.0 4.6.2 4.5
check_selection unreleased 4.7.0 4.6.2 devel
check_selection release_transition 4.7.0 4.7.0 4.7
check_selection numeric_minor_order 4.10.0 4.9.9 devel

write_fixture 4.6.0 4.6.2
printf 'devel_version: "3.24"\ndevel_version: "3.25"\nr_version_associated_with_devel: "4.6.0"\n' \
  > "$fixture_dir/config.yaml"
check_failure duplicate_config_scalar
printf 'devel_version:3.24\nr_version_associated_with_devel: "4.6.0"\n' \
  > "$fixture_dir/config.yaml"
check_failure malformed_scalar_separator
printf 'devel_version: "3.24"\nr_version_associated_with_devel: ""\n' \
  > "$fixture_dir/config.yaml"
check_failure empty_r_version
write_fixture 4.6.0 4.6.2
printf '{"version":"not-a-version"}\n' > "$fixture_dir/release.json"
check_failure malformed_release_version
write_fixture 4.6.0 4.6.2
FAIL_CONFIG=1 check_failure config_network_failure
FAIL_RELEASE=1 check_failure release_network_failure
printf 'Resolver checks completed.\n'
