#!/usr/bin/env bash
# the regression gate: format, clippy, and the whole test suite must pass.
# exits non-zero if any check fails. run from the repo root.
set -u

export PATH="$HOME/.cargo/bin:$PATH"
here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$here"

fails=0

echo "== fmt --check =="
cargo fmt --check >/dev/null 2>&1 || { echo "FAIL: formatting"; fails=$((fails+1)); }

echo "== clippy -D warnings =="
cargo clippy --workspace --all-targets -- -D warnings >/tmp/reg_clippy 2>&1 || { echo "FAIL: clippy"; fails=$((fails+1)); }

echo "== test --workspace =="
suites=$(cargo test --workspace 2>&1 | grep -cE "test result: ok")
failed=$(cargo test --workspace 2>&1 | grep -cE "FAILED|error\[")
echo "  passing suites: $suites, failing/errors: $failed"
if [ "$suites" -eq 0 ] || [ "$failed" -ne 0 ]; then
  echo "FAIL: tests"; fails=$((fails+1))
fi

echo "== repo hygiene =="
if git status --porcelain | grep -q .; then echo "FAIL: dirty tree"; fails=$((fails+1)); fi
for p in plans/ goals/ nuFor_README.md docs/research docs-site; do
  if git ls-files | grep -q "^$p"; then echo "FAIL: tracked private path $p"; fails=$((fails+1)); fi
done

if [ "$fails" -eq 0 ]; then
  echo "regression gate: PASS ($suites suites green)"
  exit 0
fi
echo "regression gate: FAIL ($fails check(s) failing)"
exit 1