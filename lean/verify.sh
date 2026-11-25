#!/usr/bin/env bash
# Verify the S1 Lean 4 proof of the Hecke–KL group layer.
# Requires Lean 4.36 on PATH (via ~/.elan/bin).  No Mathlib needed (Std only).
set -euo pipefail
export PATH="$HOME/.elan/bin:$PATH"
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

out=$(mktemp)
trap 'rm -f "$out"' EXIT
if lean group.lean >"$out" 2>&1; then
  echo "PASS  group.lean  ($(grep -cE '^(theorem|example) ' group.lean) proved declarations)"
else
  echo "FAIL  group.lean"
  cat "$out"
  exit 1
fi
if grep -qE '\b( axiom |sorry|admit|unsafe)' group.lean; then
  echo "FAIL  forbidden token found"
  grep -nE '\b( axiom |sorry|admit|unsafe)' group.lean
  exit 1
fi
echo "PASS  no axiom/sorry/admit/unsafe"
