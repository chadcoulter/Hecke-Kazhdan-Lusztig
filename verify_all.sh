#!/bin/sh
# Rebuild-from-clean + run every verification claim. One command, no screenshots.
set -e
cd "$(dirname "$0")"
echo "=== ground truth (21 anchors) ==="
python3 ground_truth.py
echo "=== Bruhat closure vs subword criterion ==="
python3 check_bruhat.py
echo "=== independent R-inversion cross-check ==="
python3 check_rinvert.py
echo "=== Hecke inverse ==="
python3 check_inverse.py
echo "=== pytest ==="
python3 -m pytest tests/ -q
echo "=== Koppa-sealed manifest ==="
python3 examples/make_manifest.py
echo "ALL VERIFICATIONS PASSED"
