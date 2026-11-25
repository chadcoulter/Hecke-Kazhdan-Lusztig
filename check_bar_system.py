#!/usr/bin/env python3
"""Regression for S2's F1: solve_kl_bar must not crash for any w, and its
recovered coefficients must match KLTable.P (the independent cross-check the
contract requires)."""
import sys
sys.path.insert(0, __import__("pathlib").Path(__file__).resolve().parent)
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl.hecke import HeckeAlgebra
from hecke_kl import poly

for typ, n in [("A", 3), ("A", 4), ("B", 2), ("G2", None)]:
    W = make(typ, n)
    T = KLTable(W)
    H = HeckeAlgebra(W, T)
    checked, bad = 0, []
    for w in W.elements:
        recovered = H.solve_kl_bar(w)  # must not crash
        for y in W.elements:
            if not W.leq(y, w) or y == w:
                continue
            d = W.length(w) - W.length(y)
            cs = [1] + [recovered.get((y, i), poly.zero()).get(0, 0)
                        for i in range(1, (d - 1) // 2 + 1)]
            expected = poly.qpoly_from_coeffs(cs)
            got = T.P(y, w)
            checked += 1
            if not poly.eq(got, expected):
                bad.append((W.describe(y), W.describe(w), poly.fmt(got), poly.fmt(expected)))
    print(f"{W.name}: checked {checked} pairs, mismatches {len(bad)}", bad[:5], flush=True)
    if bad:
        sys.exit(1)
print("PASS: solve_kl_bar runs without crashing and agrees with KLTable.P")
