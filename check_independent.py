#!/usr/bin/env python3
"""Independent KL cross-check: recover P_{y,w} by solving bar(C'_w)=C'_w
with sympy, WITHOUT the KL recursion. Compare to KLTable.P on small groups.

This is a different route: the recursion uses mu coefficients from shorter w;
the linear system uses only the bar-involution and the T-basis multiplication
table. Agreement is nontrivial evidence both are correct.
"""
import sys
sys.path.insert(0, __import__("pathlib").Path(__file__).resolve().parent)
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl.hecke import HeckeAlgebra
from hecke_kl import poly
import sympy as sp

for typ, n in [("A", 3), ("B", 2), ("I2", 4)]:
    W = make(typ, n)
    T = KLTable(W)
    H = HeckeAlgebra(W, T)
    v = sp.Symbol("v")
    bad = []
    checked = 0
    for w in W.elements:
        recovered = H.solve_kl_bar(w)
        for y in W.elements:
            if not W.leq(y, w):
                continue
            # P_{y,w}: constant term 1, higher coeffs from recovered
            cs = [1]
            d = W.length(w) - W.length(y)
            for i in range(1, (d - 1) // 2 + 1):
                cs.append(recovered.get((y, i), poly.zero()).get(0, 0))
            expected = poly.qpoly_from_coeffs(cs)
            got = T.P(y, w)
            checked += 1
            if not poly.eq(got, expected):
                bad.append((W.describe(y), W.describe(w), poly.fmt(got), poly.fmt(expected)))
    print(f"{W.name}: checked {checked} pairs, mismatches {len(bad)}", bad[:5], flush=True)
    if bad:
        sys.exit(1)
print("PASS: KL recursion agrees with independent bar-invariance linear system")
