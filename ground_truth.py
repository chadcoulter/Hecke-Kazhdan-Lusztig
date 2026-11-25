#!/usr/bin/env python3
"""Ground-truth probe — run this BEFORE any panel dispatch.

Checks the locked anchors A1-A9 directly against primary arithmetic.
Any FAIL blocks the panel.
"""
from __future__ import annotations

import sys
sys.path.insert(0, "/home/sudoshellz/Documents/hecke-kl")

from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl.hecke import HeckeAlgebra
from hecke_kl import poly

PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"{'PASS' if cond else 'FAIL'}  {name}  {detail}")


# A1: rank<=2: all KL polynomials 0 or 1
for typ, n in [("A2", 2), ("I2", 4), ("I2", 6)]:
    W = make(typ, n)
    T = KLTable(W)
    vals = T.all_polynomials()
    bad = [(W.describe(y), W.describe(w), poly.fmt(p)) for (y, w), p in vals.items()
           if p not in ({}, {0: 1})]
    check(f"A1 rank<=2 {W.name} all P in {{0,1}}", not bad, f"order={W.order} bad={bad[:3]}")

# A2: l(w)-l(y) in {0,1,2} => P=1
for typ, n in [("A", 4), ("B", 3)]:
    W = make(typ, n)
    T = KLTable(W)
    bad = []
    for w in W.elements:
        for y in W.elements:
            d = W.length(w) - W.length(y)
            if d in (0, 1, 2) and W.leq(y, w):
                if not poly.eq(T.P(y, w), poly.one()):
                    bad.append((W.describe(y), W.describe(w), poly.fmt(T.P(y, w))))
    check(f"A2 {W.name} diff 0/1/2 => P=1", not bad, f"bad={bad[:3]}")

# A3: w0 => P=1 for all y<=w
for typ, n in [("A", 4), ("B", 3), ("G2", None)]:
    W = make(typ, n)
    T = KLTable(W)
    bad = [(W.describe(y), poly.fmt(T.P(y, W.longest)))
           for y in W.elements if W.leq(y, W.longest)
           and not poly.eq(T.P(y, W.longest), poly.one())]
    check(f"A3 {W.name} w0 => all P=1", not bad, f"bad={bad[:3]}")

# A4: S4 nontrivial anchors  P_{b,bacb}=1+q  P_{ac,acbca}=1+q
W = make("A", 4)
T = KLTable(W)
cases = [
    ((1,), (1, 0, 2, 1), "P_{b,bacb}", [1, 1]),
    ((0, 2), (0, 2, 1, 2, 0), "P_{ac,acbca}", [1, 1]),
]
for yw, ww, label, expected in cases:
    y, w = W.from_word(yw), W.from_word(ww)
    got = poly.as_q_list(T.P(y, w))
    check(f"A4 S4 {label} = 1+q", got == expected,
          f"got={got} expected={expected} ylen={W.length(y)} wlen={W.length(w)}")

# A5: nonnegativity
for typ, n in [("A", 5), ("B", 3), ("G2", None)]:
    W = make(typ, n)
    T = KLTable(W)
    vals = T.all_polynomials()
    bad = [(W.describe(y), W.describe(w), p) for (y, w), p in vals.items()
           if any(c < 0 for c in p.values())]
    check(f"A5 {W.name} nonnegative coefficients", not bad, f"bad={bad[:2]}")

# A7: degree bound and constant term
for typ, n in [("A", 5), ("B", 3)]:
    W = make(typ, n)
    T = KLTable(W)
    bad = []
    for (y, w), p in T.all_polynomials().items():
        if y == w or not W.leq(y, w):
            continue
        d = W.length(w) - W.length(y)
        if poly.qdegree(p) > (d - 1) // 2 or p.get(0, 0) != 1:
            bad.append((W.describe(y), W.describe(w), poly.fmt(p)))
    check(f"A7 {W.name} degree bound + constant term 1", not bad, f"bad={bad[:3]}")

# A8: bar invariance of C'_w
for typ, n in [("A", 4), ("B", 3), ("G2", None)]:
    W = make(typ, n)
    T = KLTable(W)
    H = HeckeAlgebra(W, T)
    bad = []
    for w in W.elements:
        c = H.Cprime(w)
        if not H.equal(H.bar(c), c):
            bad.append(W.describe(w))
    check(f"A8 {W.name} bar(C'_w)=C'_w", not bad, f"bad={bad[:5]}")

# A9: R identity  T_{y^{-1}}^{-1} = sum_x bar(R_{x,y}) q^{-l(x)} T_x
for typ, n in [("A", 4), ("B", 3)]:
    W = make(typ, n)
    T = KLTable(W)
    H = HeckeAlgebra(W, T)
    bad = []
    for y in W.elements:
        lhs = H.inverse_Tw(W.inverse(y))
        rhs = {}
        for x in W.elements:
            r = T.R(x, y)
            if not r:
                continue
            br = poly.bar(r)
            shift = poly.shift(br, -2 * W.length(x))
            rhs = H.add(rhs, {x: shift})
        if not H.equal(lhs, rhs):
            bad.append(W.describe(y))
    check(f"A9 {W.name} R identity", not bad, f"bad={bad[:5]}")

# A6: Hecke associativity + quadratic relation
W = make("A", 4)
H = HeckeAlgebra(W, KLTable(W))
bad = []
for i in range(W.n):
    s = W.S[i]
    lhs = H.mul(H.T(s), H.T(s))
    rhs = H.add({W.identity: poly.q()}, H.scale(H.T(s), poly.sub(poly.q(), poly.one())))
    if not H.equal(lhs, rhs):
        bad.append(i)
check("A6 quadratic relation T_s^2 = q + (q-1)T_s", not bad, f"bad={bad}")

print(f"\nSUMMARY  PASS={len(PASS)}  FAIL={len(FAIL)}")
if FAIL:
    print("FAILED:", FAIL)
    sys.exit(1)
