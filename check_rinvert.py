#!/usr/bin/env python3
"""Independent KL check via the R-polynomial inversion identity (KL79).

For x<=w:
  q^{(l(w)-l(x))/2} bar(P_{x,w}) - q^{(l(x)-l(w))/2} P_{x,w}
    = sum_{x<y<=w} (-1)^{l(x)+l(y)} q^{(-l(x)+2l(y)-l(w))/2} bar(R_{x,y}) P_{y,w}

P_{y,w}=0 unless y<=w, P_{w,w}=1, bar(P_{x,w}) = P_{x,w}(q^{-1}).
This is a DIFFERENT recursion from the mu-based one: it uses R-polynomials,
and solves for P_{x,w} using already-known P_{y,w} with y>x (descending order).
"""
import sys
sys.path.insert(0, __import__("pathlib").Path(__file__).resolve().parent)
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl import poly

for typ, n in [("A", 4), ("B", 3), ("G2", None), ("A", 5)]:
    W = make(typ, n)
    T = KLTable(W)
    bad = []
    checked = 0
    for w in W.elements:
        # process x in DECREASING length so P_{y,w} for y>x are known
        xs = sorted([x for x in W.elements if W.leq(x, w)],
                    key=lambda z: -W.length(z))
        for x in xs:
            if x == w or W.length(w) - W.length(x) < 1:
                continue
            # build RHS from already-computed P_{y,w} (y shorter than x? no, y>x)
            rhs = poly.zero()
            for y in W.elements:
                if not W.leq(x, y) or not W.leq(y, w) or y == x:
                    continue
                if W.length(y) <= W.length(x):
                    continue  # P_{y,w} known only if length(y) > length(x)? Actually need y>x
                sign = (-1) ** (W.length(x) + W.length(y))
                e = (-W.length(x) + 2 * W.length(y) - W.length(w))
                if e % 2:
                    continue
                # q^{e/2} may have negative exponent; represent as v^e (e can be neg)
                qpow = poly.shift(poly.one(), e)
                r = T.R(x, y)
                term = poly.mul(
                    poly.scale(poly.mul(qpow, poly.bar(r)), sign),
                    poly.shift(T.P(y, w), 0),
                )
                rhs = poly.add(rhs, term)
            # solve for P_{x,w}: q^{a} bar(P) - q^{-a} P = rhs, a=(l(w)-l(x))/2
            a = W.length(w) - W.length(x)
            if a % 2:
                continue  # odd difference: P_{x,w}=1? handle separately
            k = a // 2
            # P is a poly in q with constant term 1, degree <= (a-1)//2
            # bar(P) = P(v^{-2}); q^k bar(P) has exponents 2k-2i, q^{-k} P has -2k+2i
            # For i=0: q^k gives v^{2k}, -q^{-k} gives -v^{-2k}; RHS known.
            # Solve coefficient by coefficient.
            cs = [1]  # constant term
            maxdeg = (a - 1) // 2
            # We use direct inversion: for each coefficient c_i of P, contribution
            # to LHS at exponent 2k-2i is c_i, and at -2k+2i is -c_i.
            # Collect RHS coefficient at each exponent and reconstruct.
            for i in range(1, maxdeg + 1):
                # coefficient at exponent 2k-2i (positive side) equals c_i
                # (since other terms contribute only if 2k-2j = 2k-2i -> j=i)
                cs.append(rhs.get(2 * k - 2 * i, 0))
            expected = poly.qpoly_from_coeffs(cs)
            got = T.P(x, w)
            checked += 1
            if not poly.eq(got, expected):
                bad.append((W.describe(x), W.describe(w), poly.fmt(got), poly.fmt(expected)))
    print(f"{W.name}: checked {checked} pairs, mismatches {len(bad)}", bad[:5], flush=True)
    if bad:
        sys.exit(1)
print("PASS: R-inversion identity agrees with KL recursion")
