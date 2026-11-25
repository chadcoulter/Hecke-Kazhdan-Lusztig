"""Kazhdan-Lusztig and R polynomials for a finite Coxeter group.

KL recursion (Kazhdan-Lusztig 1979, (2.2.c)): for y<=w and s in the LEFT
descents of w (sw<w):

  if sy > y:  P_{y,w} = P_{sy,w}
  else:       P_{y,w} = P_{sy,sw} + q P_{y,sw}
              - sum_{z: y<=z<sw, sz<z} mu(z,sw) q^{(l(w)-l(z))/2} P_{y,z}

mu(x,y) = coefficient of q^{(l(y)-l(x)-1)/2} in P_{x,y}, else 0.
For y not<= w: P=0. P_{w,w}=1.

R recursion (KL79 / Wikipedia): for y fixed, s with sy<y:
  R_{y,y}=1; R_{x,y}=R_{sx,sy} if sx<x else (q-1)R_{sx,y}+q R_{sx,sy}.
Identity: T_{y^{-1}}^{-1} = sum_x bar(R_{x,y}) q^{-l(x)} T_x.
"""
from __future__ import annotations

from . import poly
from .poly import Poly


class KLTable:
    def __init__(self, W):
        self.W = W
        self._P: dict[tuple, Poly] = {}
        self._R: dict[tuple, Poly] = {}
        self._mu: dict[tuple, int] = {}

    def R(self, x, y) -> Poly:
        W = self.W
        x, y = W._check(x), W._check(y)
        key = (x, y)
        if key in self._R:
            return self._R[key]
        if not W.leq(x, y):
            self._R[key] = poly.zero()
            return self._R[key]
        if x == y:
            self._R[key] = poly.one()
            return self._R[key]
        # pick s in left descent of y with sy<y, and use the KL79 recursion
        ds = W.left_descents(y)
        s = ds[0]
        sy = W.left_mult(s, y)
        sx = W.left_mult(s, x)
        qp = poly.q()
        qm1 = poly.sub(qp, poly.one())
        if W.length(sx) < W.length(x):
            res = self.R(sx, sy)
        else:
            res = poly.add(poly.mul(qm1, self.R(x, sy)),
                           poly.mul(qp, self.R(sx, sy)))
        self._R[key] = res
        return res

    def P(self, y, w) -> Poly:
        W = self.W
        y, w = W._check(y), W._check(w)
        key = (y, w)
        if key in self._P:
            return self._P[key]
        if not W.leq(y, w):
            self._P[key] = poly.zero()
            return self._P[key]
        if y == w:
            self._P[key] = poly.one()
            return self._P[key]
        ds = W.left_descents(w)
        if not ds:
            # only identity has no descent; y != identity handled by leq check
            self._P[key] = poly.zero()
            return self._P[key]
        s = ds[0]
        sw = W.left_mult(s, w)
        sy = W.left_mult(s, y)
        qp = poly.q()
        result: Poly = poly.zero()
        if W.length(sy) > W.length(y):
            result = self.P(sy, w)
        else:
            first = self.P(sy, sw)
            second = poly.mul(qp, self.P(y, sw))
            result = poly.add(first, second)
            for z in W.strictly_below_with_descents(y, sw, s):
                m = self.mu(z, sw)
                if m == 0:
                    continue
                term = poly.mul(
                    poly.scale(self.P(y, z), m),
                    poly.q_power((W.length(w) - W.length(z)) // 2),
                )
                result = poly.sub(result, term)
        self._P[key] = result
        return result

    def mu(self, x, y) -> int:
        W = self.W
        x, y = W._check(x), W._check(y)
        key = (x, y)
        if key in self._mu:
            return self._mu[key]
        d = W.length(y) - W.length(x) - 1
        if d < 0 or d % 2 != 0 or not W.leq(x, y):
            self._mu[key] = 0
            return 0
        # get P_{x,y}. Guard against recursion cycle (recursion only calls
        # mu(z,sw) with z<sw, i.e. strictly shorter second argument).
        self._mu[key] = poly.coeff_q(self.P(x, y), d // 2)
        return self._mu[key]

    def all_polynomials(self) -> dict[tuple, Poly]:
        for w in self.W.elements:
            for y in self.W.elements:
                self.P(y, w)
        return dict(self._P)
