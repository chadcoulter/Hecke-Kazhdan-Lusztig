"""Hecke algebra H(W,S) over Z[v^{±1}], q = v^2.

Basis T_w. Relations:
  T_e = 1
  T_s T_s = q + (q-1) T_s
  T_s T_w = T_{sw}          if l(sw) > l(w)
  T_s T_w = q T_{sw} + (q-1) T_w   otherwise
Bar involution: bar(v)=v^{-1}, bar(T_w) = T_{w^{-1}}^{-1}.

C'_w = v^{-l(w)} sum_{y<=w} P_{y,w}(v^2) T_y is bar-invariant.
This gives an INDEPENDENT route to P: set up bar(C'_w)=C'_w as a linear
system over Z[v^{±1}] and solve by matching coefficients of T_y.
"""
from __future__ import annotations

from collections import defaultdict
from . import poly
from .poly import Poly


class HeckeAlgebra:
    def __init__(self, W, kl=None):
        self.W = W
        self.kl = kl

    def T_scalar(self, p: Poly) -> dict:
        """The Hecke element p * T_e."""
        return {self.W.identity: dict(p)} if p else {}

    def T(self, w) -> dict:
        """Basis element T_w as a sparse map group_element -> Poly."""
        return {self.W._check(w): poly.one()}

    def add(self, a, b):
        r = defaultdict(Poly, {k: dict(v) for k, v in a.items()})
        for k, v in b.items():
            r[k] = poly.add(r.get(k, {}), v)
        return {k: v for k, v in r.items() if v}

    def scale(self, a, p):
        return {k: poly.mul(v, p) for k, v in a.items() if poly.mul(v, p)}

    def neg(self, a):
        return {k: poly.neg(v) for k, v in a.items()}

    def mul_left_s(self, i: int, a: dict) -> dict:
        """Multiply T_s_i * a on the left."""
        W = self.W
        s = W.S[i]
        qp = poly.q()
        qm1 = poly.sub(qp, poly.one())
        out: dict = {}
        for w, c in a.items():
            sw = W.mul(s, w)
            if W.length(sw) > W.length(w):
                term = {sw: c}
            else:
                term = {
                    sw: poly.mul(qp, c),
                    w: poly.mul(qm1, c),
                }
            out = self.add(out, term)
        return out

    def mul(self, a: dict, b: dict) -> dict:
        """Full product: T_u * T_v, reduced generator by generator on the right."""
        W = self.W
        out: dict = {}
        for u, cu in a.items():
            for v, cv in b.items():
                # Acc = T_v: T_{i1}...T_{ik}, so right-multiply in forward order.
                acc: dict = {W.identity: dict(cv)}
                for i in W.word(v):
                    acc = self._mul_right_s(i, acc)
                # T_u on the left: T_{j1}...T_{jl}·acc, apply reversed.
                piece: dict = {}
                for z, cz in acc.items():
                    pz: dict = {z: dict(cz)}
                    for j in reversed(W.word(u)):
                        pz = self.mul_left_s(j, pz)
                    piece = self.add(piece, pz)
                for z, pz in piece.items():
                    out = self.add(out, {z: poly.mul(pz, cu)})
        return out

    def _mul_right_s(self, i: int, a: dict) -> dict:
        W = self.W
        s = W.S[i]
        qp = poly.q()
        qm1 = poly.sub(qp, poly.one())
        out: dict = {}
        for w, c in a.items():
            ws = W.mul(w, s)
            if W.length(ws) > W.length(w):
                out = self.add(out, {ws: c})
            else:
                out = self.add(out, {ws: poly.mul(qp, c), w: poly.mul(qm1, c)})
        return out

    def inverse_Ts(self, i: int) -> dict:
        """T_s^{-1} = q^{-1} T_s + q^{-1} - 1."""
        s = self.W.S[i]
        qinv = {e: c for e, c in poly.shift(poly.one(), -2).items()}  # v^{-2}
        qminv1 = poly.sub(qinv, poly.one())
        return {s: qinv, self.W.identity: qminv1}

    def inverse_Tw(self, w) -> dict:
        w = self.W._check(w)
        acc = {self.W.identity: poly.one()}
        for i in self.W.word(w):  # inverse of T_{i1}...T_{ik} is T_{ik}^{-1}...T_{i1}^{-1}
            acc = self.mul(self.inverse_Ts(i), acc)
        return acc

    def bar_T(self, w) -> dict:
        return self.inverse_Tw(self.W.inverse(w))

    def bar(self, a: dict) -> dict:
        out: dict = {}
        for w, c in a.items():
            out = self.add(out, self.mul(self.bar_T(w), self.T_scalar(poly.bar(c))))
        return out

    def Cprime(self, w) -> dict:
        """C'_w = v^{-l(w)} sum_{y<=w} P_{y,w}(v^2) T_y."""
        if self.kl is None:
            raise ValueError("HeckeAlgebra needs a KLTable for Cprime")
        W = self.W
        w = W._check(w)
        out: dict = {}
        shift = -W.length(w)
        for y in W.elements:
            if W.leq(y, w):
                p = self.kl.P(y, w)
                if p:
                    out[y] = poly.shift(p, shift)
        return out

    def equal(self, a, b) -> bool:
        keys = set(a) | set(b)
        return all(poly.eq(a.get(k, {}), b.get(k, {})) for k in keys)

    # ------------------------------------------------- independent KL via bar-invariance
    def solve_kl_bar(self, w) -> dict[tuple, Poly]:
        """Recover P_{y,w} by solving bar(C'_w)=C'_w, WITHOUT the KL recursion.

        Unknowns c_{y,i} (P_{y,w} = sum_i c_{y,i} q^i) are found by equating
        coefficients of T_z and of each v-exponent in bar(C') = C'. Every
        equation is LINEAR in the unknowns, so we build a sparse linear system
        indexed by (T_z, exponent) and solve it with sympy's linsolve. No
        symbolic rational functions are used (the v^{-1} entries from the bar
        involution are handled as integer exponent shifts), which avoids the
        PolynomialError that a sp.Poly-on-1/v approach hits.
        """
        import sympy as sp

        W = self.W
        w = W._check(w)
        unknowns = {}
        for y in W.elements:
            if not W.leq(y, w):
                continue
            d = W.length(w) - W.length(y)
            maxdeg = 0 if d == 0 else (d - 1) // 2
            for i in range(maxdeg + 1):
                if y == w and i == 0:
                    continue  # P_{w,w}=1 is known, not an unknown
                unknowns[(y, i)] = sp.Symbol(f"c_{W.describe(y)}_{i}")
        syms = list(unknowns.values())

        # A linear form is dict[Symbol|None, int]: None is the constant slot.
        def form_const(k):
            return {None: k} if k else {}

        def form_add(a, b):
            r = dict(a)
            for k, v in b.items():
                r[k] = r.get(k, 0) + v
                if r[k] == 0:
                    del r[k]
            return r

        def form_scale(a, k):
            return {s: c * k for s, c in a.items() if c * k}

        # C'_w as map T_z -> {v-exponent -> linear form}
        C: dict = {}
        for (y, i), sym in unknowns.items():
            e = 2 * i - W.length(w)
            d = C.setdefault(y, {})
            d[e] = form_add(d.get(e, {}), {sym: 1})
        # P_{w,w} = 1 contributes v^{-l(w)} to T_w
        cw = C.setdefault(w, {})
        e0 = -W.length(w)
        cw[e0] = form_add(cw.get(e0, {}), form_const(1))

        # bar(C'): coefficient v^e -> v^{-e}, then bar(T_z).
        barC: dict = {}
        for z, coeffs in C.items():
            for e, form in coeffs.items():
                for t, cpoly in self.bar_T(z).items():
                    for e2, cint in cpoly.items():
                        d = barC.setdefault(t, {})
                        d[-e + e2] = form_add(d.get(-e + e2, {}), form_scale(form, cint))

        # Equations: for every t and exponent, barC[t][k] - C[t][k] = 0.
        rows = []
        for t in set(C) | set(barC):
            cmap, bmap = C.get(t, {}), barC.get(t, {})
            for k in set(cmap) | set(bmap):
                diff = form_add(bmap.get(k, {}), form_scale(cmap.get(k, {}), -1))
                if diff:
                    rows.append(diff)
        if not rows:
            # w = e: only unknown was already fixed; return empty (identity P=1).
            return {}

        # Build coefficient matrix over the unknowns and solve.
        A, b = [], []
        for row in rows:
            A.append([row.get(s, 0) for s in syms])
            b.append(-row.get(None, 0))
        sol = sp.linsolve((sp.Matrix(A), sp.Matrix(b)), syms)
        solution = list(sol)[0] if sol else None
        if solution is None:
            raise RuntimeError(f"no bar-invariant solution for w={W.describe(w)}")
        out = {}
        for (y, i), sym in unknowns.items():
            val = solution[syms.index(sym)]
            if not val.is_Integer:
                raise RuntimeError(f"non-integer KL coefficient at y={W.describe(y)}, i={i}: {val}")
            if int(val) != 0:
                out[(y, i)] = poly.qpoly_from_coeffs([int(val)])
        return out
