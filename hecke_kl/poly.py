"""Polynomials in v (with q = v**2) as sparse dicts: exponent -> integer coefficient.

All coefficients are exact Python ints; no floating point anywhere in the core.
"""
from __future__ import annotations
from collections import defaultdict

Poly = dict[int, int]  # v-exponent -> integer coefficient


def clean(p: Poly) -> Poly:
    return {e: c for e, c in p.items() if c != 0}


def zero() -> Poly:
    return {}


def one() -> Poly:
    return {0: 1}


def const(c: int) -> Poly:
    return {0: c} if c else {}


def var() -> Poly:
    return {1: 1}


def q() -> Poly:
    return {2: 1}


def add(a: Poly, b: Poly) -> Poly:
    r = defaultdict(int, a)
    for e, c in b.items():
        r[e] += c
    return clean(r)


def neg(a: Poly) -> Poly:
    return {e: -c for e, c in a.items()}


def sub(a: Poly, b: Poly) -> Poly:
    return add(a, neg(b))


def mul(a: Poly, b: Poly) -> Poly:
    if not a or not b:
        return {}
    r: dict[int, int] = defaultdict(int)
    for e1, c1 in a.items():
        for e2, c2 in b.items():
            r[e1 + e2] += c1 * c2
    return clean(r)


def scale(a: Poly, k: int) -> Poly:
    return {e: c * k for e, c in a.items() if c * k}


def shift(a: Poly, k: int) -> Poly:
    """Multiply by v**k (k may be negative)."""
    return {e + k: c for e, c in a.items()} if a else {}


def bar(a: Poly) -> Poly:
    """Bar involution: v -> v^{-1}."""
    return {(-e): c for e, c in a.items()}


def degree(a: Poly) -> int:
    return max(a) if a else -1


def leading_coeff(a: Poly) -> int:
    return a[degree(a)] if a else 0


def is_zero(a: Poly) -> bool:
    return not a


def eq(a: Poly, b: Poly) -> bool:
    return clean(a) == clean(b)


def qpoly_from_coeffs(coeffs: list[int]) -> Poly:
    """coeffs[i] is coefficient of q**i; returns polynomial in v (exponents 2i)."""
    return {2 * i: c for i, c in enumerate(coeffs) if c != 0}


def q_power(k: int) -> Poly:
    """q**k = v**(2k) as a polynomial; k>=0 else error."""
    if k < 0:
        raise ValueError("q_power requires k>=0")
    return {2 * k: 1}


def coeff_q(a: Poly, i: int) -> int:
    return a.get(2 * i, 0)


def qdegree(a: Poly) -> int:
    return degree(a) // 2 if a and all(e % 2 == 0 for e in a) else -1


def as_q_list(a: Poly) -> list[int]:
    if not a:
        return []
    if any(e % 2 for e in a):
        raise ValueError("not a polynomial in q (odd v-exponent present)")
    return [coeff_q(a, i) for i in range(qdegree(a) + 1)]


def fmt(a: Poly, variable: str = "q") -> str:
    """Human-readable string. If a is a polynomial in nonnegative q, print in q."""
    if not a:
        return "0"
    if all(e >= 0 and e % 2 == 0 for e in a):
        cs = as_q_list(a)
        parts = []
        for i, c in enumerate(cs):
            if c == 0:
                continue
            if i == 0:
                term = str(c)
            elif c == 1:
                term = variable
            elif c == -1:
                term = "-" + variable
            else:
                term = f"{c}{variable}"
            if i > 0:
                term += f"^{i}" if i > 1 else ""
            parts.append(term)
        if parts:
            s = parts[0]
            for t in parts[1:]:
                s += (" + " if not t.startswith("-") else " - ") + (t[1:] if t.startswith("-") else t)
            return s
    parts = []
    for e, c in sorted(a.items(), reverse=True):
        if c == 0:
            continue
        if e == 0:
            vterm = "1"
        elif e == 1:
            vterm = "v"
        elif e == -1:
            vterm = "v^-1"
        else:
            vterm = f"v^{e}"
        parts.append(f"{c}*{vterm}")
    return " + ".join(parts)
