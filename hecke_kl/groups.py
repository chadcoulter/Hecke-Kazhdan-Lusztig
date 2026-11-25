"""Cartan matrices for supported finite Coxeter types.

A_ij = 2(α_i,α_j)/(α_i,α_i); A_ii=2, A_ij<=0 for i!=j, with A_ij A_ji in
{0,1,2,3}. The integer Tits representation s_i(α_j)=α_j - A_ij α_i is faithful
for these finite groups.
"""
from __future__ import annotations

from .coxeter import CoxeterGroup


def cartan_a(n: int) -> list[list[int]]:
    if n < 1:
        raise ValueError("A_n requires n>=1")
    m = [[0] * n for _ in range(n)]
    for i in range(n):
        m[i][i] = 2
        if i > 0:
            m[i][i - 1] = -1
        if i < n - 1:
            m[i][i + 1] = -1
    return m


def cartan_b(n: int) -> list[list[int]]:
    """Type B_n/C_n (n>=2): 1 = 2 = 3 = ... = n-1 => n with A_{n-1,n}=-2."""
    if n < 2:
        raise ValueError("B_n requires n>=2")
    m = cartan_a(n)
    m[n - 2][n - 1] = -2
    return m


def cartan_d(n: int) -> list[list[int]]:
    """Type D_n (n>=4), diagram n-1 and n both attached to n-2."""
    if n < 4:
        raise ValueError("D_n requires n>=4")
    m = [[0] * n for _ in range(n)]
    for i in range(n):
        m[i][i] = 2
    # path 0-1-...-(n-3), then (n-3) branches to (n-2) and (n-1)
    for i in range(n - 3):
        m[i][i + 1] = m[i + 1][i] = -1
    m[n - 3][n - 2] = m[n - 2][n - 3] = -1
    m[n - 3][n - 1] = m[n - 1][n - 3] = -1
    return m


def cartan_g2() -> list[list[int]]:
    return [[2, -1], [-3, 2]]  # A_01=-1, A_10=-3; braid order 6


def cartan_i2(m: int) -> list[list[int]]:
    """Dihedral I_2(m), m>=3: A_01=-(2cos(pi/m)) approximated only for
    crystallographic m in {3,4,6}; for non-crystallographic m the Tits integer
    representation does not apply. We support m=3 (A2), m=4 (B2), m=6 (G2)."""
    table = {3: [[2, -1], [-1, 2]], 4: [[2, -2], [-1, 2]], 6: [[2, -3], [-1, 2]]}
    if m not in table:
        raise ValueError(f"I2(m) integer Tits rep supports only m=3,4,6; got {m}")
    return table[m]


def make(name: str, n: int | None = None) -> CoxeterGroup:
    name = name.upper()
    if name == "A":
        mat = cartan_a(n)
    elif name == "A2":
        mat = cartan_a(2)
    elif name in ("B", "C"):
        mat = cartan_b(n)
    elif name == "D":
        mat = cartan_d(n)
    elif name == "G2":
        mat = cartan_i2(6)
    elif name == "I2":
        mat = cartan_i2(n)
    else:
        raise ValueError(f"unsupported type: {name}")
    label = name if n is None else f"{name}{n}"
    return CoxeterGroup(label, mat)
