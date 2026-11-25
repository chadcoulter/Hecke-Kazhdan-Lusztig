"""Hecke–KL simulator: finite Coxeter groups, Hecke algebras, KL polynomials."""
from .poly import Poly, fmt, q, q_power, add, mul, bar, as_q_list
from .coxeter import CoxeterGroup
from .groups import make, cartan_a, cartan_b, cartan_d, cartan_i2
from .kl import KLTable
from .hecke import HeckeAlgebra
from . import types

__version__ = "0.1.0"
__all__ = [
    "Poly", "fmt", "q", "q_power", "add", "mul", "bar", "as_q_list",
    "CoxeterGroup", "make", "cartan_a", "cartan_b", "cartan_d", "cartan_i2",
    "KLTable", "HeckeAlgebra", "types", "__version__",
]
