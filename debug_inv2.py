#!/usr/bin/env python3
import sys
sys.path.insert(0, "/home/sudoshellz/Documents/hecke-kl")
from hecke_kl.groups import make
from hecke_kl.hecke import HeckeAlgebra
from hecke_kl import poly
W = make("A", 3); H = HeckeAlgebra(W)
s = W.S[0]
inv = H.inverse_Ts(0)
prod = H.mul(H.T(s), inv)
print("T_s*T_s^-1:", {W.describe(k): poly.fmt(v) for k,v in prod.items()})
# manual: T_s*(q^-1 T_s + q^-1 -1) = q^-1(q+(q-1)T_s) + (q^-1-1)T_s
# = 1 + (q^-1-1+q^-1(q-1))T_s = 1 + (q^-1-1+1-q^-1)T_s = 1
# check mul_left_s(0, {s: qinv}) manually
qinv = poly.shift(poly.one(), -2)
one = {W.S[0]: qinv}
print("mul_left_s on T_s*qinv:", {W.describe(k): poly.fmt(v) for k,v in H.mul_left_s(0, one).items()})
