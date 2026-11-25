#!/usr/bin/env python3
import sys
sys.path.insert(0, __import__("pathlib").Path(__file__).resolve().parent)
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl.hecke import HeckeAlgebra
from hecke_kl import poly

W = make("A", 2)
T = KLTable(W)
H = HeckeAlgebra(W, T)
# check mul: T_s * T_s should equal q + (q-1) T_s
s = W.S[0]
prod = H.mul(H.T(s), H.T(s))
print("prod", {W.describe(k): poly.fmt(v, "v") for k,v in prod.items()})
print("expected", {W.describe(k): poly.fmt(v, "v") for k,v in {W.identity: poly.q(), s: poly.sub(poly.q(), poly.one())}.items()})
# is mul_left_s applied correctly in mul? mul() loops over right basis b's word
# For b=T_s, word=(0,), acc={e:c}, calls _mul_right_s(0, acc)
print("right_s on T_e", {W.describe(k): poly.fmt(v, "v") for k,v in H._mul_right_s(0, H.T(W.identity)).items()})

# P[e,sts] should be 1, P[s,sts] should be 1. Check both
for yw in [(), (0,), (1,), (0,1), (1,0), (0,1,0)]:
    y = W.from_word(yw)
    print("P", yw, "sts =", poly.fmt(T.P(y, W.longest)))
