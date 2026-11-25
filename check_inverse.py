#!/usr/bin/env python3
import sys
sys.path.insert(0, __import__("pathlib").Path(__file__).resolve().parent)
from hecke_kl.groups import make
from hecke_kl.hecke import HeckeAlgebra
from hecke_kl import poly
W = make("A", 3); H = HeckeAlgebra(W)
for i in range(W.n):
    s = W.S[i]
    got = H.inverse_Tw(s)
    exp = H.inverse_Ts(i)
    print(f"s{i} inv matches:", H.equal(got, exp))
# T_w T_w^{-1} = 1 for all w
bad = []
for w in W.elements:
    prod = H.mul(H.T(w), H.inverse_Tw(w))
    if not H.equal(prod, H.T(W.identity)):
        bad.append(W.describe(w))
print("T_w T_w^-1 = 1 failures:", bad[:5], "count", len(bad))
