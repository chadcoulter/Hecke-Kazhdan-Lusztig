#!/usr/bin/env python3
import sys
sys.path.insert(0, "/home/sudoshellz/Documents/hecke-kl")
from hecke_kl.groups import make
from hecke_kl.hecke import HeckeAlgebra
from hecke_kl import poly
W = make("A", 2); H = HeckeAlgebra(W)
print("inverse_Ts", {W.describe(k): poly.fmt(v, "v") for k,v in H.inverse_Ts(0).items()})
a = H.inverse_Ts(0)
print("mul(a,a):")
try:
    r = H.mul(a, a)
    print({W.describe(k): poly.fmt(v, "v") for k,v in r.items()})
except Exception as e:
    print("ERR", type(e).__name__, e)
