#!/usr/bin/env python3
import sys
sys.path.insert(0, __import__("pathlib").Path(__file__).resolve().parent)
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl import poly
W = make("A", 2); T = KLTable(W)
s0, sts = W.S[0], W.longest
sw = W.left_mult(0, sts)
print("sw for s=0:", W.describe(sw), "len", W.length(sw))
print("left descents of sw", [W.describe(W.left_mult(i,sw)) for i in W.left_descents(sw)])
zs = W.strictly_below_with_descents(s0, sw)
print("zs", [W.describe(z) for z in zs])
for z in zs:
    print("  z", W.describe(z), "mu", T.mu(z,sw), "P[s0,z]", poly.fmt(T.P(s0,z)))
print("fresh P[s0,sts]", poly.fmt(KLTable(W).P(s0, sts)))
