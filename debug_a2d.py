#!/usr/bin/env python3
import sys
sys.path.insert(0, __import__("pathlib").Path(__file__).resolve().parent)
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl import poly
W = make("A", 2); T = KLTable(W)
e = W.identity
sts = W.longest
s = W.left_descents(sts)[0]
sw = W.left_mult(s, sts)
print("chosen s", s, "sw", W.describe(sw))
print("left descents of sw", [W.describe(W.left_mult(i,sw)) for i in W.left_descents(sw)])
zs = W.strictly_below_with_descents(e, sw)
print("zs", [W.describe(z) for z in zs])
for z in zs:
    print("  z", W.describe(z), "mu(z,sw)", T.mu(z,sw), "P[e,z]", poly.fmt(T.P(e,z)))
