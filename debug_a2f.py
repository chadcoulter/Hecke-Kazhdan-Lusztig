#!/usr/bin/env python3
import sys
sys.path.insert(0, __import__("pathlib").Path(__file__).resolve().parent)
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl import poly
W = make("A", 2); T = KLTable(W)
s0, ba = W.S[0], W.from_word((1,0))
print("P[s0,ba]", poly.fmt(T.P(s0, ba)), "lens", W.length(s0), W.length(ba))
print("leq(s0,ba)", W.leq(s0, ba))
print("s0 descents", [W.describe(W.left_mult(i,s0)) for i in W.left_descents(s0)])
# P_{s0,ba}: left descents of ba
print("ba descents", [W.describe(W.left_mult(i,ba)) for i in W.left_descents(ba)])
s = W.left_descents(ba)[0]; sw = W.left_mult(s,ba); sy = W.left_mult(s,s0)
print("s",s,"sw",W.describe(sw),"sy",W.describe(sy),"len sy",W.length(sy),"len s0",W.length(s0))
print("zs", [W.describe(z) for z in W.strictly_below_with_descents(s0, sw)])
