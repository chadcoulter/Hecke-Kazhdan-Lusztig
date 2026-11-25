#!/usr/bin/env python3
import sys
sys.path.insert(0, "/home/sudoshellz/Documents/hecke-kl")
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl import poly
W = make("A", 2); T = KLTable(W)
e, s0, s1, st = (W.identity, W.S[0], W.S[1], W.from_word((0,1)))
sts = W.longest
print("left descents of sts:", [W.describe(W.left_mult(i, sts)) for i in W.left_descents(sts)])
s0e = W.left_mult(0, e); print("s0*e", W.describe(s0e), W.length(s0e))
s0sts = W.left_mult(0, sts); print("s0*sts", W.describe(s0sts), W.length(s0sts))
print("P[s0, s0sts] =", poly.fmt(T.P(s0, s0sts)), "expected 1 (len diff 1)")
print("P[s0, sts] =", poly.fmt(T.P(s0, sts)))
print("mu(s0, sts) =", T.mu(s0, sts))
# R check
print("R[e,sts]", poly.fmt(T.R(e, sts)), "expected (q-1)^2? or q^2-2q+1")
