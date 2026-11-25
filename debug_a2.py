#!/usr/bin/env python3
"""Isolate the first divergence: check lengths, reflections, and simple cases."""
import sys
sys.path.insert(0, "/home/sudoshellz/Documents/hecke-kl")
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl.hecke import HeckeAlgebra
from hecke_kl import poly

W = make("A", 2)
print("order", W.order, "longest", W.describe(W.longest), "len", W.length(W.longest))
print("reflections", [W.describe(t) for t in W.reflections])
for w in W.elements:
    print("w", W.describe(w), "len", W.length(w), "word", W.word(w),
          "Ldesc", [i for i in range(W.n) if W.length(W.mul(W.S[i], w)) < W.length(w)])

s = W.S[0]
print("s*s == e?", W.mul(s, s) == W.identity)
print("word of s*s", W.word(W.mul(s, s)))
T = KLTable(W)
print("P[e,e]", poly.fmt(T.P(W.identity, W.identity)))
y = W.from_word((0,)); w = W.from_word((0,1,0))
print("P[s,sts]", poly.fmt(T.P(y, w)))
print("R[s,sts]", poly.fmt(T.R(y, w)))
H = HeckeAlgebra(W, T)
lhs = H.mul(H.T(s), H.T(s))
print("T_s*T_s", {W.describe(k): poly.fmt(v, "v") for k,v in lhs.items()})
rhs = H.add({W.identity: poly.q()}, H.scale(H.T(s), poly.sub(poly.q(), poly.one())))
print("q+(q-1)T_s", {W.describe(k): poly.fmt(v, "v") for k,v in rhs.items()})
