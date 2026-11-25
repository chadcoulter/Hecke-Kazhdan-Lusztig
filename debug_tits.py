#!/usr/bin/env python3
import sys
sys.path.insert(0, "/home/sudoshellz/Documents/hecke-kl")
from hecke_kl.groups import make
W = make("A", 2)
n = W.n
# Expected: order 6, lengths 0,1,1,2,2,3
print("n", n)
print("cartan", W.cartan)
for i in range(n):
    for j in range(n):
        if i != j:
            prod = W.mul(W.S[i], W.S[j])
            print(f"s{i}s{j} order via word length", W.length(prod))
print("mat S0", W.S[0])
print("mat S1", W.S[1])
# verify s_i^2 = identity
for i in range(n):
    print(f"s{i}^2=e", W.mul(W.S[i], W.S[i]) == W.identity)
# braid (s0 s1)^3 = e
w = W.identity
for _ in range(3):
    w = W.mul(w, W.mul(W.S[0], W.S[1]))
print("(s0s1)^3=e", w == W.identity, W.describe(w))
# Compare A2 algebra: T_s^2 via mul_left_s on identity should use sw>s path
# sw = s*s = e shorter -> quadratic branch. Check length(e)=0 < length(s)=1 yes
H = None
from hecke_kl.hecke import HeckeAlgebra
from hecke_kl import poly
H = HeckeAlgebra(W)
print("mul_left_s on T_s", {W.describe(k): poly.fmt(v, "v") for k,v in H.mul_left_s(0, H.T(W.S[0])).items()})
