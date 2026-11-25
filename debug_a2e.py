#!/usr/bin/env python3
import sys
sys.path.insert(0, __import__("pathlib").Path(__file__).resolve().parent)
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl import poly
W = make("A", 2); T = KLTable(W)
e, sts = W.identity, W.longest
print("fresh P[e,sts] =", poly.fmt(T.P(e, sts)))
T2 = KLTable(W)
print("fresh P[s0,sts] =", poly.fmt(T2.P(W.S[0], sts)))
T3 = KLTable(W)
print("fresh P[e,ts] =", poly.fmt(T3.P(e, W.from_word((1,0)))))
