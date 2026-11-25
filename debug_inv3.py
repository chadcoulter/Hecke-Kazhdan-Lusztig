#!/usr/bin/env python3
import sys
sys.path.insert(0, __import__("pathlib").Path(__file__).resolve().parent)
from hecke_kl.groups import make
from hecke_kl.hecke import HeckeAlgebra
from hecke_kl import poly
W = make("A", 3); H = HeckeAlgebra(W)
w = W.from_word((0,1))
inv = H.inverse_Tw(w)
print("word", W.word(w))
print("inv", {W.describe(k): poly.fmt(v) for k,v in inv.items()})
prod = H.mul(H.T(w), inv)
print("prod", {W.describe(k): poly.fmt(v) for k,v in prod.items()})
