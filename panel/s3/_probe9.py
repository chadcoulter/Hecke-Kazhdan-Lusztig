
import sys, sympy as sp
sys.path.insert(0, "/home/sudoshellz/Documents/hecke-kl")
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl.hecke import HeckeAlgebra

x = sp.Symbol("x_probe_unique")
print("BEFORE sp.solve([],[x],dict=True):", sp.solve([], [x], dict=True), flush=True)
W4=make("A",4); H4=HeckeAlgebra(W4,KLTable(W4))
try:
    H4.solve_kl_bar(W4.identity)
    print("A4 e: OK (unexpected)", flush=True)
except Exception as e:
    print("A4 e FAIL:", type(e).__name__, repr(str(e))[:120], flush=True)
print("AFTER  sp.solve([],[x],dict=True):", sp.solve([], [x], dict=True), flush=True)
W3=make("A",3); H3=HeckeAlgebra(W3,KLTable(W3))
fails=[]
for w in W3.elements:
    try:
        H3.solve_kl_bar(w)
    except Exception as e:
        fails.append((W3.describe(w),W3.length(w),type(e).__name__,str(e)[:80]))
print("A3 fails after A4 failure:",len(fails),fails[:6], flush=True)
