
import sys
sys.path.insert(0, "/home/sudoshellz/Documents/hecke-kl")
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl.hecke import HeckeAlgebra
from hecke_kl import poly
import traceback

W = make("A", 4)   # rank-4 A_4 = S5
T = KLTable(W); H = HeckeAlgebra(W, T)
print("W", W.name, "order", W.order, flush=True)
for w in W.elements[:12]:
    try:
        H.solve_kl_bar(w)
        print(W.describe(w), W.length(w), "OK", flush=True)
    except Exception as e:
        print(W.describe(w), W.length(w), "FAIL", type(e).__name__, repr(str(e))[:160], flush=True)
        traceback.print_exc()
        break

# also check A3 (rank 3 = S4) quickly
W3 = make("A", 3); H3 = HeckeAlgebra(W3, KLTable(W3))
fails = []
for w in W3.elements:
    try:
        H3.solve_kl_bar(w)
    except Exception as e:
        fails.append((W3.describe(w), W3.length(w), type(e).__name__, str(e)[:100]))
print("A3(rank3=S4) fails:", len(fails), fails[:5], flush=True)
