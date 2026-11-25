
import sys, traceback
sys.path.insert(0, "/home/sudoshellz/Documents/hecke-kl")
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl.hecke import HeckeAlgebra
W=make("A",3); H=HeckeAlgebra(W,KLTable(W))
print("order",W.order,flush=True)
for w in W.elements:
    try:
        r=H.solve_kl_bar(w)
        print("OK",W.describe(w),W.length(w),len(r),flush=True)
    except Exception as e:
        print("FAIL",W.describe(w),W.length(w),type(e).__name__,repr(str(e))[:100],flush=True)
