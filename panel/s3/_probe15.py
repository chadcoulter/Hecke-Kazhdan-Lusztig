
import sys, time
sys.path.insert(0, "/home/sudoshellz/Documents/hecke-kl")
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl.hecke import HeckeAlgebra
from hecke_kl import poly

for typ,n in [("D",4),("B",4)]:
    W=make(typ,n); H=HeckeAlgebra(W,KLTable(W))
    fails=[]; ck=0; t0=time.time()
    for w in W.elements:
        try:
            rec=H.solve_kl_bar(w); ck+=1
        except Exception as e:
            fails.append((W.describe(w),W.length(w),type(e).__name__,str(e)[:90]))
            continue
    print(f"[solve] {W.name} order={W.order} ok={ck} fails={len(fails)} {fails[:8]} t={time.time()-t0:.1f}",flush=True)
print("DONE",flush=True)
