
import sys, time
sys.path.insert(0, "/home/sudoshellz/Documents/hecke-kl")
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl.hecke import HeckeAlgebra
from hecke_kl import poly

W=make("D",4); H=HeckeAlgebra(W,KLTable(W))
fails=[]; ck=0; bad_match=[]; t0=time.time()
for w in W.elements:
    if W.length(w)>4: continue
    try:
        rec=H.solve_kl_bar(w); ck+=1
    except Exception as e:
        fails.append((W.describe(w),W.length(w),type(e).__name__,str(e)[:80])); continue
    for y in W.elements:
        if not W.leq(y,w) or y==w: continue
        d=W.length(w)-W.length(y)
        cs=[1]+[rec.get((y,i),poly.zero()).get(0,0) for i in range(1,(d-1)//2+1)]
        exp=poly.qpoly_from_coeffs(cs); got=H.kl.P(y,w)
        if not poly.eq(got,exp): bad_match.append((W.describe(y),W.describe(w),poly.fmt(got),poly.fmt(exp)))
print("D4 solve len<=4 ok",ck,"fails",len(fails),fails[:8],"mismatch",len(bad_match),bad_match[:6],"t=%.1f"%(time.time()-t0),flush=True)
print("DONE",flush=True)
