
import sys, time
sys.path.insert(0, __import__("pathlib").Path(__file__).resolve().parent)
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl.hecke import HeckeAlgebra
from hecke_kl import poly

for typ,n,wl in [("A",3,None),("A",4,None),("B",3,None)]:
    W=make(typ,n); H=HeckeAlgebra(W,KLTable(W))
    print("=== W",W.name,"order",W.order,flush=True)
    for w in W.elements:
        if W.length(w)>1: continue
        t=time.time()
        try:
            r=H.solve_kl_bar(w)
            ok = True
            # compare
            bad=[]
            for y in W.elements:
                if not W.leq(y,w): continue
                d=W.length(w)-W.length(y); cs=[1]
                for i in range(1,(d-1)//2+1): cs.append(r.get((y,i),poly.zero()).get(0,0))
                exp=poly.qpoly_from_coeffs(cs); got=KLTable(W).P(y,w)
                if not poly.eq(got,exp): bad.append((W.describe(y),poly.fmt(got),poly.fmt(exp)))
            print(" ",W.describe(w),W.length(w),"OK bad=",bad,"t=%.1f"%(time.time()-t),flush=True)
        except Exception as e:
            print(" ",W.describe(w),W.length(w),"FAIL",type(e).__name__,str(e)[:100],"t=%.1f"%(time.time()-t),flush=True)
