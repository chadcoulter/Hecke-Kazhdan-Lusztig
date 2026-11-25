
import sys, time
sys.path.insert(0, __import__("pathlib").Path(__file__).resolve().parent)
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl.hecke import HeckeAlgebra
from hecke_kl import poly

for typ,n in [("A",3),("A",4),("B",3)]:
    W=make(typ,n); H=HeckeAlgebra(W,KLTable(W))
    print("=== W",W.name,"order",W.order,flush=True)
    tested=0
    for w in W.elements:
        t=time.time()
        try:
            r=H.solve_kl_bar(w)
        except Exception as e:
            print(" FAIL",W.describe(w),"len",W.length(w),type(e).__name__,str(e)[:90],"t=%.1f"%(time.time()-t),flush=True)
            if W.length(w)<=3: tested+=1
            if tested>6: break
            continue
        bad=[]
        for y in W.elements:
            if not W.leq(y,w): continue
            d=W.length(w)-W.length(y); cs=[1]
            for i in range(1,(d-1)//2+1): cs.append(r.get((y,i),poly.zero()).get(0,0))
            exp=poly.qpoly_from_coeffs(cs); got=H.kl.P(y,w)
            if not poly.eq(got,exp): bad.append((W.describe(y),poly.fmt(got),poly.fmt(exp)))
        if bad: print(" MISMATCH",W.describe(w),"len",W.length(w),bad[:4],flush=True)
    print(" done",W.name,flush=True)
