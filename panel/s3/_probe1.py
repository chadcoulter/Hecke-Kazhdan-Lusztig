
import sys
sys.path.insert(0, "/home/sudoshellz/Documents/hecke-kl")
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl.hecke import HeckeAlgebra
from hecke_kl import poly
import sympy as sp

# 1) solve_kl_bar on S4 (A4) - cross-check requirement says at least S4,S5,B3
for typ,n in [("A",4),("A",5),("B",3)]:
    W=make(typ,n); T=KLTable(W); H=HeckeAlgebra(W,T)
    bad=[]; checked=0
    for w in W.elements:
        try:
            rec=H.solve_kl_bar(w)
        except Exception as e:
            bad.append((W.describe(w),"EXC",str(e)[:80])); continue
        for y in W.elements:
            if not W.leq(y,w): continue
            d=W.length(w)-W.length(y)
            cs=[1]
            for i in range(1,(d-1)//2+1):
                cs.append(rec.get((y,i),poly.zero()).get(0,0))
            exp=poly.qpoly_from_coeffs(cs)
            got=T.P(y,w)
            checked+=1
            if not poly.eq(got,exp):
                bad.append((W.describe(y),W.describe(w),poly.fmt(got),poly.fmt(exp)))
                if len(bad)>8: break
        if len(bad)>8: break
    print(f"[solve_kl_bar] {W.name} checked={checked} bad={len(bad)}",bad[:8],flush=True)
