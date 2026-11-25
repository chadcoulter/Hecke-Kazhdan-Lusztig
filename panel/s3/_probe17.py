
import sys, time
sys.path.insert(0, "/home/sudoshellz/Documents/hecke-kl")
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl import poly
W=make("D",4); T=KLTable(W)
t0=time.time(); vals=T.all_polynomials(); print("D4 all_polynomials",len(vals),"t=%.1f"%(time.time()-t0),flush=True)
bad=[]; nonneg=[]
for (y,w),p in vals.items():
    if any(c<0 for c in p.values()): nonneg.append((W.describe(y),W.describe(w)))
    if y!=w and W.leq(y,w):
        d=W.length(w)-W.length(y)
        if poly.qdegree(p)>(d-1)//2 or p.get(0,0)!=1:
            bad.append((W.describe(y),W.describe(w),poly.fmt(p)))
print("D4 nonneg bad",len(nonneg),nonneg[:5],"degree bad",len(bad),bad[:5],flush=True)
# longest element
badw=[(W.describe(y),poly.fmt(T.P(y,W.longest))) for y in W.elements if not poly.eq(T.P(y,W.longest),poly.one())]
print("D4 w0 all P=1 bad",len(badw),badw[:5],flush=True)
