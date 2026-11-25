
import sys, time
sys.path.insert(0, __import__("pathlib").Path(__file__).resolve().parent)
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl.hecke import HeckeAlgebra
from hecke_kl import poly

t0=time.time()
W=make("D",4); T=KLTable(W); H=HeckeAlgebra(W,T)
print("D4 order",W.order,"build %.2fs"%(time.time()-t0),flush=True)
# D4 bar invariance on all elements whose length <= 3
bad=[]; checked=0
for w in W.elements:
    if W.length(w)>3: continue
    c=H.Cprime(w); checked+=1
    if not H.equal(H.bar(c),c): bad.append(W.describe(w))
print("D4 C' bar sample checked",checked,"bad",len(bad),bad[:6],"t=%.1f"%(time.time()-t0),flush=True)
# R identity on sample
bad2=[]; ck=0
for y in W.elements:
    if W.length(y)>3: continue
    lhs=H.inverse_Tw(W.inverse(y)); rhs={}
    for x in W.elements:
        r=T.R(x,y)
        if not r: continue
        rhs=H.add(rhs,{x:poly.shift(poly.bar(r),-2*W.length(x))})
    ck+=1
    if not H.equal(lhs,rhs): bad2.append(W.describe(y))
print("D4 R identity sample checked",ck,"bad",len(bad2),bad2[:6],flush=True)
# bar involution: bar(bar(a))==a on arbitrary sparse elements (A3)
W3=make("A",3); H3=HeckeAlgebra(W3,KLTable(W3))
import random; rng=random.Random(3); bb=[]
for _ in range(20):
    a={w:poly.qpoly_from_coeffs([rng.randint(-2,2),rng.randint(-2,2)]) for w in rng.sample(W3.elements,3)}
    a={k:v for k,v in a.items() if v}
    if not H3.equal(H3.bar(H3.bar(a)),a): bb.append(1)
print("bar^2==id failures:",len(bb),"of 20",flush=True)
