
import sys, time
sys.path.insert(0, "/home/sudoshellz/Documents/hecke-kl")
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl.hecke import HeckeAlgebra
from hecke_kl import poly
import random

t0=time.time()
W=make("D",4); T=KLTable(W); H=HeckeAlgebra(W,T)
print("D4 order",W.order,"build %.2f"%(time.time()-t0),flush=True)
bad_bar=[]; checked=0
for w in W.elements:
    c=H.Cprime(w); checked+=1
    if not H.equal(H.bar(c),c): bad_bar.append(W.describe(w))
print("D4 full C' bar: checked",checked,"bad",len(bad_bar),bad_bar[:8],"t=%.1f"%(time.time()-t0),flush=True)

# negative control: corrupt one coefficient -> bar MUST fail
w=W.from_word((0,1)); c=H.Cprime(w)
bad=dict(c); bad[W.identity]=poly.add(bad.get(W.identity,{}),{0:1})
print("[NEGCTRL] corrupted C' bar-invariant?",H.equal(H.bar(bad),bad),"(expect False)",flush=True)

# bar^2 and bar multiplicative (semilinear) on A4
W3=make("A",3); H3=HeckeAlgebra(W3,KLTable(W3))
rng=random.Random(3); bb=[]; sb=[]
for _ in range(12):
    a={}; b={}
    for w in rng.sample(W3.elements,2):
        a[w]=poly.qpoly_from_coeffs([rng.randint(-2,2),rng.randint(-2,2)])
    for w in rng.sample(W3.elements,2):
        b[w]=poly.qpoly_from_coeffs([rng.randint(-2,2),rng.randint(-2,2)])
    a={k:v for k,v in a.items() if v}; b={k:v for k,v in b.items() if v}
    if not a or not b: continue
    if not H3.equal(H3.bar(H3.bar(a)),a): bb.append(1)
    if not H3.equal(H3.bar(H3.mul(a,b)),H3.mul(H3.bar(a),H3.bar(b))): sb.append(1)
print("bar^2==id bad",len(bb),"of 12; bar(ab)=bar(a)bar(b) bad",len(sb),"of 12",flush=True)

# assoc sparse
ab=[]
for _ in range(10):
    es=[]
    for _ in range(3):
        e={w:poly.qpoly_from_coeffs([rng.randint(-2,2)]) for w in rng.sample(W3.elements,3)}
        es.append(e)
    if not H3.equal(H3.mul(H3.mul(es[0],es[1]),es[2]),H3.mul(es[0],H3.mul(es[1],es[2]))): ab.append(1)
print("assoc bad",len(ab),"of 10",flush=True)
print("DONE",flush=True)
