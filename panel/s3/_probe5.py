
import sys, random
sys.path.insert(0, "/home/sudoshellz/Documents/hecke-kl")
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl.hecke import HeckeAlgebra
from hecke_kl import poly

rng = random.Random(12345)
results=[]

# ---- Test 1: bar is anti-involution on arbitrary sparse elements ----
for typ,n in [("A",4),("B",3),("G2",None),("D",4)]:
    W=make(typ,n); H=HeckeAlgebra(W, KLTable(W))
    bad_anti=[]; bad_mul=[]
    for _ in range(40):
        a={}; b={}
        for w in rng.sample(W.elements, min(4,len(W.elements))):
            a[w]=poly.qpoly_from_coeffs([rng.randint(-2,2),rng.randint(-2,2)])
        for w in rng.sample(W.elements, min(4,len(W.elements))):
            b[w]=poly.qpoly_from_coeffs([rng.randint(-2,2),rng.randint(-2,2)])
        a={k:v for k,v in a.items() if v}; b={k:v for k,v in b.items() if v}
        if not a or not b: continue
        ab=H.mul(a,b)
        lhs=H.bar(ab)
        rhs_anti=H.mul(H.bar(b),H.bar(a))
        rhs_hom=H.mul(H.bar(a),H.bar(b))
        if not H.equal(lhs,rhs_anti): bad_anti.append((W.describe(next(iter(a))),W.describe(next(iter(b)))))
        if not H.equal(lhs,rhs_hom): bad_mul.append((W.describe(next(iter(a))),W.describe(next(iter(b)))))
    results.append(f"bar_anti {W.name}: bad_anti={len(bad_anti)} bad_hom(expect fail)={len(bad_mul)}")

# ---- Test 2: mul associativity on sparse elements ----
for typ,n in [("A",4),("B",3),("D",4)]:
    W=make(typ,n); H=HeckeAlgebra(W, KLTable(W))
    bad=[]
    for _ in range(30):
        es=[]
        for _ in range(3):
            e={}
            for w in rng.sample(W.elements, min(3,len(W.elements))):
                e[w]=poly.qpoly_from_coeffs([rng.randint(-2,2)])
            es.append(e)
        l=H.mul(H.mul(es[0],es[1]),es[2]); r=H.mul(es[0],H.mul(es[1],es[2]))
        if not H.equal(l,r): bad.append(1)
    results.append(f"assoc {W.name}: bad={len(bad)}")

# ---- Test 3: Cprime bar-invariance on B4, D4 ----
for typ,n in [("B",4),("D",4)]:
    W=make(typ,n); H=HeckeAlgebra(W, KLTable(W))
    bad=[]
    for w in W.elements:
        c=H.Cprime(w)
        if not H.equal(H.bar(c),c): bad.append(W.describe(w))
    results.append(f"Cprime bar {W.name}: order={W.order} bad={len(bad)} {bad[:5]}")

for r in results: print(r, flush=True)
