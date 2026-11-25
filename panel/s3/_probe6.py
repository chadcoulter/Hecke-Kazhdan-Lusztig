
import sys, random
sys.path.insert(0, "/home/sudoshellz/Documents/hecke-kl")
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl.hecke import HeckeAlgebra
from hecke_kl import poly

rng = random.Random(12345)
for typ,n in [("A",4),("B",3),("G2",None)]:
    W=make(typ,n); H=HeckeAlgebra(W, KLTable(W))
    bad_anti=[]; bad_hom=[]
    for _ in range(30):
        a={}; b={}
        for w in rng.sample(W.elements, min(3,len(W.elements))):
            a[w]=poly.qpoly_from_coeffs([rng.randint(-2,2),rng.randint(-2,2)])
        for w in rng.sample(W.elements, min(3,len(W.elements))):
            b[w]=poly.qpoly_from_coeffs([rng.randint(-2,2),rng.randint(-2,2)])
        a={k:v for k,v in a.items() if v}; b={k:v for k,v in b.items() if v}
        if not a or not b: continue
        ab=H.mul(a,b)
        if not H.equal(H.bar(ab), H.mul(H.bar(b),H.bar(a))): bad_anti.append(1)
        if not H.equal(H.bar(ab), H.mul(H.bar(a),H.bar(b))): bad_hom.append(1)
    print(f"bar_anti {W.name}: bad_anti={len(bad_anti)} bad_hom={len(bad_hom)}",flush=True)
