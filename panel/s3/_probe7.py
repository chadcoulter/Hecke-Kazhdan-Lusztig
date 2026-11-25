
import sys
sys.path.insert(0, "/home/sudoshellz/Documents/hecke-kl")
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl.hecke import HeckeAlgebra
from hecke_kl import poly

W=make("A",4); H=HeckeAlgebra(W,KLTable(W))
# pick noncommuting u,v from words
u=W.from_word((0,)); v=W.from_word((0,1))
# actually find pair with T_u T_v != T_v T_u
found=[]
for u in W.elements[:10]:
    for v in W.elements[:10]:
        if not H.equal(H.mul(H.T(u),H.T(v)), H.mul(H.T(v),H.T(u))):
            found.append((u,v))
    if found: break
u,v=found[0]
print("u,v",W.describe(u),W.describe(v),"noncommuting",len(found))
puv=H.mul(H.T(u),H.T(v)); pvu=H.mul(H.T(v),H.T(u))
lhs=H.bar(puv)
anti=H.mul(H.bar(H.T(v)),H.bar(H.T(u)))
hom=H.mul(H.bar(H.T(u)),H.bar(H.T(v)))
print("bar(TuTv)==bar(Tv)bar(Tu) [anti]:",H.equal(lhs,anti))
print("bar(TuTv)==bar(Tu)bar(Tv) [hom] :",H.equal(lhs,hom))
print("bar(TuTv) keys:",[(W.describe(k),poly.fmt(c)) for k,c in lhs.items()])
print("anti keys       :",[(W.describe(k),poly.fmt(c)) for k,c in anti.items()])
print("hom keys        :",[(W.describe(k),poly.fmt(c)) for k,c in hom.items()])
