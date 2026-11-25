
import sys
sys.path.insert(0, __import__("pathlib").Path(__file__).resolve().parent)
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl.hecke import HeckeAlgebra
from hecke_kl import poly
W=make("A",3); H=HeckeAlgebra(W,KLTable(W))
w=W.from_word((0,1))
c=H.Cprime(w)
print("C'_{01}:",[(W.describe(k),poly.fmt(v)) for k,v in c.items()])
e=W.identity
print("coeff at e before:",poly.fmt(c.get(e,{})))
bad=dict(c)
bad[e]=poly.add(bad.get(e,{}),{0:1})
print("coeff at e after :",poly.fmt(bad.get(e,{})))
print("bar(bad)==bad ?",H.equal(H.bar(bad),bad))
print("bar(bad) at e:",poly.fmt(H.bar(bad).get(e,{})))
# try a cleaner corruption: flip leading coefficient of a non-identity term
z=W.from_word((0,))
bad2=dict(c); bad2[z]=poly.add(bad2.get(z,{}),{0:1})
print("bad2 bar-invariant?",H.equal(H.bar(bad2),bad2),"(want False)")
