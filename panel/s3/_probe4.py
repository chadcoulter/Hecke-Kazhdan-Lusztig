
import sys
sys.path.insert(0, __import__("pathlib").Path(__file__).resolve().parent)
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl.hecke import HeckeAlgebra
from hecke_kl import poly

W=make("A",5); H=HeckeAlgebra(W, KLTable(W))
bad=[]
for w in W.elements:
    c=H.Cprime(w); bar=H.bar(c)
    if not H.equal(c,bar):
        bad.append(W.describe(w))
        if len(bad)>=10: break
print("A5 bar failures:",len(bad),bad[:10],flush=True)
