
import sys, traceback
sys.path.insert(0, __import__("pathlib").Path(__file__).resolve().parent)
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl.hecke import HeckeAlgebra
from hecke_kl import poly

for typ,n in [("A",4),("B",3),("G2",None),("A",5)]:
    W=make(typ,n); T=KLTable(W); H=HeckeAlgebra(W,T)
    print("=== W",W.name,"order",W.order,flush=True)
    for w in W.elements:
        if W.length(w)==0: continue
        try:
            got=H.Cprime(w)
            bar=H.bar(got)
            if not H.equal(got,bar):
                print("BAR FAIL",W.describe(w),flush=True)
        except Exception as e:
            print("EXC",W.describe(w),type(e).__name__,str(e)[:120],flush=True)
            traceback.print_exc()
            break
    print("done",W.name,flush=True)
