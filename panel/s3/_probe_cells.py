
import sys
sys.path.insert(0, "/home/sudoshellz/Documents/hecke-kl")
from hecke_kl.groups import make
from hecke_kl.kl import KLTable

# mu-graph prototype: vertices = W, undirected edge x--y iff mu(min,max)!=0.
# Connected components approximate KL LEFT cells.
def mu_graph_cells(W, T):
    parent = {w: w for w in W.elements}
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    def union(a,b):
        ra,rb = find(a),find(b)
        if ra != rb: parent[rb] = ra
    edges = 0
    for w in W.elements:
        for y in W.elements:
            if y == w or W.length(y) > W.length(w): continue
            if not W.leq(y,w): continue
            if (W.length(w)-W.length(y)) % 2 == 0: continue
            if T.mu(y,w) != 0:
                union(y,w); edges += 1
    comps = {}
    for w in W.elements: comps.setdefault(find(w), []).append(w)
    return edges, list(comps.values())

for typ,n in [("A",3),("A",4),("B",2),("G2",None)]:
    W=make(typ,n); T=KLTable(W); T.all_polynomials()
    e,cells=mu_graph_cells(W,T)
    print(f"[mu-graph] {W.name} order={W.order} edges={e} components={len(cells)} sizes={sorted(len(c) for c in cells)}",flush=True)
