#!/usr/bin/env python3
"""S3 Hecke/module-layer adversarial probe.

Every PASS/FAIL below is an assertion on a value produced by this code.
Run:  python3 panel/s3/probe.py
"""
import subprocess, sys, time, random
sys.path.insert(0, "/home/sudoshellz/Documents/hecke-kl")
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl.hecke import HeckeAlgebra
from hecke_kl import poly

ROOT = "/home/sudoshellz/Documents/hecke-kl"
def say(*a): print(*a, flush=True)

# ---- 0. ground truth must remain 21/21 ------------------------------------
gt = subprocess.run([sys.executable, "ground_truth.py"], cwd=ROOT,
                    capture_output=True, text=True, timeout=300)
summ = [x for x in gt.stdout.splitlines() if x.startswith("SUMMARY")]
say("[GT]", summ, "exit=", gt.returncode)
assert gt.returncode == 0 and summ and "FAIL=0" in summ[0]

# ---- 1. D4 (order 192) -- NOT in ground_truth -------------------------------
W = make("D", 4); T = KLTable(W); H = HeckeAlgebra(W, T)
say("[D4] order", W.order)

# full bar-invariance of C'_w (all 192 w)
t0 = time.time(); bad_bar = []
for w in W.elements:
    c = H.Cprime(w)
    if not H.equal(H.bar(c), c): bad_bar.append(W.describe(w))
say("[D4] bar(C'_w)=C'_w: checked", len(W.elements), "bad", len(bad_bar),
    bad_bar[:6], "t=%.1fs" % (time.time()-t0))
assert not bad_bar

# full polynomial invariants: nonneg, degree bound, constant term 1, w0
t0 = time.time(); vals = T.all_polynomials()
bad_nonneg = bad_deg = bad_w0 = 0
for (y, w), p in vals.items():
    if any(c < 0 for c in p.values()): bad_nonneg += 1
    if y != w and W.leq(y, w):
        d = W.length(w) - W.length(y)
        if poly.qdegree(p) > (d-1)//2 or p.get(0, 0) != 1: bad_deg += 1
for y in W.elements:
    if not poly.eq(T.P(y, W.longest), poly.one()): bad_w0 += 1
say("[D4] polynomials", len(vals), "nonneg_bad", bad_nonneg, "deg_bad", bad_deg,
    "w0_bad", bad_w0, "t=%.1fs" % (time.time()-t0))
assert bad_nonneg == bad_deg == bad_w0 == 0

# R identity on a representative sample (all w of length <= 3)
t0 = time.time(); bad_R = []; ck_R = 0
for y in W.elements:
    if W.length(y) > 3: continue
    lhs = H.inverse_Tw(W.inverse(y)); rhs = {}
    for x in W.elements:
        r = T.R(x, y)
        if r: rhs = H.add(rhs, {x: poly.shift(poly.bar(r), -2*W.length(x))})
    ck_R += 1
    if not H.equal(lhs, rhs): bad_R.append(W.describe(y))
say("[D4] R identity checked", ck_R, "bad", len(bad_R), bad_R[:6],
    "t=%.1fs" % (time.time()-t0))
assert not bad_R

# ---- 2. solve_kl_bar independent cross-check -------------------------------
for typ, n in [("A", 3), ("A", 4), ("B", 3)]:
    Wk = make(typ, n); Hk = HeckeAlgebra(Wk, KLTable(Wk))
    t0 = time.time(); bad = []; ck = 0
    for w in Wk.elements:
        try: rec = Hk.solve_kl_bar(w)
        except Exception as e:
            bad.append((Wk.describe(w), type(e).__name__, str(e)[:60])); continue
        for y in Wk.elements:
            if not Wk.leq(y, w) or y == w: continue
            d = Wk.length(w) - Wk.length(y)
            cs = [1] + [rec.get((y, i), poly.zero()).get(0, 0)
                        for i in range(1, (d-1)//2 + 1)]
            exp = poly.qpoly_from_coeffs(cs); got = Hk.kl.P(y, w); ck += 1
            if not poly.eq(got, exp):
                bad.append((Wk.describe(y), Wk.describe(w), poly.fmt(got), poly.fmt(exp)))
    say("[solve_kl_bar]", Wk.name, "pairs", ck, "bad", len(bad), bad[:5],
        "t=%.1fs" % (time.time()-t0))
    assert not bad

# ---- 3. algebra laws on sparse elements ------------------------------------
Wk = make("A", 3); Hk = HeckeAlgebra(Wk, KLTable(Wk))
rng = random.Random(7); bad_assoc = bad_bar2 = bad_semilin = 0
for _ in range(10):
    es = []
    for _ in range(3):
        e = {z: poly.qpoly_from_coeffs([rng.randint(-2, 2)])
             for z in rng.sample(Wk.elements, 3)}
        es.append(e)
    if not Hk.equal(Hk.mul(Hk.mul(es[0], es[1]), es[2]),
                    Hk.mul(es[0], Hk.mul(es[1], es[2]))): bad_assoc += 1
    if not Hk.equal(Hk.bar(Hk.bar(es[0])), es[0]): bad_bar2 += 1
    if not Hk.equal(Hk.bar(Hk.mul(es[0], es[1])),
                    Hk.mul(Hk.bar(es[0]), Hk.bar(es[1]))): bad_semilin += 1
say("[algebra] assoc_bad", bad_assoc, "bar^2_bad", bad_bar2,
    "bar(ab)=bar(a)bar(b)_bad", bad_semilin, "of 10")
assert bad_assoc == bad_bar2 == bad_semilin == 0

# ---- 4. bar is SEMILINEAR HOMOMORPHISM, not anti-involution ----------------
u = Wk.from_word((0,)); v = Wk.from_word((1,))
puv = Hk.mul(Hk.T(u), Hk.T(v))
is_hom = Hk.equal(Hk.bar(puv), Hk.mul(Hk.bar(Hk.T(u)), Hk.bar(Hk.T(v))))
is_anti = Hk.equal(Hk.bar(puv), Hk.mul(Hk.bar(Hk.T(v)), Hk.bar(Hk.T(u))))
say("[convention] bar(T_u T_v)=bar(T_u)bar(T_v):", is_hom,
    " anti-version:", is_anti, "(standard KL79 = homomorphism)")
assert is_hom and not is_anti

# ---- 5. negative controls --------------------------------------------------
w = Wk.from_word((0, 1)); c = Hk.Cprime(w)
corrupt = dict(c); z = Wk.from_word((0,))
corrupt[z] = poly.add(corrupt.get(z, {}), {0: 1})
say("[NEGCTRL] corrupted C' is bar-invariant?", Hk.equal(Hk.bar(corrupt), corrupt),
    "(expected False)")
assert not Hk.equal(Hk.bar(corrupt), corrupt)

class BadKL(KLTable):
    def P(self, y, w):
        p = super().P(y, w)
        if y != w: p = poly.add(p, {0: 1})   # corrupt constant term
        return p
Wb = make("A", 3); Hb = HeckeAlgebra(Wb, BadKL(Wb))
wb = Wb.from_word((0, 1))
say("[NEGCTRL] corrupted P breaks bar(C')?",
    not Hb.equal(Hb.bar(Hb.Cprime(wb)), Hb.Cprime(wb)), "(expected True)")
assert not Hb.equal(Hb.bar(Hb.Cprime(wb)), Hb.Cprime(wb))

# ---- 6. mu-graph left-cell prototype ---------------------------------------
def mu_graph_cells(W, T):
    par = {w: w for w in W.elements}
    def find(a):
        while par[a] != a: par[a] = par[par[a]]; a = par[a]
        return a
    edges = 0
    for w in W.elements:
        for y in W.elements:
            if y == w or W.length(y) > W.length(w) or not W.leq(y, w): continue
            if (W.length(w)-W.length(y)) % 2 == 0: continue
            if T.mu(y, w) != 0:
                edges += 1; ra, rb = find(y), find(w)
                if ra != rb: par[rb] = ra
    comps = {}
    for w in W.elements: comps.setdefault(find(w), []).append(w)
    return edges, sorted(len(c) for c in comps.values())
for typ, n in [("A", 3), ("A", 4), ("B", 2), ("G2", None)]:
    Wm = make(typ, n); Tm = KLTable(Wm); Tm.all_polynomials()
    e, sizes = mu_graph_cells(Wm, Tm)
    say("[mu-graph]", Wm.name, "order", Wm.order, "edges", e,
        "components", len(sizes), "sizes", sizes)
say("DONE -- all assertions passed")
