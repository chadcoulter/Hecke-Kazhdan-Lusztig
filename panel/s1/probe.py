#!/usr/bin/env python3
"""S1 adversarial audit of the FROZEN Hecke-KL group layer.

Attacks Cartan/Tits conventions, BFS canonical words, left/right descents,
Bruhat closure and reflection enumeration using an INDEPENDENT geometric
(signed-permutation / root-system) model of each finite Coxeter group.

All assertions are on values actually computed by code.  Usage:

    python3 probe.py            # full audit, includes A6 (slow build ~6 min)
    python3 probe.py --quick    # skip A6 (development / fast feedback)
"""
from __future__ import annotations

import itertools
import math
import os
import random
import subprocess
import sys
import time
from fractions import Fraction as F

ROOT = "/home/sudoshellz/Documents/hecke-kl"
sys.path.insert(0, ROOT)

from hecke_kl.coxeter import CoxeterGroup
from hecke_kl.groups import cartan_a, cartan_d, make

QUICK = "--quick" in sys.argv
RNG = random.Random(20261125)
RESULTS: list[tuple[str, bool, str]] = []
FINDINGS: list[tuple[str, str]] = []


def check(name: str, cond: bool, detail: str = "") -> bool:
    cond = bool(cond)
    RESULTS.append((name, cond, detail))
    print(f"{'PASS' if cond else 'FAIL'}  {name}  {detail}", flush=True)
    return cond


def finding(name: str, detail: str) -> None:
    """Document a latent edge defect without aborting the probe (the frozen
    core cannot be edited from this seat; the parent needs the reproducer)."""
    FINDINGS.append((name, detail))
    print(f"FINDING  {name}  {detail}", flush=True)


# --------------------------------------------------------------------------
# Independent geometric model.  Elements act on the ambient e_basis by signed
# permutations g(e_i)=sign_i*e_{p_i}, stored as tuple((sign, image_index),...).
# Composition gcomp(g,h) applies h first, then g (matches core g @ h).
# --------------------------------------------------------------------------
def geom_generators(typ: str, n: int) -> tuple[list[tuple], int]:
    if typ == "A":
        m = n + 1
        gens = []
        for i in range(n):
            g = [(1, j) for j in range(m)]
            g[i], g[i + 1] = (1, i + 1), (1, i)
            gens.append(tuple(g))
        return gens, m
    if typ == "B":
        m = n
        gens = []
        for i in range(n - 1):
            g = [(1, j) for j in range(m)]
            g[i], g[i + 1] = (1, i + 1), (1, i)
            gens.append(tuple(g))
        g = [(1, j) for j in range(m)]
        g[n - 1] = (-1, n - 1)                        # short simple root e_{n-1}
        gens.append(tuple(g))
        return gens, m
    if typ == "D":
        m = n
        gens = []
        for i in range(n - 1):
            g = [(1, j) for j in range(m)]
            g[i], g[i + 1] = (1, i + 1), (1, i)
            gens.append(tuple(g))
        g = [(1, j) for j in range(m)]
        g[n - 2], g[n - 1] = (-1, n - 1), (-1, n - 2)  # root e_{n-2}+e_{n-1}
        gens.append(tuple(g))
        return gens, m
    raise ValueError(typ)


def gcomp(g, h):
    return tuple((h[i][0] * g[h[i][1]][0], g[h[i][1]][1]) for i in range(len(g)))


def gident(m):
    return tuple((1, i) for i in range(m))


def ginverse(g):
    inv = []
    for p in range(len(g)):
        slot = next(i for i, (_, q) in enumerate(g) if q == p)
        inv.append((g[slot][0], slot))   # each sign is its own inverse
    return tuple(inv)


class GeomGroup:
    def __init__(self, typ, n):
        self.typ, self.n = typ, n
        self.gens, self.m = geom_generators(typ, n)
        e = gident(self.m)
        word_of = {e: ()}
        order = [e]
        qi = 0
        while qi < len(order):
            g = order[qi]
            qi += 1
            for i, s in enumerate(self.gens):
                ng = gcomp(g, s)
                if ng not in word_of:
                    word_of[ng] = word_of[g] + (i,)
                    order.append(ng)
        self.elements = order
        self.word_of = word_of
        self.length_of = {g: len(w) for g, w in word_of.items()}
        self.inv = {g: ginverse(g) for g in order}
        self.longest = max(order, key=lambda g: self.length_of[g])
        self.nreflections = len(self.reflections())

    def reflections(self):
        out = set()
        for w in self.elements:
            wi = self.inv[w]
            for s in self.gens:
                out.add(gcomp(gcomp(w, s), wi))
        return out

    def left_descents(self, g):
        return tuple(i for i, s in enumerate(self.gens)
                     if self.length_of[gcomp(s, g)] < self.length_of[g])

    def right_descents(self, g):
        return tuple(i for i, s in enumerate(self.gens)
                     if self.length_of[gcomp(g, s)] < self.length_of[g])

    def lex_least_word(self, g):
        """Independent greedy reconstruction: the smallest left descent is
        the first letter of the lex-least reduced word; recurse on s_i*g."""
        v = g
        out = []
        e = gident(self.m)
        L = self.length_of[g]
        while v != e:
            ds = self.left_descents(v)
            assert ds, "non-identity element with no left descent"
            i = min(ds)
            out.append(i)
            v = gcomp(self.gens[i], v)
            assert self.length_of[v] == L - len(out)
        return tuple(out)

    def rank_matrix(self, g):
        one = [p for _, p in g]
        N = len(one)
        return [[sum(1 for i in range(p) if one[i] <= q) for q in range(N)]
                for p in range(1, N + 1)]

    def bruhat_rank(self, u, w):
        """Type-A criterion: u <= w iff r_u >= r_w entrywise."""
        ru, rw = self.rank_matrix(u), self.rank_matrix(w)
        return all(ru[i][j] >= rw[i][j]
                   for i in range(len(ru)) for j in range(len(ru[0])))

    def subword_down(self, w):
        """Bruhat down-set via the subword property (B/D).  Idempotent
        right-multiplication closure of {e} under the letters of a reduced
        word for w; deduplication keeps each set <= |W|."""
        cur = {gident(self.m)}
        for i in self.word_of[w]:
            s = self.gens[i]
            cur |= {gcomp(p, s) for p in cur}
        return frozenset(cur)


# Change-of-basis matrices from ambient e coordinates to simple-root
# coordinates (used only for the reflection-matrix cross-check).
def root_change_of_basis(typ, n, m):
    """Return (P, L) with P[:,j] = alpha_j in e-coordinates and
    L = (P^T P)^{-1} P^T, so coordinates in the simple-root basis are L*e."""
    alphas = []
    if typ == "A":
        for i in range(n):
            a = [0] * m; a[i] = 1; a[i + 1] = -1; alphas.append(a)
    if typ == "B":
        for i in range(n - 1):
            a = [0] * m; a[i] = 1; a[i + 1] = -1; alphas.append(a)
        # cartan_b uses A_{n-1,n}=-2, i.e. alpha_n = 2*e_{n-1} (long-root
        # convention for B_n).  Must scale accordingly or the change of
        # basis is wrong (this bit me once -- see report.md).
        a = [0] * m; a[n - 1] = 2; alphas.append(a)
    elif typ == "D":
        for i in range(n - 1):
            a = [0] * m; a[i] = 1; a[i + 1] = -1; alphas.append(a)
        a = [0] * m; a[n - 2] = 1; a[n - 1] = 1; alphas.append(a)
    P = [[F(alpha[j]) for alpha in alphas] for j in range(m)]
    A = [[sum(P[k][i] * P[k][j] for k in range(m)) for j in range(n)]
         for i in range(n)]
    M = [row[:] + [F(int(i == j)) for j in range(n)] for i, row in enumerate(A)]
    for c in range(n):
        piv = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[piv] = M[piv], M[c]
        d = M[c][c]
        M[c] = [x // d if isinstance(x, int) else x / d for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [M[r][j] - f * M[c][j] for j in range(2 * n)]
    Ainv = [row[n:] for row in M]
    L = [[sum(Ainv[i][k] * P[j][k] for k in range(n)) for j in range(m)]
         for i in range(n)]
    return P, L


def e_matrix(r):
    """Ambient e-basis integer matrix of signed permutation r."""
    return tuple(tuple(r[j][0] if k == r[j][1] else 0 for j in range(len(r)))
                 for k in range(len(r)))


def known_order(typ, n):
    return {"A": math.factorial(n + 1),
            "B": 2 ** n * math.factorial(n),
            "D": 2 ** (n - 1) * math.factorial(n)}[typ]


# --------------------------------------------------------------------------
# Audit one group against the independent model.
# --------------------------------------------------------------------------
def audit(typ, n, expected_w0len, full_bruhat, bruhat_samples, word_samples):
    tag = f"{typ}{n}"
    t0 = time.time()
    W = make(typ, n)
    G = GeomGroup(typ, n)
    print(f"  [{tag}] core order={W.order} geometric order={len(G.elements)} "
          f"build={time.time()-t0:.2f}s", flush=True)

    check(f"{tag} order == known finite order",
          W.order == known_order(typ, n) == len(G.elements),
          f"core={W.order} known={known_order(typ,n)} geometric={len(G.elements)}")
    check(f"{tag} longest length == #positive roots",
          W.length(W.longest) == expected_w0len == G.length_of[G.longest],
          f"core={W.length(W.longest)} expected={expected_w0len} "
          f"geometric={G.length_of[G.longest]}")

    # Align by canonical word, not discovery index, so a word mismatch
    # cannot silently corrupt the length/descent comparisons below.
    core_by_word = {W.word(g): g for g in W.elements}
    geo_by_word = {G.word_of[h]: h for h in G.elements}
    check(f"{tag} canonical-word key sets agree",
          set(core_by_word) == set(geo_by_word) and
          len(core_by_word) == len(geo_by_word) == W.order,
          f"core_keys={len(core_by_word)} geo_keys={len(geo_by_word)} "
          f"symdiff={len(set(core_by_word) ^ set(geo_by_word))}")
    pairs = [(core_by_word[word], geo_by_word[word])
             for word in core_by_word if word in geo_by_word]

    bad_len, bad_desc = [], []
    for g, h in pairs:
        if W.length(g) != G.length_of[h]:
            bad_len.append((W.describe(g), W.length(g), G.length_of[h]))
        lcore, rcore = W.left_descents(g), W.right_descents(g)
        lgeo, rgeo = G.left_descents(h), G.right_descents(h)
        if (lcore, rcore) != (lgeo, rgeo):
            bad_desc.append((W.describe(g), (lcore, rcore), (lgeo, rgeo)))
    check(f"{tag} all lengths agree", not bad_len,
          f"mismatches={len(bad_len)} first={bad_len[:2]}")
    check(f"{tag} all left/right descents agree", not bad_desc,
          f"mismatches={len(bad_desc)} first={bad_desc[:2]}")

    sample = pairs if word_samples is None else RNG.sample(pairs, word_samples)
    bad_lex = []
    for g, h in sample:
        lx = G.lex_least_word(h)
        if W.word(g) != lx:
            bad_lex.append((W.describe(g), W.word(g), lx))
    check(f"{tag} canonical word == greedy lex-least reduced word", not bad_lex,
          f"checked={len(sample)} mismatches={len(bad_lex)} first={bad_lex[:2]}")

    check(f"{tag} reflection count == #positive roots",
          len(W.reflections) == expected_w0len == G.nreflections,
          f"core={len(W.reflections)} expected={expected_w0len} "
          f"geometric={G.nreflections}")

    # Strong reflection check: conjugate EVERY simple reflection by EVERY
    # element in both models, convert geometric matrices into simple-root
    # coordinates, and compare the two reflection sets exactly.
    P, L = root_change_of_basis(typ, n, G.m)
    def to_root(M):
        MP = [[sum(M[i][kk] * P[kk][j] for kk in range(G.m)) for j in range(n)]
              for i in range(G.m)]
        return tuple(tuple(int(sum(L[i][kk] * MP[kk][j]
                                 for kk in range(G.m)))
                           for j in range(n)) for i in range(n))
    gset = {to_root(e_matrix(r)) for r in G.reflections()}
    cset = set(W.reflections)
    check(f"{tag} reflection matrices == geometric conjugates",
          cset == gset,
          f"core={len(cset)} geo={len(gset)} core_only={len(cset-gset)} "
          f"geo_only={len(gset-cset)}")
    # Also verify the core's S_i themselves match the geometric simple
    # reflections after the same change of basis.
    simple_bad = [i for i in range(W.n)
                  if to_root(e_matrix(G.gens[i])) != W.S[i]]
    check(f"{tag} simple reflection matrices agree", not simple_bad,
          f"bad_indices={simple_bad}")

    # Bruhat: full pair check on small groups, sampled on large ones.
    if full_bruhat:
        if typ == "A":
            expected = {(g1, g2): G.bruhat_rank(h1, h2)
                        for g1, h1 in pairs for g2, h2 in pairs}
        else:
            downs = {h2: G.subword_down(h2) for _, h2 in pairs}
            expected = {(g1, g2): (h1 in downs[h2])
                        for g1, h1 in pairs for g2, h2 in pairs}
        test_pairs = [(g1, g2) for g1, _ in pairs for g2, _ in pairs]
    else:
        test_pairs = [(pairs[i][0], pairs[j][0])
                      for i, j in ((RNG.randrange(W.order), RNG.randrange(W.order))
                                   for _ in range(bruhat_samples))]
        expected = {}
        for g1, g2 in test_pairs:
            h1 = geo_by_word[W.word(g1)]
            h2 = geo_by_word[W.word(g2)]
            expected[(g1, g2)] = (G.bruhat_rank(h1, h2) if typ == "A"
                                  else h1 in G.subword_down(h2))
    gbads = [(W.describe(g1), W.describe(g2), W.leq(g1, g2), expected[(g1, g2)])
             for g1, g2 in test_pairs
             if W.leq(g1, g2) != expected[(g1, g2)]]
    check(f"{tag} Bruhat agrees with independent criterion", not gbads,
          f"pairs={len(test_pairs)} mismatches={len(gbads)} first={gbads[:3]}")

    check(f"{tag} core longest maps to geometric longest",
          G.word_of[G.longest] == W.word(W.longest),
          f"core_word={W.word(W.longest)} geo_word={G.word_of[G.longest]}")
    return W, G, pairs


# --------------------------------------------------------------------------
# Negative controls: a deliberately corrupted Cartan must be DETECTED.
# --------------------------------------------------------------------------
def negative_control():
    bad = cartan_a(4)
    bad[0][1] = bad[1][0] = 0          # A3 -> A1 x A3 (order 48, not 120)
    Wbad = CoxeterGroup("A4_CORRUPT", bad)
    detected = (Wbad.order != 120) or (Wbad.length(Wbad.longest) != 6)
    check("NEGCTRL corrupted Cartan (A4 edge removed) FAILS probe", detected,
          f"corrupt order={Wbad.order} w0len={Wbad.length(Wbad.longest)} "
          f"(good A4: order=120 w0len=6)")

    bad2 = cartan_d(4)
    bad2[0][1] = bad2[1][0] = 0        # D4 -> A1 x A3 (order 48, not 192)
    Wbad2 = CoxeterGroup("D4_CORRUPT", bad2)
    detected2 = (Wbad2.order != 192) or (Wbad2.length(Wbad2.longest) != 12)
    check("NEGCTRL corrupted D4 Cartan FAILS probe", detected2,
          f"corrupt order={Wbad2.order} w0len={Wbad2.length(Wbad2.longest)} "
          f"(good D4: order=192 w0len=12)")


# --------------------------------------------------------------------------
# Re-run ground truth from inside the probe.
# --------------------------------------------------------------------------
def ground_truth():
    p = subprocess.run([sys.executable, os.path.join(ROOT, "ground_truth.py")],
                       capture_output=True, text=True, timeout=600, cwd=ROOT)
    summary = [ln for ln in p.stdout.splitlines() if ln.startswith("SUMMARY")]
    ok = p.returncode == 0 and any("PASS=21" in s and "FAIL=0" in s for s in summary)
    check("ground_truth.py re-run: 21/21 PASS", ok,
          f"rc={p.returncode} summary={summary} stderr_tail={p.stderr[-200:]!r}")
    return p.stdout


# --------------------------------------------------------------------------
# Edge-case / boundary audit.  Degenerate ranks, identity/longest endpoints,
# input validation, Laurent-polynomial corners, and Hecke arithmetic at T_e.
# --------------------------------------------------------------------------
def edge_cases():
    from hecke_kl import poly
    from hecke_kl.kl import KLTable
    from hecke_kl.hecke import HeckeAlgebra

    # 1. Degenerate rank A1: group of order 2, one reflection.
    W1 = make("A", 1)
    check("EDGE A1 order=2 l(w0)=1 one reflection",
          W1.order == 2 and W1.length(W1.longest) == 1
          and len(W1.reflections) == 1 and W1.word(W1.longest) == (0,),
          f"order={W1.order} w0len={W1.length(W1.longest)} "
          f"refl={len(W1.reflections)} w0={W1.word(W1.longest)}")

    # 2. B2 == I2(4) and G2 == I2(6): aliases must be the same group.
    B2, I4 = make("B", 2), make("I2", 4)
    check("EDGE B2 is I2(4): Cartan + order + word set agree",
          B2.cartan == I4.cartan and B2.order == I4.order == 8
          and {B2.word(g) for g in B2.elements} == {I4.word(g) for g in I4.elements},
          f"cartan_equal={B2.cartan == I4.cartan} orders=({B2.order},{I4.order})")
    G2, I6 = make("G2"), make("I2", 6)
    check("EDGE G2 is I2(6): order 12 + word set agree",
          G2.order == I6.order == 12 and
          {G2.word(g) for g in G2.elements} == {I6.word(g) for g in I6.elements},
          f"orders=({G2.order},{I6.order})")

    # 3. cartan_g2() vs make('G2') are transposes but give order-12 groups.
    from hecke_kl.groups import cartan_g2
    cg = cartan_g2()
    Gt = CoxeterGroup("G2_transpose", cg)
    check("EDGE cartan_g2 transpose variant still finite order 12",
          G2.cartan == tuple(tuple(cg[j][i] for j in range(2)) for i in range(2))
          and Gt.order == 12,
          f"make cartan={G2.cartan} cartan_g2={cg} transpose_order={Gt.order}")

    # 4. Invalid parameters must raise, never silently produce a group.
    from hecke_kl.groups import cartan_a, cartan_b, cartan_i2
    rejects = []
    for fn, args in [(cartan_i2, (5,)), (cartan_i2, (7,)),
                     (cartan_a, (0,)), (cartan_b, (1,)), (make, ("D", 3))]:
        try:
            fn(*args)
            rejects.append((fn.__name__, args, "NO ERROR"))
        except ValueError:
            pass
    check("EDGE invalid Cartan parameters rejected (i2 5/7, A0, B1, D3)",
          not rejects, f"unexpected={rejects}")

    # 5. Identity/longest structural boundaries on every tested group.
    edge_bad = []
    for typ, n in [("A", 1), ("A", 2), ("B", 2), ("D", 4), ("B", 4)]:
        W = make(typ, n)
        e, w0 = W.identity, W.longest
        if not (W.left_descents(e) == W.right_descents(e) == ()
                and W.inverse(e) == e and W.inverse(w0) == w0):
            edge_bad.append(W.name)
        if any(W.right_descents(g) != W.left_descents(W.inverse(g))
               for g in W.elements):
            edge_bad.append(W.name + ":desc-inverse")
        if any(not W.leq(g, g) or (W.leq(g, h) and W.leq(h, g) and g != h)
               for g in W.elements for h in W.elements):
            edge_bad.append(W.name + ":bruhat-reflexivity")
        if not all(W.mul(t, t) == e for t in W.reflections):
            edge_bad.append(W.name + ":refl-involution")
        try:
            W.from_word((W.n,))
            edge_bad.append(W.name + ":bad-index-accepted")
        except ValueError:
            pass
        try:
            W.from_word((-1,))
            edge_bad.append(W.name + ":negative-index-accepted")
        except ValueError:
            pass
    check("EDGE identity/longest descents, inverse, Bruhat reflexivity, "
          "reflection involutions, from_word validation",
          not edge_bad, f"bad={edge_bad}")

    # 6. Laurent-polynomial corners: bar involution, clean, q_power guard,
    #    as_q_list rejects odd exponents.
    check("EDGE poly.bar involutive on negative exponents",
          poly.bar(poly.bar({-3: 2, 0: 5, 4: -1})) == {-3: 2, 0: 5, 4: -1},
          f"bar(bar(p)) == p")
    try:
        poly.q_power(-1)
        check("EDGE q_power(-1) rejected", False, "no error")
    except ValueError:
        check("EDGE q_power(-1) rejected", True, "ValueError")
    try:
        poly.as_q_list({1: 2})
        check("EDGE as_q_list rejects odd v-exponent", False, "no error")
    except ValueError:
        check("EDGE as_q_list rejects odd v-exponent", True, "ValueError")
    check("EDGE poly.zero()/bar/shift on empty poly",
          poly.bar({}) == {} and poly.shift({}, -5) == {}
          and poly.mul({0: 1}, {}) == {},
          "empty polynomials stable under bar/shift/zero-mul")

    # 7. KL/Hecke at identity: P_{e,s}=1, P_{s,e}=0, mu(e,s)=1, R(e,s)=0,
    #    solve_kl_bar(e) = {} (no unknowns).
    W = make("A", 2)
    T, H = KLTable(W), None
    H = HeckeAlgebra(W, T)
    e, s0 = W.identity, W.S[0]
    check("EDGE KL/Hecke identity values",
          poly.fmt(T.P(e, s0)) == "1" and poly.eq(T.P(s0, e), poly.zero())
          and T.mu(e, s0) == 1 and poly.eq(T.R(s0, e), poly.zero())
          and poly.eq(T.R(e, e), poly.one()) and H.solve_kl_bar(e) == {},
          f"P(e,s)={poly.fmt(T.P(e,s0))} P(s,e)={poly.fmt(T.P(s0,e))} "
          f"mu(e,s)={T.mu(e,s0)} R(e,e)={poly.fmt(T.R(e,e))} "
          f"solve(e)={H.solve_kl_bar(e)}")
    check("EDGE Hecke T_e unit and empty-element arithmetic",
          H.equal(H.mul(H.T(e), H.T(e)), H.T(e))
          and H.mul(H.T(e), {}) == {} and H.mul({}, H.T(e)) == {}
          and H.bar({}) == {} and H.equal(H.bar(H.T(e)), H.T(e)),
          "T_e*T_e=T_e, T_e*0=0*T_e=0, bar(0)=0, bar(T_e)=T_e")

    # 8. strictly_below_with_descents vs brute-force definition on small
    #    groups (identity/low-rank boundary especially).
    sb_bad = []
    for typ, n in [("A", 1), ("A", 2), ("B", 2), ("A", 3), ("D", 4), ("B", 3)]:
        W = make(typ, n)
        for w in W.elements:
            for s in W.left_descents(w):
                sw = W.left_mult(s, w)
                brute = {z for z in W.elements
                         if W.leq(w, z) and z != sw and W.leq(z, sw)
                         and W.length(W.left_mult(s, z)) < W.length(z)}
                got = set(W.strictly_below_with_descents(w, sw, s))
                if got != brute:
                    sb_bad.append((W.name, W.describe(w), s))
    check("EDGE strictly_below_with_descents == brute-force definition",
          not sb_bad, f"mismatches={len(sb_bad)} first={sb_bad[:3]}")

    # 8. LATENT DEFECT documented (not edited: frozen core).
    #    poly.mul(p, p) where p is a non-clean dict returns a fresh dict with
    #    the same entries -- harmless algebraically, but shows poly helpers
    #    do not preserve input identity; included so the parent can decide if
    #    aliasing matters downstream.
    p = {0: 1, 1: 2}
    finding("poly.mul/clean return new dict (no input aliasing)",
            f"poly.mul(p,p) is p = {poly.mul(p,p) is p}; "
            f"poly.clean(p) is p = {poly.clean(p) is p}; both equal p")


def main():
    print("=== S1 group-layer adversarial audit ===", flush=True)
    negative_control()
    edge_cases()
    audit("D", 4, 12, full_bruhat=True, bruhat_samples=0, word_samples=None)
    audit("B", 4, 16, full_bruhat=True, bruhat_samples=0, word_samples=None)
    audit("A", 5, 15, full_bruhat=False, bruhat_samples=4000, word_samples=None)
    if not QUICK:
        audit("A", 6, 21, full_bruhat=False, bruhat_samples=4000,
              word_samples=400)
    gt = ground_truth()
    fails = [r for r in RESULTS if not r[1]]
    print(f"\n=== SUMMARY: {len(RESULTS)-len(fails)} PASS / {len(fails)} FAIL ===")
    for name, ok, detail in fails:
        print(f"FAILED: {name}  {detail}")
    if gt:
        print("--- ground_truth.py output ---")
        print(gt)
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()