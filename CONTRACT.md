# Hecke–KL simulator: shared contract (read-only for panel seats)

## Scope
Finite Coxeter groups (types A_n via permutations; B_n/C_n, D_n, G_2, I_2(m) via
Cartan/reflection matrices), their Hecke algebras H(W,S) over Z[v^{±1}] with
q = v^2, R-polynomials, Kazhdan–Lusztig polynomials P_{y,w}(q), the KL basis
C'_w, bar involution, Bruhat order, and KL cells (W-graphs where feasible).

## Canonical conventions (all modules MUST agree)
- Generators indexed 0..n-1; simple reflections s_i.
- Group element = canonical reduced word (tuple of generator indices), lexicographically
  smallest among reduced expressions; multiplication = word concatenation reduced via
  the group signature (permutation for A, integer/sympy reflection matrix otherwise).
- Hecke relation: T_s T_s = q + (q-1) T_s, equivalently (T_s+1)(T_s-q)=0.
  If l(sw)>l(w): T_s T_w = T_{sw}; else T_s T_w = (q-1)T_w + q T_{sw}.
  Right-handed version: T_w T_s = T_{ws} if l(ws)>l(w), else (q-1)T_w + q T_{ws}.
- Polynomials stored as dict[int,int]: v-exponent -> integer coefficient.
  q = v^2. KL polynomials live in Z[q] (even v-exponents only).
- mu(x,y) = coefficient of v^{l(y)-l(x)-1} in P_{x,y}; 0 unless x<y and
  l(y)-l(x) is odd.
- KL recursion (KL79 (2.2.c)), for y<=w, s in left descents of w (sw<w):
  if sy > y: P_{y,w} = P_{sy,w};
  else: P_{y,w} = P_{sy,sw} + q P_{y,sw}
          - sum_{z: y<=z<sw, sz<z} mu(z,sw) q^{(l(w)-l(z))/2} P_{y,z}.
  Base: P_{w,w}=1; P_{y,w}=0 if y not<= w.
- R recursion (Wikipedia/KL79), for s with sy<y:
  R_{y,y}=1; if sx<x: R_{x,y}=R_{sx,sy}; if sx>x: R_{x,y}=(q-1)R_{sx,y}+q R_{sx,sy}.
  Identity to verify: T_{y^{-1}}^{-1} = sum_x bar(R_{x,y}) q^{-l(x)} T_x.
- KL basis: C'_w = v^{-l(w)} sum_{y<=w} P_{y,w}(v^2) T_y; must satisfy bar(C'_w)=C'_w.
- Bar involution: bar(v)=v^{-1}; bar(T_w)=T_{w^{-1}}^{-1} (anti-involution).

## Ground truth anchors (must hold, no exceptions)
A1. Rank <= 2 Coxeter groups: P_{y,w}=1 for all y<=w, 0 otherwise.
A2. l(w)-l(y) in {0,1,2} => P_{y,w}=1 when y<=w.
A3. w = longest element of finite W => P_{y,w}=1 for all y<=w.
A4. S4 (A3) with generators a=s0,b=s1,c=s2 (a,c commute):
    P_{b,bacb} = 1+q and P_{ac,acbca} = 1+q. Nontrivial: NOT all ones.
A5. Coefficients of KL polynomials for finite W are nonnegative integers.
A6. Hecke: associativity, quadratic relation, braid relations, T_e identity.
A7. Degree bound: for y<w, deg_q P_{y,w} <= (l(w)-l(y)-1)/2, constant term 1.
A8. Bar invariance of C'_w for every w in every supported group.
A9. R identity above for every y.

## Verification requirements
- Every claimed pass must come from an assertion on a value produced by the code.
- Negative controls must fire: a deliberately broken recursion / wrong q-normalization
  must FAIL the invariants; a corrupted P table must be detected.
- Cross-check the same KL values by an independent method (bar-invariance linear
  system or geometric/Soergel-style reconstruction) on at least S4, S5, B3.
- Clean-venv install test from outside the source tree.

## Deliverable layout
~/Documents/hecke-kl/
  hecke_kl/{__init__,poly,coxeter,hecke,kl,cells,cli,types}.py
  tests/ (pytest)
  examples/demo.py
  BRIEF.md (this file is CONTRACT.md; seats read both)
  VERIFY.md, verify_all.sh
