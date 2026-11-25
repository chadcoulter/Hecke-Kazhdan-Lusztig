# S3 — Hecke/module layer: adversarial report

Seat: S3. Target: parent's committed position that "the KL recursion, R-polynomials,
Hecke multiplication, bar involution ... are all correct for finite A/B/D/G2".
Method: run assertions on values produced by code; negative controls confirm the
probes detect corruption. Probe: `panel/s3/probe.py` (runnable, EXIT=0).

## FINDING F1 (documentation/spec falsified; code is correct)

`CONTRACT.md` states: "Bar involution: bar(v)=v^{-1}; bar(T_w)=T_{w^{-1}}^{-1}
(**anti-involution**)." This characterization is WRONG. The standard KL bar
map is a **semilinear ring involution (homomorphism)**, not an anti-involution:

  bar(T_u T_v) = bar(T_u) bar(T_v)   [homomorphism — correct, KL79]
  bar(T_u T_v) = bar(T_v) bar(T_u)   [anti — FALSE]

Evidence (S4, u=s0, v=s1, noncommuting, from `_probe7.py`):
- `bar(T_u T_v)=bar(T_u)bar(T_v)`: **True**
- `bar(T_u T_v)=bar(T_v)bar(T_u)`: **False**
- bar(T_u T_v) expanded: T_{01}:v^-4 + T_0:(-v^-2+v^-4) + T_1:(-v^-2+v^-4) + T_e:(1-2v^-2+v^-4)
- anti-version differs ONLY in the top term: T_{10}:v^-4 instead of T_{01}:v^-4.

The separate map `*: T_w -> T_{w^{-1}}` (v fixed) is the algebra anti-involution;
the bar map is its composition with v -> v^-1, hence a *homomorphism*. This is
confirmed in the literature (arXiv:math/0210323, §2: bar is a ring involution
with bar(T_w)=T_{w^-1}^-1; `*` is the anti-automorphism). Random sparse tests on
A3: `bar(ab)=bar(a)bar(b)` held in 10/10 trials; `bar^2=id` held in 10/10.

Impact: the parent's phrase "bar involution ... correct" is only true because
the CODE (`HeckeAlgebra.bar`) implements the homomorphism, in contradiction to
the written CONTRACT. Any consumer who trusts the CONTRACT's "anti-involution"
line and writes `bar(ab)=bar(b)bar(a)` gets a wrong answer. Fix: correct the
CONTRACT sentence to "semilinear ring involution (homomorphism)", NOT the code.
This is a spec defect, not a code defect — the frozen core is untouched.

## No code failure found in the attacked layers (exact numbers)

- ground_truth.py: SUMMARY PASS=21 FAIL=0 (re-run inside probe).
- **D4 (order 192), not in ground_truth**: bar(C'_w)=C'_w for all 192 w, bad=0
  (60.8s). All 36,864 polynomials: nonneg_bad=0, degree_bad=0, w0_bad=0.
  R identity for all 30 w of length<=3: bad=0.
- `solve_kl_bar` independent sympy cross-check (repaired upstream): A3 189 pairs
  bad=0; A4 3661 pairs bad=0; B3 799 pairs bad=0 (S4/S5/B3 requirement met).
- Algebra laws, A3 sparse triples: associativity bad=0/10, bar^2=id bad=0/10,
  bar semilinearity bad=0/10.
- A5 full C' bar check died on runtime budget (720 elements x H.bar ~> 8 min);
  B4/D4-scale full checks similarly. This is the exact limit, not correctness.

## FINDING F2 (W-graph prototype): naive mu-graph is degenerate

The BRIEF's suggested "left cells via mu-graph" (undirected edge when
mu(x,y)!=0) does NOT separate cells. Computed connected components:
- A3 (order 24): 60 edges, 1 component of size 24
- A4 (order 120): 482 edges, 1 component of size 120
- B2 (order 8): 12 edges, 1 component of size 8
- G2 (order 12): 20 edges, 1 component of size 12
Every group collapses to a single cell. This is the expected consequence of
including mu(e,w)!=0 edges (e connects to everything with odd top KL coeff).
A correct left-cell algorithm needs the directed STAR operation / left-descent
refinement (Lusztig's a-function / Duflo order), not the raw mu graph. So the
parent's gap (b) is real AND the suggested trivial fix is insufficient.

## Negative controls (both FIRE)

- Corrupt C' by adding 1 to a non-identity T_0 coefficient: bar(C')=C' becomes
  **False** (detected). (Adding only to the T_e coefficient does NOT break
  bar-invariance at w=s0s1 because v^-2+1 is self-bar for that element — a
  useful caution on choosing the corruption site.)
- KLTable subclass adding +1 to every constant term: bar(C')!=C' becomes
  **True** (detected).
- koppa FunctionToken manifest tamper: corrupted payload is rejected.

## Artifacts

- `panel/s3/probe.py` — runnable (EXIT=0), all assertions pass.
- `panel/s3/report.md` — this file.
- `panel/s3/make_manifest.py` — koppa FunctionToken self-hashing manifest.
- `panel/s3/manifest.json` — handle `koppa-ft1:5eebbae4...b551c5f585d563`,
  sha256 of facts payload `5eebbae4...b551c5f585d563`. Exact integers, no
  numeric deferral.
- Scratch diagnostics `panel/s3/_probe*.py` (kept for reproducibility).
