# Hecke–KL simulator — adversarial panel brief

## Frozen core (READ-ONLY — do not edit these files)
`/home/sudoshellz/Documents/hecke-kl/hecke_kl/{poly,coxeter,groups,kl,hecke,types}.py`

Contract: `CONTRACT.md` in the same directory. Ground truth: `ground_truth.py`
(21/21 PASS), `check_bruhat.py`, `check_rinvert.py` (50,656 pairs, 0 mismatches).
Koppa is installed (`import koppa`, v0.1.1) — use `koppa.token.FunctionToken`
for any self-hashing manifest; do NOT use numeric deferral for exact integer
coefficients (there is nothing to defer — no floats in the core).

## Parent's committed position (ATTACK IT — do not confirm)
"I believe the KL recursion, R-polynomials, Hecke multiplication, bar involution,
Bruhat order, and the published S4 1+q anchors are all correct for finite A/B/D/G2,
and that the only real gaps are (a) no KL cell/W-graph computation and (b) no
parabolic/affine types."

Each seat must either (i) FALSIFY a specific piece with a minimal reproducer
that produces a wrong value/crash, or (ii) report the exact test that killed its
best attempt. A seat that returns "the premise is wrong" without a run has
wasted its budget. No "cannot be done" sections.

## Disjoint seats
- **S1 — group layer**: attack Cartan/Tits conventions, BFS canonical words,
  Bruhat closure, reflection enumeration. Try to find a group where lengths,
  descents, Bruhat order, or longest element are wrong. Check D4, B4, A6.
  Required artifact: `panel/s1/report.md` + `panel/s1/probe.py` (runnable).
- **S2 — KL/R layer**: attack the KL recursion and R recursion with cases
  beyond the anchors. Check S5, B4, D4; verify mu/degree/positivity; try
  non-crystallographic G2 weight functions (should be rejected, not silently
  wrong). Required artifact: `panel/s2/report.md` + `panel/s2/probe.py`.
- **S3 — Hecke/module layer**: attack `mul`, `bar`, `Cprime`, `inverse_Tw`,
  associativity, the R identity, and the `solve_kl_bar` sympy fallback.
  Also build a KL cell / W-graph prototype (left cells via mu-graph) if time.
  Required artifact: `panel/s3/report.md` + `panel/s3/probe.py`.

## Rules
- Every claimed PASS/FAIL must come from an assertion on a value produced by code.
- Include a negative control: a deliberately corrupted variant must FAIL your probe.
- Report raw numbers. No hedging, no "honest limit", no framing.
- Do not edit the frozen core. Write only under `panel/<seat>/`.
- Re-run `python3 ground_truth.py` in your probe (must stay 21/21).
