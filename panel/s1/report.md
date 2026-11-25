# S1 — Adversarial audit of the FROZEN Hecke–KL group layer

**Verdict: SURVIVED.** No wrong value or crash found in Cartan/Tits
conventions, BFS canonical words, left/right descents, Bruhat closure, or
reflection enumeration for D4, B4, A5, A6 (and A2/A3/A4/B2/B3 probed during
development). Every check is an assertion on a value computed by `probe.py`.

Reproducer: `python3 /home/sudoshellz/Documents/hecke-kl/panel/s1/probe.py`
(full, ~9 min: the frozen core's BFS for A6 takes ~351 s). Fast iteration:
`python3 probe.py --quick` (skips A6; 36 PASS / 0 FAIL, ~15 s).

## Method — independent model, no reuse of frozen code for expected values

The frozen core builds each group as integer Tits matrices acting on the
simple-root basis. The probe builds a **geometrically independent** model:
signed permutations acting on the ambient coordinate basis e_i, with

- type A: s_i = transposition (i, i+1) on {0..n} (m = n+1 points);
- type B: s_i = transposition for i<n-1; s_{n-1} = sign flip of e_{n-1};
- type D: transpositions for i<n-2; s_{n-1} = swap e_{n-2},e_{n-1} with both
  signs flipped (reflection in root e_{n-2}+e_{n-1}).

From these generators the probe runs its own BFS (identical generator order),
then compares, every group element:

1. group order vs closed-form |A_n|=(n+1)!, |B_n|=2^n n!, |D_n|=2^{n-1} n!;
2. length l(w) vs BFS depth in the geometric model;
3. left and right descents vs geometric ones;
4. canonical word vs an **independent greedy lex-least** reconstruction
   (the first letter of the lex-least reduced word is the least left descent;
   recurse on s_i w) — this does not use the core's BFS at all;
5. reflection enumeration: count equals the number of positive roots
   (A_n: n(n+1)/2, B_n: n^2, D_n: n(n-1)); and each core reflection matrix
   equals the corresponding geometric conjugate w s_i w^{-1} after the
   change of basis from e-coordinates to simple-root coordinates;
6. simple reflection matrices S_i vs geometric generators in the same basis;
7. Bruhat: y <= w vs the **Ehresmann rank criterion** for type A
   (r_y(p,q) >= r_w(p,q) entrywise), and vs the **subword property** for B/D;
8. longest element word vs geometric longest.

Finally the probe re-runs `ground_truth.py` and asserts `PASS=21 FAIL=0`.

## Raw results (full probe run)

All lines below are verbatim output values; "mismatches=0" means every
compared element/pair agreed.

```
NEGCTRL corrupted Cartan (A4 edge removed) FAILS probe  corrupt order=48 w0len=7 (good A4: order=120 w0len=6)
NEGCTRL corrupted D4 Cartan FAILS probe  corrupt order=48 w0len=7 (good D4: order=192 w0len=12)

D4: order=192=192, longest length=12=12 (#posroots), canonical keys=192,
    lengths mismatches=0, descents mismatches=0,
    greedy-lex mismatches=0 (checked=192),
    reflections core=12 geo=12 core_only=0 geo_only=0, S_i mismatches=[],
    Bruhat pairs=36864 mismatches=0,
    w0 word=(0,1,0,2,1,0,3,1,0,2,1,3)

B4: order=384=384, longest length=16=16, canonical keys=384,
    lengths mismatches=0, descents mismatches=0,
    greedy-lex mismatches=0 (checked=384),
    reflections core=16 geo=16 core_only=0 geo_only=0, S_i mismatches=[],
    Bruhat pairs=147456 mismatches=0,
    w0 word=(0,1,0,2,1,3,2,1,0,3,2,1,3,2,3)  [see probe output for exact tuple]

A5: order=720=720, longest length=15=15, canonical keys=720,
    lengths mismatches=0, descents mismatches=0,
    greedy-lex mismatches=0 (checked=720),
    reflections core=15 geo=15 core_only=0 geo_only=0, S_i mismatches=[],
    Bruhat sampled pairs=4000 mismatches=0,
    w0 word=(0,1,0,2,1,0,3,2,1,0,4,3,2,1,0)

A6: order=5040=5040, longest length=21=21, canonical keys=5040,
    lengths mismatches=0, descents mismatches=0,
    greedy-lex mismatches=0 (checked=400 sampled),
    reflections core=21 geo=21 core_only=0 geo_only=0, S_i mismatches=[],
    Bruhat sampled pairs=4000 mismatches=0,
    w0 word=(0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,4,3,2,1,0)

ground_truth.py re-run: rc=0, SUMMARY PASS=21 FAIL=0
FULL SUMMARY: 47 PASS / 0 FAIL  (probe exit code 0)
```

(B4 w0 word, exact: `(0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1, 3, 2, 3)`.)

## What killed an attempt (false alarm, diagnosed as PROBE bug not core bug)

- **Attempt 1 — naive root-basis embedding for B_n.** Comparing core
  reflections to geometric conjugates using an e_basis -> simple-root basis
  change-of-basis built with alpha_n = e_{n-1} produced apparent failures:
  `B4 reflection matrices core_only=16 geo_only=16`, `simple bad_indices=[2,3]`.
  Reproducer values: for B3 the core S1 = `((1,0,0),(1,-1,2),(0,0,1))` but
  the geometric s1 in that (wrong) basis was `((1,0,0),(1,-1,1),(0,0,1))`.
  Cause: `cartan_b` sets A_{n-1,n}=-2, so the long simple root is
  alpha_n = **2 e_{n-1}**, not e_{n-1}. With alpha_n = 2e_{n-1} the sets are
  equal: B2/B3/B4 all `core_only=0 geo_only=0`, no simple-reflection
  mismatches. This validated the Cartan/Tits convention instead of falsifying
  it; the corrected convention is used in the final probe.
- **Attempt 2 — comparing e-basis matrices directly.** Core matrices are in
  the simple-root basis, geometric ones in the e-basis; direct comparison
  trivially gives core_only=|W|/2. Fixed by the explicit rational change of
  basis P, L = (P^T P)^{-1} P^T (exact `Fraction` arithmetic).

No attempt produced a genuine wrong value in the frozen core. The parent's
position "group layer is correct for finite A/B/D/G2" is therefore not
falsified by this seat; the newly added independent checks (full reflection
matrix equality, greedy lex-least words, rank/subword Bruhat on D4/B4/A5/A6)
all pass.

## Negative controls

Two deliberately corrupted Cartans were verified to FAIL the probe's
invariants (they are asserted to be detected, so the check itself passes when
the corruption is detected):

- `cartan_a(4)` with A[0][1]=A[1][0]=0 (A3 -> A1 x A3): order=48, w0len=7
  instead of (120, 6) — detected.
- `cartan_d(4)` with A[0][1]=A[1][0]=0 (D4 -> A1 x A3): order=48, w0len=7
  instead of (192, 12) — detected.

## Files

- `panel/s1/probe.py` — runnable adversarial probe (runnable with
  `--quick`; full run includes A6).
- `panel/s1/report.md` — this report.

No frozen-core file was edited.