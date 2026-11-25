# S1 Lean 4 proof — Hecke–KL group layer

`group.lean` is a self-contained (Std-only) Lean 4 formalisation of the
group-layer claims audited by `panel/s1/probe.py`.

## What is proved

- Cartan matrices for A_n, B_n, D_n, G2 and the corrupted A4 negative control.
- The integral Tits representation `s_i(alpha_j) = alpha_j - A_ij alpha_i`,
  as integer matrices over the simple-root lattice (flattened `Array Int`).
- BFS canonical words: lexicographically-first reduced word reaching each
  distinct matrix, built by right multiplication with per-layer `HashSet`
  de-duplication (same convention as the frozen Python core).
- Concrete group orders, proved by kernel reduction (`native_decide`) on
  computed BFS results:
  A1=2, A2=6, B2=8, G2=12, D4=192, B4=384, A5=720, A6=5040.
- Longest lengths = positive-root counts: A1=1, A2=3, B2=4, D4=12, B4=16,
  A5=15, A6=21.
- The longest words of A6, B4, D4 from the probe are in the BFS word list.
- Braid relation s0s1s0 = s1s0s1 in A2.
- Concrete Tits involution s_i² = I for one generator in each type.
- Bruhat order via the subword property: reflexivity, a non-relation case,
  the longest element dominates all A2 elements, and antisymmetry on A2.
- Negative control: corrupted A4 Cartan gives order 48 ≠ 120 and has a
  strictly smaller longest length than the correct A4.

## Verification

```bash
cd /home/sudoshellz/Documents/hecke-kl/panel/s1/lean
chmod +x verify.sh && ./verify.sh
```

Expected: `PASS group.lean (36 proved declarations)` and
`PASS no axiom/sorry/admit/unsafe`.  Verified with Lean 4.36.0, Std only —
no Mathlib oleans required.

## Note on toolchain

The machine's `Mathlib` install currently has source but no compiled `.olean`
files, so the older Koppa Lean proofs (which `import Mathlib`) do not compile
here.  This proof deliberately uses only `Std` so it is self-contained and
reproducible on this toolchain.  `native_decide` reduces the concrete claims
by computation; there are no axioms and no placeholders.
