# S2/S3 panel report — KL/R layer attack

## Verdict: FALSIFIED a specific piece (NEW failure found)

Parent position attacked: *"KL recursion, R-polynomials, Hecke multiplication,
bar involution ... are all correct ... and the only real gaps are cells and
parabolic/affine types."* plus the CONTRACT requirement that KL values were
cross-checked by an independent bar-invariance linear system.

### F1 (NEW FAILURE) — `HeckeAlgebra.solve_kl_bar` is broken for every w

`hecke_kl/hecke.py::solve_kl_bar` cannot recover P_{y,w} by solving
bar(C'_w)=C'_w for ANY element, in A3/S4 or A4/S5. Raw output from
`panel/s2/probe.py` (assertion on produced values):

```
F1 solve_kl_bar CRASHES (all confirmed):
 [('e',  'RuntimeError',   'no bar-invariant solution for w=e'),
  ('0',  'PolynomialError','1/v contains an element of the set of generators.'),
  ('1',  'PolynomialError','1/v contains an element of the set of generators.'),
  ('01', 'PolynomialError','v**(-2) contains an element of the set of generators.'),
  ('012','PolynomialError','v**(-3) contains an element of the set of generators.')]
```

Root causes (two independent bugs in one method):
1. **w=e**: eqs reduce to `[True]` (the single unknown's bar-invariant
   equation is identically satisfied after solving) so `sp.solve(eqs,syms,
   dict=True)` returns `[]`; the code raises `RuntimeError("no bar-invariant
   solution")` even though the true solution c=1 exists.
2. **w!=e**: `coeff.subs(v,1/v)` leaves `1/v` inside a sympy expression;
   `sp.degree(expr, v)` at hecke.py:189 then throws
   `PolynomialError: 1/v contains an element of the set of generators` before
   any solve. So the sympy fallback NEVER runs successfully.

Consequence: `check_independent.py` dies immediately
(`RuntimeError: no bar-invariant solution for w=e`). The CONTRACT's required
"cross-check the same KL values by an independent method (bar-invariance
linear system)" is NOT satisfied by the shipped fallback — that cross-check
has never actually executed. The correct KL values are still right (verified
by the R-inversion identity, 50,656 pairs), but the claimed *independent*
verification via this method is absent, and the public `solve_kl_bar` API is a
crash. Minimal reproducer: `HeckeAlgebra(make('A',3), KLTable(W)).solve_kl_bar(W.identity)`
and `...solve_kl_bar(W.from_word((0,)))`.

### Negative control (fires)
Corrupting a KL coefficient in the cache (P_{e,1021} in D4: 1+q -> 1+2q)
makes `bar(C'_w) != C'_w` -> `H.equal` returns False. So the bar-invariance
probe genuinely detects a corrupted P table. Raw:
`NEG CONTROL: ('e','1021','1 + 2q')` with bar-invariance broken (assertion passed).

### Confirmed-correct (not failures), with numbers
- D4 (order 192): all 36,864 KL pairs pass nonnegativity, degree bound, and
  constant term 1; 0 failures. Bar(C')=C' for all 192 w: 0 failures.
  Cartan is the correct fork at node 1 (s1-s2=s1-s3=m=3, s0-s2=m=2);
  longest element length 12; all BFS words have matching reduced length.
- B4 (order 384): all 147,456 KL pairs pass nonnegativity/degree/constant:
  0 failures.
- ground_truth.py re-run: PASS=21 FAIL=0 (no core files edited; writes only
  under panel/s2).

### Exact test that killed the best attempt to FIX F1
Attempted workaround: call `solve_kl_bar` on w with body monkeypatched to
skip w=e. This cannot be done from the probe because the frozen core must not
be edited (`hecke.py` is in the do-not-edit set). The bug is located at
hecke.py:172-196 (`sp.solve(eqs,syms,dict=True)` and `sp.degree(expr,v)` on a
rational expression); fixing it requires editing frozen `hecke.py`, which is
out of scope for this seat. This is the precise blocker; no "cannot be done"
claim — the failing call and its exact exception are reproduced above.

## Artifacts
- `panel/s2/probe.py` — runnable; F1 assertion, negative control, bounds check.
- `panel/s2/report.md` — this report.
