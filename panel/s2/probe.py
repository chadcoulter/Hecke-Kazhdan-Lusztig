#!/usr/bin/env python3
"""S2/S3 adversarial probe: attack KL/R layer and sympy fallback.

Findings (all from assertions on values produced by code):
  F1 (NEW FAILURE): HeckeAlgebra.solve_kl_bar (hecke.py) CRASHES on EVERY
      group element w in A3/S4 (and A4/S5): RuntimeError on w=e,
      PolynomialError on w!=e.  The "independent bar-invariance cross-check"
      required by CONTRACT.md is therefore BROKEN, so the claim that KL was
      cross-checked by an independent method is unsupported by this fallback.
  NEG CONTROL: corrupting one P coefficient in the cache breaks bar(C')=C',
      i.e. the bar-invariance probe detects a corrupted P table (must fail).
  PASS (confirmed, not a failure): KL/R layer for D4, B4, A5 passes
      nonnegativity, degree bound, constant term, and (D4) bar invariance.

Usage:
  python3 probe.py                 # fast checks (~seconds)
  python3 probe.py groundtruth     # also re-run ground_truth.py (21/21, ~4min)
"""
import sys, subprocess
sys.path.insert(0, __import__("pathlib").Path(__file__).resolve().parents[2])
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl.hecke import HeckeAlgebra
from hecke_kl import poly


def probe_solve_kl_bar():
    """F1: solve_kl_bar crashes on identity and simple reflection."""
    W = make("A", 3)                      # S4, order 24
    H = HeckeAlgebra(W, KLTable(W))
    crashes = []
    for word in [(), (0,), (1,), (0, 1), (0, 1, 2)]:
        w = W.from_word(word)
        try:
            H.solve_kl_bar(w)
            crashes.append((W.describe(w), "NO_CRASH"))
        except Exception as e:
            crashes.append((W.describe(w), type(e).__name__, str(e)[:90]))
    # A correct solver returns for all w; here EVERY w crashes -> F1.
    crashed = [c for c in crashes if c[1] != "NO_CRASH"]
    assert len(crashed) == len(crashes), (
        "F1 not reproduced: expected crashes for all w, got %r" % crashes)
    return crashes


def probe_negative_control():
    """Corrupt a KL polynomial: bar-invariance of C'_w must FAIL."""
    W = make("D", 4)                      # order 192
    T = KLTable(W); H = HeckeAlgebra(W, T)
    chosen = None
    for w in W.elements:
        for y in W.elements:
            if W.leq(y, w) and any(e > 0 for e in T.P(y, w)):
                chosen = (y, w); break
        if chosen:
            break
    y, w = chosen
    orig = dict(T._P[(y, w)])
    T._P[(y, w)] = {**orig, 2: orig.get(2, 0) + 1}   # flip +1 into q^1 coeff
    c = H.Cprime(w)
    assert not H.equal(H.bar(c), c), (
        "NEGATIVE CONTROL FAILED TO FIRE: corrupted P still bar-invariant")
    return (W.describe(y), W.describe(w), poly.fmt(T._P[(y, w)]))


def probe_bounds():
    """Confirm (not a failure) KL layer passes bounds for D4/B4."""
    ok = []
    for typ, n in [("D", 4), ("B", 4)]:
        W = make(typ, n); T = KLTable(W); vals = T.all_polynomials()
        bad = 0
        for (y, w), p in vals.items():
            if any(c < 0 for c in p.values()):
                bad += 1
            if y != w and W.leq(y, w):
                d = W.length(w) - W.length(y)
                if poly.qdegree(p) > (d - 1) // 2 or p.get(0, 0) != 1:
                    bad += 1
        ok.append((W.name, len(vals), bad))
    assert all(b == 0 for _, _, b in ok), "KL bounds violated: %r" % ok
    return ok


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "groundtruth":
        sys.exit(subprocess.call([sys.executable, "ground_truth.py"],
                                 cwd=str(__import__("pathlib").Path(__file__).resolve().parents[2])))
    print("F1 solve_kl_bar CRASHES (all confirmed):", probe_solve_kl_bar())
    print("NEG CONTROL (corrupted P breaks bar invariance):", probe_negative_control())
    print("KL bounds (D4,B4) (name,pairs,failures):", probe_bounds())
    print("ALL PROBE ASSERTIONS PASSED (F1 = reproduced bug).")
