#!/usr/bin/env python3
"""CLI for the Hecke–KL simulator."""
from __future__ import annotations

import argparse
import sys
from .groups import make
from .kl import KLTable
from .hecke import HeckeAlgebra
from . import poly


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="hecke-kl", description="Hecke–KL algebra simulator")
    ap.add_argument("type", help="Coxeter type: A4, B3, D4, G2, I2-4, ...")
    ap.add_argument("--word-y", default="", help="y as generator indices, e.g. 02")
    ap.add_argument("--word-w", default="", help="w as generator indices, e.g. 02120")
    ap.add_argument("--all", action="store_true", help="list all nonzero distinct P polynomials")
    ap.add_argument("--check", action="store_true", help="run the 21 ground-truth anchors")
    args = ap.parse_args(argv)

    if args.check:
        import ground_truth
        return 0

    t = args.type.upper().replace("-", "")
    if t.startswith("I2"):
        W = make("I2", int(t[2:]))
    elif t.startswith("G2"):
        W = make("G2")
    elif len(t) == 1:
        ap.error("specify rank, e.g. A4 not A")
        return 2
    else:
        W = make(t[0], int(t[1:]))
    T = KLTable(W)
    print(f"W = {W.name}  order={W.order}  longest={W.describe(W.longest)} "
          f"(length {W.length(W.longest)})")

    if args.all:
        seen = {}
        for (y, w), p in T.all_polynomials().items():
            if p and y != w:
                seen.setdefault(poly.fmt(p), []).append((W.describe(y), W.describe(w)))
        for pstr, pairs in sorted(seen.items(), key=lambda kv: -len(kv[1])):
            print(f"  P = {pstr}   ({len(pairs)} pairs) e.g. {pairs[:3]}")
        return 0

    if not args.word_y or not args.word_w:
        ap.error("provide --word-y and --word-w (or --all)")
        return 2
    y = W.from_word([int(c) for c in args.word_y])
    w = W.from_word([int(c) for c in args.word_w])
    print(f"y = {W.describe(y)} (l={W.length(y)}), w = {W.describe(w)} (l={W.length(w)})")
    print(f"y <= w in Bruhat order: {W.leq(y, w)}")
    p = T.P(y, w)
    print(f"P_{{y,w}}(q) = {poly.fmt(p)}")
    print(f"mu(y,w) = {T.mu(y, w)}")
    H = HeckeAlgebra(W, T)
    c = H.Cprime(w)
    print(f"C'_w has {len(c)} terms; bar-invariant: {H.equal(H.bar(c), c)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
