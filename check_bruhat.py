#!/usr/bin/env python3
"""Cross-check the fast Bruhat closure against the slow subword criterion.

Independent method: subword property (y <= w iff some subword of a reduced
word for w multiplies to y). If both agree on every pair in small groups,
the bitset DP is validated on the same inputs before we use it for A5.
"""
import sys
from itertools import combinations
sys.path.insert(0, "/home/sudoshellz/Documents/hecke-kl")
from hecke_kl.groups import make

for typ, n in [("A", 3), ("A", 4), ("B", 3), ("I2", 6)]:
    W = make(typ, n)
    bad = []
    for w in W.elements:
        ww = W.word(w)
        for y in W.elements:
            ly = W.length(y)
            sub = any(W.from_word(c) == y for c in combinations(ww, ly))
            fast = W.leq(y, w)
            if sub != fast:
                bad.append((W.describe(y), W.describe(w), sub, fast))
    print(f"{W.name}: {len(bad)} mismatches vs subword criterion", bad[:3], flush=True)
    if bad:
        sys.exit(1)
print("PASS: Bruhat closure == subword criterion on all small groups")
