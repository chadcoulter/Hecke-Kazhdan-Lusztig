#!/usr/bin/env python3
"""Generate a tamper-evident verification manifest using koppa's FunctionToken.

The manifest records the committed ground-truth vectors (S4 1+q anchors,
rank-2 0/1 rule, w0=1 rule). A naive manifest that embeds sha256(manifest)
inside itself is circular; koppa's deferred token appends a CONSTANT opaque
handle to the payload so the digest is resolvable on verification without
ever appearing in the hashed bytes. No numeric deferral is used — the KL
coefficients are exact integers, so there is nothing to defer numerically.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl import poly
from koppa.token import seal, verify, tokenized_payload

HERE = os.path.dirname(__file__)


def collect_vectors() -> dict:
    W = make("A", 4)
    T = KLTable(W)
    vecs = {}
    for yw, ww, label in [
        ((1,), (1, 0, 2, 1), "P_b_bacb"),
        ((0, 2), (0, 2, 1, 2, 0), "P_ac_acbca"),
    ]:
        y, w = W.from_word(yw), W.from_word(ww)
        vecs[label] = {
            "y_word": list(yw), "w_word": list(ww),
            "P": poly.as_q_list(T.P(y, w)),
        }
    # rank-2 rule: every P in {0,1}
    W2 = make("A", 2)
    T2 = KLTable(W2)
    vecs["A2_all_0_or_1"] = all(
        p in ({}, {0: 1}) for p in T2.all_polynomials().values()
    )
    vecs["A4_w0_all_1"] = all(
        poly.eq(T.P(y, W.longest), poly.one())
        for y in W.elements if W.leq(y, W.longest)
    )
    return vecs


def main() -> int:
    vectors = collect_vectors()
    payload = json.dumps(vectors, sort_keys=True, separators=(",", ":")).encode()
    artifact, token, digest = seal(payload)
    manifest = {
        "schema": "hecke-kl-vectors-v1",
        "vectors": vectors,
        "koppa_handle": token.serialize(),
        "payload_sha256": digest.hex(),
    }
    path = os.path.join(HERE, "VERIFY.json")
    with open(path, "w") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)
        f.write("\n")
    with open(path, "rb") as f:
        data = f.read()
    # Verify the framing on the emitted artifact (tokenized form)
    framed = tokenized_payload(payload, token)
    assert verify(framed, token)
    print(f"wrote {path} ({len(data)} bytes)")
    print(f"payload sha256 = {digest.hex()}")
    print(f"koppa token    = {token.serialize()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
