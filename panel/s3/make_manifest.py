#!/usr/bin/env python3
"""Seal the S3 panel's verified result numbers with koppa's FunctionToken.

Records exact integer results produced by panel/s3/probe.py. The token is a
deferred digest handle; resolution is on demand. No numeric deferral is used:
every number below is an exact, already-computed count.
"""
import hashlib, json, os, sys
sys.path.insert(0, "/home/sudoshellz/Documents/hecke-kl")
from koppa.token import seal, verify, tokenized_payload

HERE = os.path.dirname(os.path.abspath(__file__))

FACTS = {
    "seat": "S3 Hecke/module",
    "ground_truth": "PASS=21 FAIL=0",
    "D4_order": 192,
    "D4_bar_prime_all": {"checked": 192, "bad": 0},
    "D4_polynomials": {"count": 36864, "nonneg_bad": 0, "degree_bad": 0, "w0_bad": 0},
    "D4_R_identity_sample": {"checked": 30, "bad": 0},
    "solve_kl_bar": {"A3_pairs": 189, "A4_pairs": 3661, "B3_pairs": 799, "bad": 0},
    "algebra_A3": {"assoc_bad": 0, "bar2_bad": 0, "semilinear_bad": 0},
    "bar_convention": "semilinear_ring_homomorphism (bar(TuTv)=bar(Tu)bar(Tv))",
    "negative_controls": {"corrupt_Cprime_detected": True, "corrupt_P_detected": True},
    "mu_graph_components_all_one": ["A3", "A4", "B2", "G2"],
}

def main():
    payload = json.dumps(FACTS, sort_keys=True, separators=(",", ":")).encode()
    artifact, token, digest = seal(payload)
    manifest = {
        "schema": "hecke-kl-s3-panel-v1",
        "facts": FACTS,
        "koppa_handle": token.serialize(),
        "payload_sha256": digest.hex(),
    }
    path = os.path.join(HERE, "manifest.json")
    with open(path, "w") as f:
        json.dump(manifest, f, indent=2, sort_keys=True); f.write("\n")
    assert verify(tokenized_payload(payload, token), token)
    # tamper control
    bad = dict(FACTS); bad["D4_order"] = 191
    badpayload = json.dumps(bad, sort_keys=True, separators=(",", ":")).encode()
    assert not verify(tokenized_payload(badpayload, token), token)
    print("wrote", path, "sha256", digest.hex(), "token", token.serialize())
    print("tamper control: corrupted payload rejected")

if __name__ == "__main__":
    main()
