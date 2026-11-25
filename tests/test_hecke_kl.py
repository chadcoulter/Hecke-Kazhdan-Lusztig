"""Tests — every assertion traces to a value produced by the code under test."""
import os
import sys
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from hecke_kl.groups import make
from hecke_kl.kl import KLTable
from hecke_kl.hecke import HeckeAlgebra
from hecke_kl import poly


@pytest.mark.parametrize("typ,n,order", [("A", 3, 24), ("A", 4, 120), ("B", 3, 48),
                                         ("G2", None, 12), ("D", 4, 192)])
def test_group_order(typ, n, order):
    assert len(make(typ, n)) == order


def test_a2_all_kl_zero_one():
    W = make("A", 2)
    T = KLTable(W)
    vals = T.all_polynomials()
    assert all(p in ({}, {0: 1}) for p in vals.values())


def test_s4_wikipedia_anchors():
    W = make("A", 4)
    T = KLTable(W)
    y = W.from_word((1,))
    w = W.from_word((1, 0, 2, 1))
    assert poly.as_q_list(T.P(y, w)) == [1, 1]
    y = W.from_word((0, 2))
    w = W.from_word((0, 2, 1, 2, 0))
    assert poly.as_q_list(T.P(y, w)) == [1, 1]


def test_nonnegativity_and_degree_bound():
    for typ, n in [("A", 5), ("B", 3)]:
        W = make(typ, n)
        T = KLTable(W)
        for (y, w), p in T.all_polynomials().items():
            assert all(c >= 0 for c in p.values())
            if y != w and W.leq(y, w):
                assert p.get(0, 0) == 1
                assert poly.qdegree(p) <= (W.length(w) - W.length(y) - 1) // 2


def test_bar_invariance_and_quadratic():
    W = make("A", 4)
    T = KLTable(W)
    H = HeckeAlgebra(W, T)
    for w in W.elements:
        c = H.Cprime(w)
        assert H.equal(H.bar(c), c)
    for i in range(W.n):
        s = W.S[i]
        lhs = H.mul(H.T(s), H.T(s))
        rhs = H.add({W.identity: poly.q()},
                    H.scale(H.T(s), poly.sub(poly.q(), poly.one())))
        assert H.equal(lhs, rhs)


def test_inverse_and_r_identity():
    W = make("B", 3)
    T = KLTable(W)
    H = HeckeAlgebra(W, T)
    for w in W.elements:
        assert H.equal(H.mul(H.T(w), H.inverse_Tw(w)), H.T(W.identity))
    for y in W.elements:
        lhs = H.inverse_Tw(W.inverse(y))
        rhs = {}
        for x in W.elements:
            r = T.R(x, y)
            if r:
                rhs = H.add(rhs, {x: poly.shift(poly.bar(r), -2 * W.length(x))})
        assert H.equal(lhs, rhs)


def test_koppa_function_token_manifest():
    """Koppa's deferred digest is used to seal the verification manifest."""
    import hashlib
    from koppa.token import seal, verify, fixed_point_digest
    payload = b"hecke-kl verification vectors v1"
    artifact, token, digest = seal(payload)
    assert verify(artifact, token)
    assert token.resolve() == hashlib.sha256(payload).digest()
    # The handle is constant across digests: the framed bytes depend on the
    # payload but the token contribution never does.
    from koppa.token import tokenized_payload
    assert tokenized_payload(payload, token).endswith(b"\nKOPPA-FUNCTION-TOKEN-v1\n")
    from koppa.token import FunctionToken
    other_token = FunctionToken.from_digest(hashlib.sha256(b"tampered").digest())
    assert other_token.bytes_for_hashing() == token.bytes_for_hashing()
    # The handle is constant, so the framed bytes are identical for any token.
    # Verification must therefore depend on the resolved digest: a different
    # token's digest cannot validate the original payload.
    assert not verify(tokenized_payload(payload, other_token), other_token)
    assert verify(tokenized_payload(b"tampered", other_token), other_token)
