"""Type A3 (S4) conformance fixtures taken verbatim from published examples.

Source: Wikipedia "Kazhdan-Lusztig polynomial", Examples section:
  S = {a, b, c} with a and c commuting (the A3 Coxeter diagram a-b-c):
      P_{b, bacb} = 1 + q
      P_{ac, acbca} = 1 + q
Generator convention here: a=s0, b=s1, c=s2.
Also: rank<=2 groups all KL polynomials are 0/1; l-difference in {0,1,2}
gives P=1; P_{y,w0}=1 for finite W; all coefficients nonnegative.
"""
from __future__ import annotations

# words as generator-index tuples, a=s0,b=s1,c=s2
A3_P_b_bacb = (1,)
A3_WORD_bacb = (1, 0, 2, 1)
A3_P_ac_acbca = (1,)
A3_WORD_ac = (0, 2)
A3_WORD_acbca = (0, 2, 1, 2, 0)

A3_NONTRIVIAL = [
    # (y_word, w_word, expected coeffs of P_{y,w} in q, source)
    (A3_WORD_bacb[:1], A3_WORD_bacb, [1, 1], "Wikipedia: P_{b,bacb}=1+q"),
    (A3_WORD_ac, A3_WORD_acbca, [1, 1], "Wikipedia: P_{ac,acbca}=1+q"),
]

# Sanity anchors independent of any KL implementation
ANCHORS = {
    "rank_le_2_all_0_or_1": True,
    "length_diff_0_1_2_P_is_1": True,
    "longest_element_all_P_are_1": True,
    "coefficients_nonnegative": True,
    "S4_nontrivial_count_at_least_2": 2,
}
