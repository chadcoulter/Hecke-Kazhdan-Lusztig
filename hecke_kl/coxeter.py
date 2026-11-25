"""Exact finite Coxeter groups via the integral reflection (Tits) representation.

For a crystallographic Cartan matrix A with A_ii = 2 and A_ij = 2(α_i,α_j)/(α_i,α_i),
the simple reflection s_i acts on the simple-root coordinate lattice by

    s_i(α_j) = α_j - A_ij α_i .

Thus s_i is an integer matrix S_i with column j equal to e_j - A_ij e_i.
Group elements are multiplication tables of these matrices; the representation is
faithful for finite Coxeter groups. Canonical reduced words are the
lexicographically least reduced word found by BFS (which visits length-k words in
lexicographic order).
"""
from __future__ import annotations
from itertools import combinations
from collections import deque


class CoxeterGroup:
    def __init__(self, name: str, cartan: list[list[int]]):
        self.name = name
        self.cartan = tuple(tuple(int(x) for x in row) for row in cartan)
        self.n = len(self.cartan)
        self._build()
        self._bruhat_cache: dict[tuple[tuple[int, ...], tuple[int, ...]], bool] = {}

    # ---------------------------------------------------------------- building
    def _simple_matrix(self, i: int) -> tuple[tuple[int, ...], ...]:
        n = self.n
        m = [[0] * n for _ in range(n)]
        for j in range(n):
            for k in range(n):
                m[k][j] = 1 if k == j else 0
            m[i][j] -= self.cartan[i][j]
        return tuple(tuple(row) for row in m)

    @staticmethod
    def _mat_mul(a, b):
        n = len(a)
        return tuple(
            tuple(sum(a[i][k] * b[k][j] for k in range(n)) for j in range(n))
            for i in range(n)
        )

    def _build(self):
        S = [self._simple_matrix(i) for i in range(self.n)]
        self.S = S
        n = self.n
        ident = tuple(tuple(1 if i == j else 0 for j in range(n)) for i in range(n))
        self.identity = ident
        # BFS over right multiplication; words visited in lex order per layer.
        word_of = {ident: ()}
        elements = [ident]
        by_length: list[list] = [[ident]]
        queue = deque([(ident, ())])
        while queue:
            m, w = queue.popleft()
            for i in range(self.n):
                nm = self._mat_mul(m, S[i])
                if nm not in word_of:
                    nw = w + (i,)
                    word_of[nm] = nw
                    elements.append(nm)
                    while len(by_length) <= len(nw):
                        by_length.append([])
                    by_length[len(nw)].append(nm)
                    queue.append((nm, nw))
        self.elements = elements
        self.word_of = word_of
        self.by_length = by_length
        self.length_of = {m: len(w) for m, w in word_of.items()}
        self.index_of = {m: i for i, m in enumerate(elements)}
        self.longest = max(elements, key=lambda m: self.length_of[m])
        self.order = len(elements)
        self._bruhat_cache: dict[tuple[tuple[int, ...], tuple[int, ...]], bool] = {}
        self._build_bruhat()

    def _build_bruhat(self):
        """Bruhat order from covering relations, not exponential subwords.

        For a reflection t = w s_i w^{-1}, u < u·t is a cover exactly when
        length(u·t) = length(u)+1. The Bruhat order is the transitive closure
        of these covers. We compute down-sets by a length-descending bitset DP:
            down[w] = {w} ∪ union over covers c of w of down[c].
        """
        # all reflections: conjugates w s_i w^{-1}
        refl = set()
        for w in self.elements:
            winv = self.inverse(w)
            for i in range(self.n):
                refl.add(self._mat_mul(self._mat_mul(w, self.S[i]), winv))
        self.reflections = list(refl)
        idx = self.index_of
        # down-set bitsets, computed by INCREASING length (covers are shorter).
        down = [0] * self.order
        for layer in self.by_length:
            for w in layer:
                bits = 1 << idx[w]
                lw = self.length_of[w]
                for t in self.reflections:
                    c = self._mat_mul(w, t)  # c = w t
                    if c in idx and self.length_of[c] == lw - 1:
                        bits |= down[idx[c]]
                down[idx[w]] = bits
        self._down = down

    def bruhat_leq(self, y, w) -> bool:
        y, w = self._check(y), self._check(w)
        return bool((self._down[self.index_of[w]] >> self.index_of[y]) & 1)

    # ------------------------------------------------------------------ access
    def __len__(self):
        return self.order

    def __iter__(self):
        return iter(self.elements)

    def word(self, w):
        return self.word_of[self._check(w)]

    def length(self, w):
        return self.length_of[self._check(w)]

    def _check(self, w):
        if isinstance(w, tuple) and len(w) == self.n and all(
            isinstance(row, tuple) for row in w
        ):
            if w in self.word_of:
                return w
        raise ValueError(f"not a group element of {self.name}: {w!r}")

    def from_word(self, word) -> tuple:
        """Product of simple reflections given by generator indices."""
        m = self.identity
        for i in word:
            if not (0 <= i < self.n):
                raise ValueError(f"generator index out of range: {i}")
            m = self._mat_mul(m, self.S[i])
        return m

    def mul(self, u, v):
        return self._mat_mul(self._check(u), self._check(v))

    def inverse(self, w):
        w = self._check(w)
        # Inverse matrix found by lookup (group is finite and represented faithfully).
        for g in self.elements:
            if self._mat_mul(w, g) == self.identity:
                return g
        raise RuntimeError("inverse not found")

    def right_descents(self, w):
        w = self._check(w)
        return tuple(i for i in range(self.n)
                     if self.length(self.mul(w, self.S[i])) < self.length(w))

    def left_descents(self, w):
        w = self._check(w)
        return tuple(i for i in range(self.n)
                     if self.length(self.mul(self.S[i], w)) < self.length(w))

    def right_mult(self, w, i):
        return self.mul(w, self.S[i])

    def left_mult(self, i, w):
        return self.mul(self.S[i], w)

    def leq(self, y, w) -> bool:
        """Bruhat order via covering-relation transitive closure (bitset DP)."""
        return self.bruhat_leq(y, w)

    def strictly_below_with_descents(self, y, sw, s=None):
        """z with y <= z < sw and s z < z (used by the KL recursion)."""
        sw = self._check(sw)
        out = []
        for k in range(self.length(y), self.length(sw)):
            for z in self.by_length[k]:
                if self.leq(y, z):
                    if s is None:
                        out.append(z)
                    elif self.length(self.left_mult(s, z)) < self.length(z):
                        out.append(z)
        return out

    def elements_between(self, y, w):
        return [g for g in self.elements if self.leq(y, g) and self.leq(g, w)]

    def is_longest(self, w):
        return self._check(w) == self.longest

    def describe(self, w):
        return "e" if w == self.identity else "".join(str(i) for i in self.word(w))
