#!/usr/bin/env python3
import time, sys
sys.path.insert(0, '/home/sudoshellz/Documents/hecke-kl')
from hecke_kl.groups import make
W = make('A', 5)
# Check subword enumeration cost for pairs
T = W.word(W.longest)
print('longest', W.describe(W.longest), 'len', len(T))
from itertools import combinations
from collections import Counter
cnt = Counter()
for y in W.elements:
    ly = W.length(y)
    cnt[ly] += 1
print('length dist', sorted(cnt.items()))
# time one worst-case leq
y = W.from_word(tuple(range(5)) + (0,1,0,2) )
t = time.time()
ok = W.leq(y, W.longest)
print('one leq', ok, time.time()-t)
