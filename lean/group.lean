import Std

set_option maxHeartbeats 4000000
set_option maxRecDepth 30000

-- Flattened row-major matrices: entry (i,j) = data.getD (i*n+j) 0.
abbrev MatKey := Array Int

namespace M

def id (n : Nat) : MatKey :=
  Array.ofFn (fun k : Fin (n*n) => if k.val / n = k.val % n then 1 else 0)

def entry (n : Nat) (M : MatKey) (i j : Nat) : Int := M.getD (i*n + j) 0

def mk (n : Nat) (f : Fin n -> Fin n -> Int) : MatKey :=
  Array.ofFn (fun k : Fin (n*n) =>
    let i : Fin n := ⟨k.val / n, Nat.div_lt_of_lt_mul k.isLt⟩
    let j : Fin n := ⟨k.val % n, Nat.mod_lt k.val (Nat.pos_of_mul_pos_left (Nat.lt_of_le_of_lt (Nat.zero_le k.val) k.isLt))⟩
    f i j)

def mul (n : Nat) (a b : MatKey) : MatKey :=
  Array.ofFn (fun k : Fin (n*n) =>
    let i := k.val / n; let j := k.val % n
    (List.finRange n).foldl (fun acc t => acc + entry n a i t.val * entry n b t.val j) 0)

end M

def cartanA (n : Nat) (i j : Fin n) : Int :=
  if i = j then 2
  else if (j.val + 1 = i.val) ∨ (i.val + 1 = j.val) then -1 else 0

def cartanB (n : Nat) (i j : Fin n) : Int :=
  if i = j then 2
  else if i.val + 1 = j.val then (if j.val + 1 = n then -2 else -1)
  else if j.val + 1 = i.val then -1
  else 0

def cartanD (n : Nat) (i j : Fin n) : Int :=
  if i = j then 2
  else if (i.val < n - 3 ∧ j.val = i.val + 1) ∨ (j.val < n - 3 ∧ i.val = j.val + 1) then -1
  else if (i.val = n - 3 ∧ j.val = n - 2) ∨ (j.val = n - 3 ∧ i.val = n - 2) then -1
  else if (i.val = n - 3 ∧ j.val = n - 1) ∨ (j.val = n - 3 ∧ i.val = n - 1) then -1
  else 0

def cartanG2 : Fin 2 -> Fin 2 -> Int := fun i j =>
  if i = j then 2 else if i.val = 0 ∧ j.val = 1 then -1 else -3

def cartanA4Corrupt : Fin 4 -> Fin 4 -> Int := fun i j =>
  if i = j then 2
  else if (i.val = 0 ∧ j.val = 1) ∨ (i.val = 1 ∧ j.val = 0) then 0
  else if (j.val + 1 = i.val) ∨ (i.val + 1 = j.val) then -1 else 0

def simpleMat {n} (A : Fin n -> Fin n -> Int) (i : Fin n) : MatKey :=
  M.mk n (fun j k => (if k = j then 1 else 0) - A i j * (if k = i then 1 else 0))

def fromWord {n} (A : Fin n -> Fin n -> Int) (w : List (Fin n)) : MatKey :=
  w.foldl (fun acc i => M.mul n acc (simpleMat A i)) (M.id n)

def bfsWords {n} (A : Fin n -> Fin n -> Int) (fuel : Nat) : List (List (Fin n)) :=
  let rec extend (len : Nat) (cur : List (List (Fin n))) (seen : Std.HashSet MatKey) :=
    if _ : len >= fuel then [] else
      let step (st : List (List (Fin n)) × Std.HashSet MatKey) (w : List (Fin n)) :=
        let m := fromWord A w
        if st.2.contains m then st else (st.1 ++ [w], st.2.insert m)
      let (fresh, seen') := (cur.flatMap (fun w => (List.finRange n).map (fun i => w ++ [i]))).foldl step ([], seen)
      fresh ++ extend (len+1) fresh seen'
  let initial : List (List (Fin n)) := [[]]
  let seen0 : Std.HashSet MatKey := {M.id n}
  initial ++ extend 1 initial seen0

def bfsCount {n} (A : Fin n -> Fin n -> Int) (fuel : Nat) : Nat := (bfsWords A fuel).length

def numPosRootsA (n : Nat) : Nat := n * (n + 1) / 2
def numPosRootsB (n : Nat) : Nat := n * n
def numPosRootsD (n : Nat) : Nat := n * (n - 1)

-- Group orders and positive-root counts.
example : bfsCount (cartanA 1) 5 = 2 := by native_decide
example : bfsCount (cartanA 2) 6 = 6 := by native_decide
example : bfsCount (cartanB 2) 7 = 8 := by native_decide
example : bfsCount cartanG2 9 = 12 := by native_decide
example : bfsCount (cartanD 4) 14 = 192 := by native_decide
example : bfsCount (cartanB 4) 18 = 384 := by native_decide
example : bfsCount (cartanA 5) 17 = 720 := by native_decide
example : bfsCount (cartanA 6) 23 = 5040 := by native_decide

example : numPosRootsA 1 = 1 := by native_decide
example : numPosRootsA 2 = 3 := by native_decide
example : numPosRootsB 2 = 4 := by native_decide
example : numPosRootsD 4 = 12 := by native_decide
example : numPosRootsB 4 = 16 := by native_decide
example : numPosRootsA 5 = 15 := by native_decide
example : numPosRootsA 6 = 21 := by native_decide

-- Negative control.
example : bfsCount cartanA4Corrupt 11 = 48 := by native_decide
example : bfsCount (cartanA 4) 13 = 120 := by native_decide
theorem corrupted_is_not_correct : bfsCount cartanA4Corrupt 11 ≠ bfsCount (cartanA 4) 13 := by native_decide
theorem corrupted_has_wrong_longest : bfsCount cartanA4Corrupt 11 < bfsCount (cartanA 4) 13 := by native_decide

-- Braid relations and longest words.
theorem A2_braid : fromWord (cartanA 2) [0,1,0] = fromWord (cartanA 2) [1,0,1] := by native_decide

def w0_A6 : List (Fin 6) := [0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,4,3,2,1,0]
example : w0_A6.length = 21 := by native_decide
example : w0_A6 ∈ bfsWords (cartanA 6) 23 := by native_decide

def w0_B4 : List (Fin 4) := [0,1,0,2,1,0,3,2,1,0,3,2,1,3,2,3]
example : w0_B4.length = 16 := by native_decide
example : w0_B4 ∈ bfsWords (cartanB 4) 18 := by native_decide

def w0_D4 : List (Fin 4) := [0,1,0,2,1,0,3,1,0,2,1,3]
example : w0_D4.length = 12 := by native_decide
example : w0_D4 ∈ bfsWords (cartanD 4) 14 := by native_decide

-- Bruhat order via subword property.  sublistsOf implements the list of all
-- subsequences (delete zero or more elements), which is sufficient because
-- the subword criterion uses subsequences of a reduced word.
def sublistsOf {α} : List α -> List (List α)
  | [] => [[]]
  | a::as => let rest := sublistsOf as; rest ++ rest.map (fun l => a :: l)

def bruhatLeq {n} (A : Fin n -> Fin n -> Int) (y : MatKey) (w : List (Fin n)) : Bool :=
  (sublistsOf w).any (fun sub => fromWord A sub == y)

example : bruhatLeq (cartanA 2) (fromWord (cartanA 2) []) [] = true := by native_decide
example : bruhatLeq (cartanA 2) (fromWord (cartanA 2) [0]) [0] = true := by native_decide
example : bruhatLeq (cartanA 2) (fromWord (cartanA 2) [1]) [0] = false := by native_decide
example : bruhatLeq (cartanA 2) (fromWord (cartanA 2) [0]) [0,1] = true := by native_decide
example : (bfsWords (cartanA 2) 6).all (fun w => bruhatLeq (cartanA 2) (fromWord (cartanA 2) w) [0,1,0]) = true := by native_decide
example : ((bfsWords (cartanA 2) 6).all (fun y => (bfsWords (cartanA 2) 6).all (fun w =>
    !(bruhatLeq (cartanA 2) (fromWord (cartanA 2) y) w) ||
    !(bruhatLeq (cartanA 2) (fromWord (cartanA 2) w) y) ||
    (fromWord (cartanA 2) y == fromWord (cartanA 2) w)))) = true := by native_decide

-- Concrete Tits involution checks on each Cartan type.
example : M.mul 2 (simpleMat (cartanA 2) 0) (simpleMat (cartanA 2) 0) = M.id 2 := by native_decide
example : M.mul 2 (simpleMat (cartanB 2) 1) (simpleMat (cartanB 2) 1) = M.id 2 := by native_decide
example : M.mul 4 (simpleMat (cartanD 4) 3) (simpleMat (cartanD 4) 3) = M.id 4 := by native_decide
example : M.mul 6 (simpleMat (cartanA 6) 5) (simpleMat (cartanA 6) 5) = M.id 6 := by native_decide
