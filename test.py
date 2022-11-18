from permclass import *


# create a signed permutation
p = SignedPerm([3,2,1,4,5,6,7], [-1,-1,-1,1,1,1,1])
print(p)


# create the grid class of a single signed permutation
c = GridClass(p)
print(c.new_set)  # all compact permutations contained by the original set (original set = {p})


# create the grid class of permutations which can be sorted in at most 2 prefix reversals
c = GridClass.prefix_reversal(2)
print(c.original_set)  # the set of signed permutations we used to define this grid class


# find all elements of S_4 that can be sorted in at most 1 prefix reversal
c = GridClass.prefix_reversal(1)
print(c.members(4))


# find the unique compact signed permutation which is filled by p
compact_perm, vector = p.fills()
assert compact_perm.inflate(vector) == p
print(str(compact_perm)+" inflated by "+str(vector)+" is "+str(p))


print()


# the main purpose of this algorithm: generating functions and polynomials
c = GridClass.prefix_reversal(3)
print("Generating Function: "+str(c.genfcn()))
print("Enumerating Polynomial: "+str(c.polynomial()))
print("First 10 terms: "+str(c.sequence(10)))
