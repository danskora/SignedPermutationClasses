from newpolypermclass import *


flips = 7
C = PolyClass.prefix_reversal(flips)


#print(C.genfcn())
print(C.sequence(10))
print(C.polynomial())


# S = PegPerm.sort_pr(4)
# for s in S:
#     print(s)