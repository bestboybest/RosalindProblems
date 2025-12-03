import itertools
from DnaToolkit import *

n = 6

nums = list(range(1, n + 1))
perms = list(itertools.permutations(nums))

f = open(strongDir + "19solution.txt", "w")
f.write(str(len(perms)) + "\n")
for p in perms:
    f.write(" ".join([str(i) for i in p]))
    f.write("\n")