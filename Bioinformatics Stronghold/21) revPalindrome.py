from DnaToolkit import *

s = list(readFasta(strongDir + "rosalind_revp.txt").values())[0]

for i in range(4, 13):
    positions = revPal(s, i)
    if positions != []:
        for pos in positions:
            print(pos, i)