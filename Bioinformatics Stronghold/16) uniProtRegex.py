from DnaToolkit import *

with open(strongDir + "rosalind_mprt.txt") as f:
    proteins = [prot.strip() for prot in f.readlines()]

for prot in proteins:
    seq = list(getUniprot(prot).values())[0]
    positions = findMotifs(seq, r"N[^P][ST][^P]")
    if positions != []:
        print(prot)
        print(*positions)
    
