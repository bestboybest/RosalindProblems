from DnaToolkit import *

file = strongDir + "rosalind_grph.txt"

seqs = readFasta(file)

for idi, vali in seqs.items():
    for idj, valj in seqs.items():
        if idi != idj:
            if vali[-3:] == valj[:3]:
                print(idi, idj)
