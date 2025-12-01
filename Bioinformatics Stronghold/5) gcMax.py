from DnaToolkit import *

file = strongDir + "rosalind_gc.txt"

seqs = readFasta(file)

gcCont = 0
for id, seq in seqs.items():
    if gcContent(seq) > gcCont:
        gcCont = gcContent(seq)
        maxId = id

print(maxId)
print(gcCont)