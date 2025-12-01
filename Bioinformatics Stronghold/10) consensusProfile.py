from DnaToolkit import *

file = strongDir + "rosalind_cons.txt"

seqs = readFasta(file)

avSeq = consensus(seqs)

seqs = list(seqs.values())
a, c, g, t = [], [], [], []
for i in range(len(seqs[0])):
    tempSeq = [seqs[j][i] for j in range(len(seqs))]
    a.append(tempSeq.count('A'))
    c.append(tempSeq.count('C'))
    g.append(tempSeq.count('G'))
    t.append(tempSeq.count('T'))

print(avSeq)
print("A: ", *a)
print("C: ", *c)
print("G: ", *g)
print("T: ", *t)


