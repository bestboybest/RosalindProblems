from DnaToolkit import *

dna = list(readFasta(strongDir + "rosalind_orf.txt").values())[0]

prots = translateORFs(transcribe(dna))

for prot in set(prots):
    print(prot)
    