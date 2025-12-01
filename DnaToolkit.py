Nucleotides = ["A", "C", "G", "T"]
Complements = {"A":"T", "C":"G", "G":"C", "T":"A"}

strongDir = "./Bioinformatics Stronghold/Input/"

def readFasta(file):
    with open(file) as f:
      first = 0
      finalSeqs = {}
      seqs = []
      for line in f:
        if line.startswith(">"):
            if (first != 0):
                finalSeqs[head] = "".join(seqs)
            head = line[1:].strip()
            seqs = []
            first = 1
        else:
            seqs.append(line.strip())
      finalSeqs[head] = "".join(seqs)
    return finalSeqs

def validateSeq(seq):
    for nuc in seq:
        if nuc.upper() not in Nucleotides:
            return False
    return seq.upper()

def countNucs(seq):
    return {nuc: seq.count(nuc) for nuc in Nucleotides}

def transcribe(seq):
    return seq.replace("T", "U")

def complement(seq):
    return "".join(Complements[nuc] for nuc in seq)

def revComplement(seq):
    return complement(seq)[::-1]

def gcContent(seq):
    return ((seq.count("C") + seq.count("G"))/len(seq)) * 100
    
