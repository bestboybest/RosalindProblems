Nucleotides = ["A", "C", "G", "T"]
Complements = {"A":"T", "C":"G", "G":"C", "T":"A"}

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

