Nucleotides = ["A", "C", "G", "T"]

def validateSeq(seq):
    for nuc in seq:
        if nuc.upper() not in Nucleotides:
            return False
    return seq.upper()

def len(seq):
    return len(seq)

def countNucs(seq):
    count = {}
    for nuc in Nucleotides:
        count[nuc] = seq.count(nuc)
    return count
