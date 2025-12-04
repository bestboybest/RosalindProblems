Nucleotides = ["A", "C", "G", "T"]
ComplementsD = {"A":"T", "C":"G", "G":"C", "T":"A"}
ComplementsR = {"A":"U", "C":"G", "G":"C", "U":"A"}
Aminos = {"UUU": "F", "UUC": "F", "UUG": "L", "UUA": "L", "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*", "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W", "CUU": "L", "CUA": "L", "CUC": "L", "CUG": "L", "CCC": "P", "CCU": "P", "CCA": "P", "CCG": "P", "CAU": "H", "CAC": "H", "CAG": "Q", "CAA": "Q", "CGU": "R", "CGG":"R", "CGA":"R", "CGC":"R", "AUU": "I", "AUC": "I", "AUA": "I", "AUG" : "M", "ACU":"T", "ACC":"T", "ACG":"T", "ACA":"T", "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K", "AGU": "S", "AGC":"S", "AGA":"R", "AGG":"R", "GUU":"V", "GUA":"V", "GUG":"V", "GUC":"V", "GCU":"A", "GCC":"A", "GCG":"A", "GCA":"A", "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E", "GGU":"G", "GGG":"G", "GGA":"G", "GGC":"G"}
Residues={"A":71.03711,"C":103.00919,"D":115.02694,"E":129.04259,"F":147.06841,"G":57.02146,"H":137.05891,"I":113.08406,"K":128.09496,"L":113.08406,"M":131.04049,"N":114.04293,"P":97.05276,"Q":128.05858,"R":156.10111,"S":87.03203,"T":101.04768,"V":99.06841,"W":186.07931,"Y":163.06333}
H20 = 18.01056

strongDir = "./Bioinformatics Stronghold/Input/"

import requests
import re

#functionality to read from fasta string or file
def readFasta(file):
    if not file.startswith(">"):
        with open(file, "r") as f:
            fasta = f.readlines()
    else:
        fasta = file.splitlines()
    first = 0
    finalSeqs = {}
    seqs = []
    for line in fasta:
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
    if "T" in seq:
        return "".join(ComplementsD[nuc] for nuc in seq)
    else:
        return "".join(ComplementsR[nuc] for nuc in seq)

def revComplement(seq):
    return complement(seq)[::-1]

def gcContent(seq):
    return ((seq.count("C") + seq.count("G"))/len(seq)) * 100
    
def hamming(seq1, seq2):
    if (len(seq1) != len(seq2)):
        return False
    distance = 0
    for i in range(len(seq1)):
        if (seq1[i] != seq2[i]):
            distance += 1
    return distance

#Works for any dna string from its start
def translate(seq):
    seq = seq[0:len(seq) - (len(seq) % 3)]
    return "".join([Aminos[seq[i:i+3]] for i in range(0, len(seq), 3)])

#Finding with regex lookahead
def proteins(seq):
    return [prot[:-1] for prot in re.findall(r"(?=(M[^*]*\*))", seq)]

#functionality to find from string/regex 
def findMotifs(seq, motif): 
    positions = []
    for match in re.finditer(motif, seq):
        positions.append(match.start() + 1)
    return positions

def consensus(seqs):
    avSeq = []
    seqs = list(seqs.values())
    for i in range(len(seqs[0])):
        tempSeq = [seqs[j][i] for j in range(len(seqs))]
        avSeq.append(max(set(tempSeq), key = tempSeq.count))
    return "".join(avSeq)

def longestMotif(seqs):
    seqs = list(seqs.values())
    short = min(seqs, key = len)
    x = len(short)
    while True:
        if x == 0:
            return ""
        for i in range(len(short) - x + 1):
            check = True
            for seq in seqs:
                if short[i: i + x] not in seq:
                    check = False
                    break
            if check:
                return short[i:i + x]
        x -= 1

def getUniprot(Uniprot_id):
    Uniprot_id = re.match(r"^[^_]+", Uniprot_id).group()
    url = f'http://www.uniprot.org/uniprot/{Uniprot_id}.fasta'
    res = requests.get(url)
    if res.status_code != 200:
        return None
    fasta = res.text
    return readFasta(fasta)

def translatePar(seq):
    proteins1 = proteins(translate(seq))
    proteins2 = proteins(translate(seq[1:]))
    proteins3 = proteins(translate(seq[2:]))
    return proteins1+proteins2+proteins3

def translateORFs(seq):
    proteinss1 = translatePar(seq)
    proteinss2 = translatePar(revComplement(seq))
    return proteinss1 + proteinss2

def protMass(seq, full = 0):
    base = sum(Residues[amino] for amino in seq)
    return (base + H20) if full else base

def revPal(seq, length):
    positions = []
    for i in range(len(seq)-length+1):
        if seq[i:i+length] == revComplement(seq[i:i+length]):
            positions.append(i+1)
    return positions