Nucleotides = ["A", "C", "G", "T"]
Complements = {"A":"T", "C":"G", "G":"C", "T":"A"}
Aminos = {"UUU": "F", "UUC": "F", "UUG": "L", "UUA": "L", "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*", "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W", "CUU": "L", "CUA": "L", "CUC": "L", "CUG": "L", "CCC": "P", "CCU": "P", "CCA": "P", "CCG": "P", "CAU": "H", "CAC": "H", "CAG": "Q", "CAA": "Q", "CGU": "R", "CGG":"R", "CGA":"R", "CGC":"R", "AUU": "I", "AUC": "I", "AUA": "I", "AUG" : "M", "ACU":"T", "ACC":"T", "ACG":"T", "ACA":"T", "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K", "AGU": "S", "AGC":"S", "AGA":"R", "AGG":"R", "GUU":"V", "GUA":"V", "GUG":"V", "GUC":"V", "GCU":"A", "GCC":"A", "GCG":"A", "GCA":"A", "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E", "GGU":"G", "GGG":"G", "GGA":"G", "GGC":"G"}

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
    return "".join(Complements[nuc] for nuc in seq)

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

def translate(seq):
    if (len(seq) % 3 != 0):
        return False
    return "".join([Aminos[seq[i:i+3]] for i in range(0, len(seq), 3)])

def proteins(seq):
    return [protein for protein in seq.split("*") if protein != ""]

#functionality to find from string/regex 
def findMotifs(seq, motif):
    positions = []
    for match in re.compile(motif).finditer(seq):
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