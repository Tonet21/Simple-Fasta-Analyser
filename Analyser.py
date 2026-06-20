from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqUtils import gc_fraction
import pandas as pd
import matplotlib.pyplot as plt

def n_cp(seq):
    A = seq.count("A")
    C =seq.count("C")
    G = seq.count("G")
    T = seq.count("T")
    AP = round(A / len(seq) * 100, 2)
    CP = round(C / len(seq) * 100, 2)
    GP = round(G / len(seq) * 100, 2)
    TP = round(T / len(seq) * 100, 2)
    return {"A": (A, AP),
            "C": (C, CP),
            "G": (G, GP),
            "T": (T, TP)}


sequences = []

with open ('Sequences.fasta') as handle:
    for s_id, seq in SimpleFastaParser(handle):
        seq = Seq(seq)
        sequences.append((s_id,seq))

stats = [] 

for s_id, seq in sequences:
    length = len(seq)
    nucle = n_cp(seq)
    GC = gc_fraction(seq)
    RC = seq.reverse_complement()
    rna = seq.transcribe()
    stats.append({"ID": s_id,
                  "length": length,
                  "nucle": nucle,
                  "GC": GC,
                  "RC": RC,
                  "rna": rna})
    

summary_table = pd.DataFrame(stats, columns = ["ID","length", "GC"])

def template(item): 
    ## The items are the elements from the stats list. For instance, stats[0] which would be the first dictionary.
    
    s_id, l, cp, gc, rc, rna = item["ID"], item["length"], item["nucle"], item["GC"], item["RC"], item["rna"]
    a = cp["A"]
    c = cp["C"]
    g = cp["G"]
    t = cp["T"]
    return f""" 
    {s_id}
    --------------    
    Length: {l} bp
    --------------
    A: {a[0]}, Percentage of A: {a[1]}%
    C: {c[0]}, Percentage of C: {c[1]}%
    G: {g[0]}, Percentage of G: {g[1]}%
    T: {t[0]}, Percentage of T: {t[1]}%
    --------------
    GC content: {round(gc * 100,2)}%
    --------------
    Reverse complement: {rc}
    --------------
    RNA: {rna}

"""
print(summary_table)

summary_table["GC"].hist()
plt.show()

summary_table["length"].hist()
plt.show()

