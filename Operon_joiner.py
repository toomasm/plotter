import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


input_16S = {"rrsA": '/home/toomas/git/plot_generation/plotter/sequences/rrsA_extra_20.fasta',
          "rrsB": '/home/toomas/git/plot_generation/plotter/sequences/rrsB_extra_20.fasta',
          "rrsC": '/home/toomas/git/plot_generation/plotter/sequences/rrsC_extra_20.fasta',
          "rrsD": '/home/toomas/git/plot_generation/plotter/sequences/rrsD_extra_20.fasta',
          "rrsE": '/home/toomas/git/plot_generation/plotter/sequences/rrsE_extra_20.fasta',
          "rrsG": '/home/toomas/git/plot_generation/plotter/sequences/rrsG_extra_20.fasta',
          "rrsH": '/home/toomas/git/plot_generation/plotter/sequences/rrsH_extra_20.fasta'}


input_23S =  {"rrlA": '/home/toomas/git/plot_generation/plotter/sequences/rrlA_extended_20_compensated.fasta',
          "rrlB": '/home/toomas/git/plot_generation/plotter/sequences/rrlB_extra_20.fasta',
          "rrlC": '/home/toomas/git/plot_generation/plotter/sequences/rrlC_extra_20.fasta',
          "rrlD": '/home/toomas/git/plot_generation/plotter/sequences/rrlD_extra_20.fasta',
          "rrlE": '/home/toomas/git/plot_generation/plotter/sequences/rrlE_extra_20.fasta',
          "rrlG": '/home/toomas/git/plot_generation/plotter/sequences/rrlG_extra_20.fasta',
          "rrlH": '/home/toomas/git/plot_generation/plotter/sequences/rrlH_extra_20.fasta'}


IUPAC_dict = {"A": "A", "C": "C", "G": "G", "T": "T", 
              "R": "AG", "Y": "CT", "S": "CG", "W": "AT",
              "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
              "H": "ACT", "V": "ACG", "N": "ACGT"}


seq_dict_16S = {}
seq_dict_23S = {}
seq_list_16S = ["rrsA", "rrsB", "rrsC", "rrsD", "rrsE", "rrsG", "rrsH"]
seq_list_23S = ["rrlA", "rrlB", "rrlC", "rrlD", "rrlE", "rrlG", "rrlH"]

new_16S = ""
new_23S = ""

#parse 16S rRNA-s
for seq in seq_list_16S:
    for seq_dict_16S[seq] in SeqIO.parse(input_16S[seq], "fasta"):
        pass
#parse 23S rRNA-s
for seq in seq_list_23S:
    for seq_dict_23S[seq] in SeqIO.parse(input_23S[seq], "fasta"):
        pass

for index in range(1, 1542+1+40):
    nucleotide_set = set([])
    for seq in seq_list_16S:
        nucleotide_set.add((seq_dict_16S[seq]).seq[index-1])
        
    nucleotide_string = ''.join(sorted((str(i) for i in nucleotide_set)))
    
    for key, value in IUPAC_dict.items():
        if value == nucleotide_string:
            new_16S += key

for index in range(1, 2905+40):
    nucleotide_set = set([])
    for seq in seq_list_23S:
        nucleotide_set.add((seq_dict_23S[seq]).seq[index-1])
    
    nucleotide_string = ''.join(sorted((str(i) for i in nucleotide_set)))
    
    for key, value in IUPAC_dict.items():
        if value == nucleotide_string:
            new_23S += key
            
print new_16S
print len(new_16S)
print new_23S
print len(new_23S)

rec_16S = SeqRecord(Seq(new_16S, IUPAC.IUPACAmbiguousDNA),
                 id="E_coli_16S_7_ribosomes_joined",
                 description="E. coli 16S ribosome seq composed of all 7 operons (IUPAC)")

output_handle = open("16S_new.fasta", "w")
SeqIO.write(rec_16S, output_handle, "fasta")
output_handle.close()

rec_23S = SeqRecord(Seq(new_23S, IUPAC.IUPACAmbiguousDNA),
                 id="E_coli_23S_7_ribosomes_joined",
                 description="E. coli 23S ribosome seq composed of all 7 operons (IUPAC)")

output_handle = open("23S_new.fasta", "w")
SeqIO.write(rec_23S, output_handle, "fasta")
output_handle.close()