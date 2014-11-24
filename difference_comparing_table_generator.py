import pandas as pd
from Bio import SeqIO



ref_genome_fasta_16S = '/home/toomas/git/plot_generation/plotter/16S_new.fasta'
ref_genome_fasta_23S = '/home/toomas/git/plot_generation/plotter/23S_new.fasta'
#position = [j for j in range(1, 1543)]
position_16S = []
position_23S = []
nucleotide_16S = []
nucleotide_23S = []


for fasta in SeqIO.parse(ref_genome_fasta_16S, "fasta"):
    pass

for index in range(1-20, 1582 + 1-20):
    position_16S.append(index)
    nucleotide_16S.append(str(fasta.seq[index-4+19:index+3+19]))

for fasta in SeqIO.parse(ref_genome_fasta_23S, "fasta"):
        pass

for index in range(1-20, 2944 + 1-20):
     position_23S.append(index)
     nucleotide_23S.append(str(fasta.seq[index-4+19:index+3+19]))


     
df_16S = pd.DataFrame(index=position_16S, columns=["nucleotide", "MazF_none/MG_log_none", "MqsR_none/MG_log_none",
                                               "MG_stats_none/MG_log_none", "d3_none/MG_log_none", "MG_log_PNK/MG_log_none", 
                                               "MazF_PNK/MazF_none", "MqsR_PNK/MqsR_none", "MG_stats_PNK/MG_stats_none",
                                               "d3_PNK/d3_none"])

df_16S["nucleotide"] = nucleotide_16S
print df_16S

df_23S = pd.DataFrame(index=position_23S, columns=["nucleotide", "MazF_none/MG_log_none", "MqsR_none/MG_log_none",
                                               "MG_stats_none/MG_log_none", "d3_none/MG_log_none", "MG_log_PNK/MG_log_none", 
                                               "MazF_PNK/MazF_none", "MqsR_PNK/MqsR_none", "MG_stats_PNK/MG_stats_none",
                                               "d3_PNK/d3_none"])

df_23S["nucleotide"] = nucleotide_23S
print df_23S


df_16S.to_csv("16S_table.csv", index=True ,header=True)
df_23S.to_csv("23S_table.csv", index=True ,header=True)