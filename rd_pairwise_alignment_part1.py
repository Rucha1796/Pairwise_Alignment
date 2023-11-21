from Bio import pairwise2
from Bio.Align import substitution_matrices
from Bio import SeqIO

# Reading sequences directly from files using SeqIO
seq1_records = list(SeqIO.parse("./seq1.txt", "fasta"))
seq2_records = list(SeqIO.parse("./seq2.txt", "fasta"))

# Extracting the sequence from the first (and only) record.
seq1 = str(seq1_records[0].seq)
seq2 = str(seq2_records[0].seq)

# Loading the BLOSUM62 matrix
blosum62 = substitution_matrices.load("BLOSUM62")
    
#Global alignment using BLOSUM62 matrix
global_alignments = pairwise2.align.globalds(seq1, seq2, blosum62, -10, -0.5)
# Local alignment using BLOSUM62 matrix
local_alignments = pairwise2.align.localds(seq1, seq2, blosum62, -10, -0.5)

print("\nGlobal Alignments:\n")
for alignment in global_alignments:
    print(pairwise2.format_alignment(*alignment))

print("\nLocal Alignments:\n")
for alignment in local_alignments:
    print(pairwise2.format_alignment(*alignment))

