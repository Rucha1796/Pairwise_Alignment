from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices

# Reading sequences directly from files using SeqIO
seq1_records = list(SeqIO.parse("./seq1.txt", "fasta"))
seq2_records = list(SeqIO.parse("./seq2.txt", "fasta"))

# Extracting sequences
seq1 = str(seq1_records[0].seq)
seq2 = str(seq2_records[0].seq)

# Create an aligner object
aligner = PairwiseAligner()

# Load the BLOSUM62 matrix
blosum62 = substitution_matrices.load("BLOSUM62")
aligner.substitution_matrix = blosum62

best_score = float('-inf')  # Initialize to negative infinity
best_alignment = None
best_params = (None, None, None)

# Experimenting with parameters
for gap_open in [-5, -10, -15]:
    for gap_extend in [-0.5, -1, -1.5]:
        for penalize in [True, False]:
            aligner.open_gap_score = gap_open
            aligner.extend_gap_score = gap_extend
            
            # Adjust end gap behavior
            if penalize:
                aligner.target_end_gap_score = gap_open + gap_extend
                aligner.query_end_gap_score = gap_open + gap_extend
            else:
                aligner.target_end_gap_score = 0
                aligner.query_end_gap_score = 0
            
            alignment = aligner.align(seq1, seq2)
            
            # Assuming the first alignment is the best one (there could be multiple optimal alignments)
            score = alignment[0].score
            
            if score > best_score:
                best_score = score
                best_alignment = alignment[0]  # Take the first alignment
                best_params = (gap_open, gap_extend, penalize)

print(f"Best alignment score: {best_score}")
print(f"Best parameters: Gap open: {best_params[0]}, Gap extend: {best_params[1]}, Penalize end gaps: {best_params[2]}")
print(best_alignment)
