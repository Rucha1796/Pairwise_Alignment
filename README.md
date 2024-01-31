# rd_pairwise_alignment
Utilizes the Biopython package for conducting pairwise sequence alignment, aiming to assess the similarity or relatedness between sequences from two distinct organisms.
This process is essential for understanding the evolutionary relationships between sequences from different organisms. The script's key features are:
1) Reading Sequences: It reads two sequences from files ('seq1.txt' and 'seq2.txt') in FASTA format using SeqIO.parse.
Extracting Sequences: The script extracts the sequences from these files and converts them into strings for alignment.
2) Pairwise Alignment Setup: It initializes a PairwiseAligner object from Biopython, which is used for the alignment process.
3) Substitution Matrix: The script uses the BLOSUM62 matrix, a common substitution matrix for protein alignments, loaded via substitution_matrices.load.
4) Parameter Experimentation: The script iterates over different values of gap opening (gap_open) and gap extension (gap_extend) penalties, as well as toggling the penalization of gaps at the ends of alignments. These are parameters that influence how the alignment is scored and optimized.
5) Alignment and Scoring: For each set of parameters, the script performs an alignment and scores it. The assumption is that the first alignment in the returned set is the best one (in cases of multiple optimal alignments).
6) Identifying the Best Alignment: It compares the scores from different parameter configurations to find the best alignment score and the corresponding parameters.
7) Output: Finally, the script prints the best alignment score, the parameters that led to this score, and the alignment itself.
