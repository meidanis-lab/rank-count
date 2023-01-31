Repository with code and data to compute the number of intermediate genomes and 
the number of scenarios to go from a genome to another by operations of rank 
distance 1 (cuts and joins) or 2 (double cut-and-joins).

Reference: J. P. Pereira Zanetti, L. P. Oliveira, J. Meidanis, L. Chindelevitch, 
"Counting Sorting Scenarios and Intermediate Genomes for the Rank Distance"; 
to appear in IEEE/ACM Transactions on Computational Biology and Bioinformatics.

The workflow used to produce Table 1 in the paper above is found in RankCount.R,
using the mainDriver() function. Please note that it relies on three packages:
combinat, igraph, and matrixStats, which need to be installed before using it.