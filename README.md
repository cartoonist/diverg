DiVerG
======
Scalable Distance Index for Validation of Paired-End Alignments in Sequence Graphs
----------------------------------------------------------------------------------

This is an implementation of DiVerG method introduced in:

> A. Ghaffaari, T. Marschall, and A. Schoenhuth.
> "DiVerG: Scalable Distance Index for Validation of Paired-End Alignments in Sequence Graphs", 2024
>
> [manuscript](./doc/diverg-manuscript.pdf)

### Introduction

DiVerG is an indexing scheme based on [PairG](https://github.com/ParBLiSS/PairG)
that efficiently determines whether there is a path between any two loci in a
sequence graph that falls within a particular, statistically well-motivated
range $[d_1, d_2]$. This problem is motivated by the need in paired-end read
mapping to sequence graph to verify distance constraints between candidate
alignments.

It introduces a compact data structure for representing Boolean sparse matrices,
as well as a fast and scalable algorithm for computing matrix-matrix
multiplication and addition using the compressed representation on CUDA and
OpenMP backends.

It overcomes the limitations of PairG by exploiting the
extensive potential for improvements in terms of scalability and space
efficiency. As a consequence, it can process substantially larger datasets, such
as whole human genomes, which are unmanageable by PairG.
DiVerG offers faster index construction time and consistently faster query time
with gains proportional to the size of the underlying compact data structure.

The constructed index that can be used for efficient distance querying between
two loci for constraint validation.
