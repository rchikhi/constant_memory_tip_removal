 This software remove tips and short isolated vertices from a de bruijn graph that resides in a file, in constant memory
 input: set of nodes (fasta)
 output: filtered set of nodes (fasta)

 definitions:

 - a tip is a node such that:
   * it has no in-neighbors xor no out-neighbors, and
   * mean k-mer coverage < coverage_threshold, and 
   * node length < max_tip_length
 - an isolated vertex is one which as no in-neighbor and no in-neighbor and has length < max_tip_length
 - max_tip_length is set to 2k, in accordance with Velvet 
 - coverage_threshold is a user-defined parameter

 For convenience, the default memory usage is set to be exactly (2*number of distinct k-mers) bits

 known bugs:
 - if a node is larger than (genome_size/4) nucleotides then this software will use more than (2*genome_size) memory

 license: CECILL 
 (uses parts of Minia)

 author: Rayan Chikhi
 date: 09/2013
