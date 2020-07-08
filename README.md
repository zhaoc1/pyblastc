# pyblastc

## BLAST command

The following blast command reproduce the same results with online blastn interface:

```
blastn  -query {input.contigs} \
        -db {params.db} \
        -outfmt '6 qseqid sseqid pident qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore' \
        -num_threads {threads} -evalue 1e-10 \
        -max_target_seqs 5000 \
        -out {output}
```

- When building custom database, e.g. the sequences of a gene family of interests (BSH genes), the blast database need to be non-redundant.

## Interpretation

### Pitfalls

BLAST is a local alignment tool. When the input query sequences is longer, multiple regions of the query sequences can be aligned to regions of the subject sequneces (the one in the database) separately, aka HSPs. The online BLAST interface sort the blast hits by *MAX_SCORE*, which is the best aligned HSP for each `qseqid|sseqid`. We want to reproduce the online blast results and visualization.

### Current approach

Therefore when the query sequences are longer sequence, the best local HSP alignment doesn't guarantee the best hits for the full contig sequences, and therefore we need to assemble all the HSPs together, calculate the sequence identity (`pident`), total alignment length (`total_aln_len`) and perhaps total bitscore (`total_bitscore`).


**constraint optimization problem**: for each query contig Q and for each reference sequence S, extend the highest-scoring fragment aligning Q to S with as many other non-overlapping fragments as possible, to form a collection of high scoring fragments HSP(Q, S) that do not overlap in Q.  Then, for each Q, output the S* with highest sum of bitscores over HSP(Q, S*).


### Alternative approach 1

For each `qseqid|sseqid`, we get the union of the `query_hsp_union_ranges` and `sub_hsp_union_ragnes`, extract the sequences `qseqs` and `sseqs`, calcuate the pairwise alignments (ANI), and generate the virsulization similar to [MAUVE](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC442156/). To choose the best `sseqid`, we can calcuate the ANI, query coverage, bitscore, and keep the top hits. (outfmt 5 is needed for this purpose).

### Alternative approache 2

Instead of blast the whole contig, split the contigs into non-overlapping 500 bps fragments and separately blast each fragment. To some extent, this can be indiciative or robust to recombination.

- `min_piden`: 60%
- `min_qcov`: 80%

## Q & A

### Why do we need the sequencd similarity filter?
A: We are not interested in low-similar alignment, and without filter those HSPs out, we would overestimate the total_alignmen_length.

### Why do we need the query coverage filter?
A: Database incomplete or contigs mis-assembly. Therefore, to define the **best hits**, we not only to assemble HSPs from the same `qseqid|sseqid`, but also need enough `qcov`.

### Structure of the algorithm

Give a list of HSPs (blast results):

1. Go through every HSP and perform a per-HSP filter based on sequence identity (`pident`)
2. Group HSPs by qseqid-sseqid pairs.
3. Within each group (HSPs are alread sorted by HSP-Bitscore by blast): keep a subset of disjoint HSPs that max the total `bit_score` while minimizing the overlaps.
- definition of 'disjoint'
4. Compute a score for each group.
5. Group sseqid by qseqid, and within each group, produce the winner.
