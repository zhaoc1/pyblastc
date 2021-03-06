#!/usr/bin/env python3

# set to -1 to unlimit
MAX_LINES =  -1 #1 * 1000 * 1000 #-1

MIN_OVERLAP = 0.1
MIN_PIDENT = 90
MIN_QCOV = 80
MAX_TARGET_SEQS = 2500

# FROM IDseq for NT only
MIN_ALIGNMENT_LENGTH = 36

blast_outfmt6_schema = {
    "qseqid": str,
    "sseqid": str,
    "pident": float,
    "qlen": int,
    "slen": int,
    "hsplen": int,
    "mismatch": int,
    "gapopen": int,
    "qstart": int,
    "qend": int,
    "sstart": int,
    "send": int,
    "evalue": float,
    "bitscore": float,
}

ranked_blast_output_schema = dict(blast_outfmt6_schema)
ranked_blast_output_schema.update({
    "qcov": float,
    "hsp_count": int
})
