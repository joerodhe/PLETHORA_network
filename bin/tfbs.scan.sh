#!/bin/bash

FASTA=$(mktemp)

echo \"${2/txt/'ft'}\"

rsat retrieve-seq\
 -org $1\
 -feattype gene\
 -type upstream\
 -format fasta\
 -label id\
 -noorf\
 -ids_only\
 -i $2\
 -o $FASTA

rsat matrix-scan\
 -v 1\
 -matrix_format transfac\
 -m $3\
 -pseudo 1\
 -decimals 1\
 -2str\
 -origin end\
 -bgfile $CONDA_PREFIX/share/rsat/public_html/data/genomes/$1/oligo-frequencies/2nt_upstream-noorf_$1-ovlp-1str.freq.gz\
 -bg_pseudo 0.01\
 -return sites\
 -lth score 1\
 -uth pval 1e-4\
 -i $FASTA\
 -seq_format fasta\
 -n score\
 -o \"${2/txt/'ft'}\"
