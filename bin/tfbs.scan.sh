#!/bin/bash

###############################################################################
# Description:
#   This script uses RSAT (Regulatory Sequence Analysis Tools) to perform:
#   1) `retrieve-seq`: extracting upstream sequences for a list of gene IDs
#   2) `matrix-scan`: scanning those extracted sequences with a given motif 
#      (in TRANSFAC format).
#   The user must supply the organism name, a file containing gene IDs, and 
#   a motif file. The script outputs the scan results in an RSAT ".ft" format 
#   file, suitable for further processing or analysis.
###############################################################################

# Create a temporary file to store the extracted FASTA sequences
FASTA=$(mktemp)

# Echo the name of the final output file (transforming 'txt' to 'ft' in the input)
echo "${2/txt/'ft'}"

###############################################################################
# Step 1: Retrieve upstream sequences
###############################################################################
# -org      : The organism name recognized by RSAT (e.g., "Arabidopsis_thaliana")
# -feattype : The feature type (gene)
# -type     : "upstream" to retrieve promoter/upstream regions
# -format   : Format of output sequences (fasta)
# -label id : Label sequences using IDs only
# -noorf    : Exclude coding region from the retrieval
# -ids_only : Do not retrieve sequences for IDs not found
# -i        : Input file containing gene IDs
# -o        : Output file for the retrieved sequences
###############################################################################
rsat retrieve-seq \
 -org "$1" \
 -feattype gene \
 -type upstream \
 -format fasta \
 -label id \
 -noorf \
 -ids_only \
 -i "$2" \
 -o "$FASTA"

###############################################################################
# Step 2: Run matrix-scan
###############################################################################
# -v 1            : Verbosity level
# -matrix_format  : Format of the motif file (TRANSFAC)
# -m <motif_file> : Path to motif file
# -pseudo 1       : Pseudocount for matrix scoring
# -decimals 1     : Number of decimal places in output
# -2str           : Scan both strands
# -origin end     : Defines how to handle the coordinate system
# -bgfile         : Path to background frequencies file 
# -bg_pseudo 0.01 : Pseudocount for background
# -return sites   : Return only significant sites
# -lth score 1    : Lower threshold (score >= 1)
# -uth pval 1e-4  : Upper threshold on p-value (<= 1e-4)
# -i <fasta>      : Input sequences in FASTA format
# -seq_format     : Format of input sequences (fasta)
# -n score        : Sorting or returning top results by score
# -o <file.ft>    : Output file with ".ft" extension
###############################################################################
rsat matrix-scan \
 -v 1 \
 -matrix_format transfac \
 -m "$3" \
 -pseudo 1 \
 -decimals 1 \
 -2str \
 -origin end \
 -bgfile "$CONDA_PREFIX/share/rsat/public_html/data/genomes/$1/oligo-frequencies/2nt_upstream-noorf_$1-ovlp-1str.freq.gz" \
 -bg_pseudo 0.01 \
 -return sites \
 -lth score 1 \
 -uth pval 1e-4 \
 -i "$FASTA" \
 -seq_format fasta \
 -n score \
 -o "${2/txt/'ft'}"
