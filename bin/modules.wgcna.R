###############################################################################
# This script reads previously computed WGCNA RData files (query and subject),
# a table of orthologous gene IDs, and a table of PLT gene IDs. It matches 
# query-subject gene pairs to their WGCNA module assignments, then identifies 
# those subject genes that are PLTs.
###############################################################################

# Suppress startup messages from the optparse library
suppressPackageStartupMessages(library(optparse))

# Define command-line options
option_list <- list(
    make_option(c("-q", "--query"), type = "character", default = NULL,
              help = "Query protein fasta file.", metavar = "path"),
    make_option(c("-s", "--subject"), type = "character", default = NULL,
              help = "Subject protein fasta file.", metavar = "path"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Path to save output files.", metavar = "path"),
    make_option(c("-c", "--cores"), type = "integer", default = 1,
              help = "Cores to use in BLASTp parallel processing.", metavar = "integer")     
)

# Parse command-line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# For demonstration/testing, manually set file paths
opt$subject <- 'temp/Arabidopsis_thaliana/WGCNA/WGCNA.RData'
opt$query <- 'temp/Cucumis_sativus/WGCNA/WGCNA.RData'
opt$orts <- 'output/Cucumis_sativus/BLASTp/orthologs.tsv'
opt$plts <- 'others/plt_ids.txt'

# Read a table of PLT gene IDs (two columns: Name, ID)
plts.table <- read.delim(opt$plts, col.names = c('Name', 'ID'), header = FALSE)

# Read orthologs table (assuming first 2 columns: query_id, subject_id)
orts.modules <- read.delim(opt$orts)
orts.modules <- orts.modules[,1:2]

# Load the WGCNA results from the "query" species
load(opt$query)

# Assign module colors to gene names for the query
names(dynamicColors) <- gene.names

# Filter the orthologs to only those query IDs present in this WGCNA dataset
orts.modules <- orts.modules[
  which(orts.modules$query_id %in% intersect(orts.modules$query_id, gene.names)),
]

# Map the query modules onto the orthologs data frame
orts.modules$query_module <- dynamicColors[orts.modules$query_id]

# Load the WGCNA results from the "subject" species
load(opt$subject)

# Assign module colors to gene names for the subject
names(dynamicColors) <- gene.names

# Filter the orthologs to only those subject IDs present in the subject dataset
orts.modules <- orts.modules[
  which(orts.modules$subject_id %in% intersect(orts.modules$subject_id, gene.names)),
]

# Map the subject modules onto the orthologs data frame
orts.modules$subject_module <- dynamicColors[orts.modules$subject_id]

# Finally, see which subject genes are in the PLT table
orts.modules[ which(orts.modules$subject_id %in% plts.table$ID), ]
