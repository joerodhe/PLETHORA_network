###############################################################################
# Purpose of this script:
#   This script reads a PLT-regulated gene compendium (Santuari et al. data),
#   intersects it with a set of Meristem DEGs, and outputs the subset of genes
#   that are both PLT-regulated and Meristem-upregulated. It then filters
#   a FASTA file to include only these intersecting genes and writes both
#   a list of selected gene IDs and the corresponding FASTA sequences.
###############################################################################

# Suppress library startup messages
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(seqinr))

# Define command-line options
option_list = list(
    make_option(c("-c", "--compendium"), type = "character", default = NULL,
                help = "Santurari's PLT-genes-regulated compendium", metavar = "character"),
    make_option(c("-d", "--degs"), type = "character", default = NULL,
                help = "MZ DEGs list", metavar = "character"),
    make_option(c("-f", "--fasta"), type = "character", default = NULL,
                help = "Fasta file to filter.", metavar = "character"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Output file.", metavar = "character"),     
    make_option(c("-v", "--fastaOutput"), type = "character", default = NULL,
                help = "Output fasta file.", metavar = "character")       
)

# Parse arguments
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# Read the Santuari compendium of PLT-regulated genes
compendium <- read_excel(
  opt$compendium,
  sheet = "Compendium", 
  range = "A8:AB20534", 
  na = "NA", 
  col_names = c(
    "AGI", 
    "ProbeID", 
    "Gene Symbol", 
    "Description", 
    "PLT1_FC", 
    "PLT2_FC", 
    "PLT3_FC", 
    "PLT4_FC", 
    "PLT5_FC", 
    "PLT7_FC", 
    "QC_1h_FC", 
    "QC_4h_FC", 
    "PLT1_p", 
    "PLT2_p", 
    "PLT3_p", 
    "PLT4_p", 
    "PLT5_p", 
    "PLT7_p", 
    "QC_1h_p", 
    "QC_4h_p", 
    "PLT1", 
    "PLT2", 
    "PLT3", 
    "PLT4", 
    "PLT5", 
    "PLT7", 
    "QC_1h_plt", 
    "QC_4h_plt"
  )
)

# Identify which genes are regulated (activated or repressed) by any PLT
compendium.regulated.gen <- apply(compendium[,21:26], MARGIN = 1, function(gen)
  1 %in% gen | -1 %in% gen
)
compendium.plt.targets <- compendium$AGI[compendium.regulated.gen]

# Read a list of Meristem DEGs
mz.degs <- readLines(opt$degs)

# Intersect the compendium PLT targets with the Meristem DEGs
mz.compendium.plt.targets <- intersect(compendium.plt.targets, mz.degs)

# Write out the intersecting genes to a text file
writeLines(mz.compendium.plt.targets, opt$output)

# Read the full FASTA file and filter sequences belonging to the intersecting genes
fasta <- read.fasta(opt$fasta, seqtype = "AA")
fasta.annot <- sapply(fasta, attr, "Annot")
fasta.locus <- sapply(fasta.annot, function(annot) strsplit(annot, "locus=")[[1]][2])

# Subset only those sequences corresponding to the intersecting gene set
fasta.mz.compendium <- fasta[fasta.locus %in% mz.compendium.plt.targets]

# Prepare names for writing to new FASTA
ids <- sapply(fasta.mz.compendium, attr, "Annot")
ids <- gsub(ids, pattern = ">", replacement = "")

# Write filtered sequences to the specified FASTA output
write.fasta(
  fasta.mz.compendium,
  ids,
  opt$fastaOutput
)
