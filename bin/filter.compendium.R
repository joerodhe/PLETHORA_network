suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(seqinr))

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

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

compendium <- read_excel(opt$compendium,
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
    "QC_4h_plt"))

compendium.regulated.gen <- apply(compendium[,21:26], MARGIN = 1, function(gen) 1%in%gen | -1%in%gen)
compendium.plt.targets <- compendium$AGI[compendium.regulated.gen]

mz.degs <- readLines(opt$degs)

mz.compendium.plt.targets <- intersect(compendium.plt.targets,mz.degs)

writeLines(mz.compendium.plt.targets, opt$output)

fasta <- read.fasta("input/Arabidopsis_thaliana/GENOME/TAIR10.protein.fa", seqtype = "AA")

fasta.annot <- sapply(fasta, attr, "Annot")
fasta.locus <- sapply(fasta.annot, function(annot) strsplit(annot, "locus=")[[1]][2])

fasta.mz.compendium <- fasta[fasta.locus %in% mz.compendium.plt.targets]
ids <- sapply(fasta.mz.compendium, attr, "Annot")
ids <- gsub(ids, pattern = ">", replacement = "")

write.fasta(
    fasta.mz.compendium,
    ids,
    opt$fastaOutput
)