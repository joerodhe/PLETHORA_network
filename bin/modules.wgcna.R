suppressPackageStartupMessages(library(optparse))

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

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


opt$subject <- 'temp/Arabidopsis_thaliana/WGCNA/WGCNA.RData'
opt$query <- 'temp/Cucumis_sativus/WGCNA/WGCNA.RData'
opt$orts <- 'output/Cucumis_sativus/BLASTp/orthologs.tsv'
opt$plts <- 'others/plt_ids.txt'

plts.table <- read.delim(opt$plts, col.names = c('Name', 'ID'), header = FALSE)

orts.modules <- read.delim(opt$orts)
orts.modules <- orts.modules[,1:2]

load(opt$query)

names(dynamicColors) <- gene.names
orts.modules <- orts.modules[which(orts.modules$query_id %in% intersect(orts.modules$query_id, gene.names)),]
orts.modules$query_module <- dynamicColors[orts.modules$query_id]

load(opt$subject)

names(dynamicColors) <- gene.names
orts.modules <- orts.modules[which(orts.modules$subject_id %in% intersect(orts.modules$subject_id, gene.names)),]
orts.modules$subject_module <- dynamicColors[orts.modules$subject_id]

orts.modules[which(orts.modules$subject_id %in% plts.table$ID), ]
