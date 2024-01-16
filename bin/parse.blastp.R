suppressPackageStartupMessages(library(optparse))



option_list = list(
    make_option(c("-q", "--query"), type = "character", default = NULL,
              help = "Query protein fasta file.", metavar = "path"),
    make_option(c("-s", "--subject"), type = "character", default = NULL,
              help = "Subject protein fasta file.", metavar = "path"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Best reciprocal BLASTp hit table file.", metavar = "path"),
    make_option(c("-c", "--cores"), type = "integer", default = 1,
              help = "Cores to use in BLASTp parallel processing.", metavar = "integer")     
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

opt$orgs <- 'Arabidopsis thaliana,Solanum lycopersicum'

orgs <- strsplit(opt$orgs, ',')[[1]]

sol <- read.delim('output/Solanum_lycopersicum/BLASTp/orthologs.tsv')
sol <- sol[,c('query_id', 'subject_id')]
sol <- sol[!duplicated(sol$subject_id),]
cuc <- read.delim('output/Cucumis_sativus/BLASTp/orthologs.tsv')
cuc <- cuc[,c('query_id', 'subject_id')]
gly <- read.delim('output/Glycine_max/BLASTp/orthologs.tsv')
gly <- gly[,c('query_id', 'subject_id')]

x <- merge(sol, cuc, by = 'subject_id')
y <- merge(x, gly, by = 'subject_id')

library(tidyverse)

#ERROR: dependencies ‘googledrive’, ‘googlesheets4’ are not available for package ‘tidyverse’