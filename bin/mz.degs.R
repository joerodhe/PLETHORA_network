###############################################################################
# Description:
#   This script takes two DEA result tables: one comparing MZ vs. DZ, and 
#   another comparing MZ vs. EZ. It identifies genes significantly upregulated 
#   in both comparisons (log2 fold-change > log2(1.5), adjusted p-value < 0.05), 
#   merges their log2 fold-change and p-value information, and writes out a 
#   combined table of these 'meristem upregulated' genes.
###############################################################################

###############################################################################
# Load required library (suppressing startup messages) and parse command-line
###############################################################################
suppressPackageStartupMessages(library(optparse))

option_list = list(
    make_option(c("-d", "--mzdz"), type = "character", default = NULL,
                help = "MZ vs DZ DEA results table.", metavar = "character"),
    make_option(c("-e", "--mzez"), type = "character", default = NULL,
                help = "MZ vs EZ DEA results table.", metavar = "character"),  
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Output results table.", metavar = "character")           
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

###############################################################################
# Step 1: Read DEA results tables for MZ-DZ and MZ-EZ comparisons
###############################################################################
dea.mzdz <- read.delim(opt$mzdz)
dea.mzez <- read.delim(opt$mzez)

###############################################################################
# Step 2: Filter for upregulated genes (log2FoldChange > log2(1.5), padj < 0.05)
###############################################################################
dea.mzdz.up <- subset(dea.mzdz, log2FoldChange > log2(1.5) & padj < 0.05)
dea.mzdz.up <- dea.mzdz.up[order(dea.mzdz.up$padj),]

dea.mzez.up <- subset(dea.mzez, log2FoldChange > log2(1.5) & padj < 0.05)
dea.mzez.up <- dea.mzez.up[order(dea.mzez.up$padj),]

###############################################################################
# Step 3: Identify genes upregulated in both MZ-DZ and MZ-EZ
###############################################################################
genes.mzdz.up <- rownames(dea.mzdz.up)
genes.mzez.up <- rownames(dea.mzez.up)
genes.mz <- intersect(genes.mzdz.up, genes.mzez.up)

###############################################################################
# Step 4: Gather log2FC and p-value information for the overlapping genes
###############################################################################
l2fc.mzez <- sapply(genes.mz, function(gen){
  dea.mzez.up$log2FoldChange[rownames(dea.mzez.up) == gen]
})
l2fc.mzdz <- sapply(genes.mz, function(gen){
  dea.mzdz.up$log2FoldChange[rownames(dea.mzdz.up) == gen]
})
padj.mzez <- sapply(genes.mz, function(gen){
  dea.mzez.up$padj[rownames(dea.mzez.up) == gen]
})
padj.mzdz <- sapply(genes.mz, function(gen){
  dea.mzdz.up$padj[rownames(dea.mzdz.up) == gen]
})

###############################################################################
# Step 5: Combine into a final data frame, including mean log2FC & mean padj
###############################################################################
meristem <- data.frame(
  'gene_id' = genes.mz,
  'log2FoldChange_MeristemElong' = l2fc.mzez,
  'log2FoldChange_MeristemDif'   = l2fc.mzdz,
  'padj_MeristemElong'           = padj.mzez,
  'padj_MeristemDif'             = padj.mzdz,
  'log2FoldChange_mean'          = rowMeans(data.frame(l2fc.mzdz, l2fc.mzez)),
  'padj_mean'                    = rowMeans(data.frame(padj.mzdz, padj.mzez))
)

meristem <- meristem[order(meristem$log2FoldChange_mean, decreasing = TRUE), ]

###############################################################################
# Step 6: Write final table of upregulated meristem genes to an output file
###############################################################################
write.table(
  meristem, 
  file = opt$output, 
  row.names = FALSE, 
  quote = FALSE, 
  sep = '\t'
)
