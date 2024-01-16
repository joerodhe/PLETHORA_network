suppressPackageStartupMessages(library(optparse))


option_list = list(
    make_option(c("-d", "--mzdz"), type = "character", default = NULL,
              help = "MZ vz DZ DEA results table.", metavar = "character"),
    make_option(c("-e", "--mzez"), type = "character", default = NULL,
              help = "MZ vz EZ DEA results table.", metavar = "character"),  
    make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output results table.", metavar = "character")           
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

dea.mzdz <- read.delim(opt$mzdz)
dea.mzez <- read.delim(opt$mzez)

dea.mzdz.up <- subset(dea.mzdz, log2FoldChange > log2(1.5))
dea.mzdz.up <- subset(dea.mzdz.up, padj<.05)
dea.mzdz.up <- dea.mzdz.up[order(dea.mzdz.up$padj),]
dea.mzez.up <- subset(dea.mzez, log2FoldChange > log2(1.5))
dea.mzez.up <- subset(dea.mzez.up, padj<.05)
dea.mzez.up <- dea.mzez.up[order(dea.mzez.up$padj),]

genes.mzdz.up <- rownames(dea.mzdz.up)
genes.mzez.up <- rownames(dea.mzez.up)

genes.mz <- intersect(genes.mzdz.up, genes.mzez.up)

l2fc.mzez <- sapply(genes.mz, function(gen){
  dea.mzez.up$log2FoldChange[rownames(dea.mzez.up) ==  gen]})
l2fc.mzdz <- sapply(genes.mz, function(gen){
  dea.mzdz.up$log2FoldChange[rownames(dea.mzdz.up) ==  gen]})
padj.mzez <- sapply(genes.mz, function(gen){
  dea.mzez.up$padj[rownames(dea.mzez.up) ==  gen]})
padj.mzdz <- sapply(genes.mz, function(gen){
  dea.mzdz.up$padj[rownames(dea.mzdz.up) ==  gen]})

meristem <- data.frame(
  'gene_id' = genes.mz,
  'log2FoldChange_MeristemElong' = l2fc.mzez,
  'log2FoldChange_MeristemDif' = l2fc.mzdz,
  'padj_MeristemElong' = padj.mzez,
  'padj_MeristemDif' = padj.mzdz,
  'log2FoldChange_mean' = rowMeans(data.frame(l2fc.mzdz, l2fc.mzez)),
  'padj_mean' = rowMeans(data.frame(padj.mzdz, padj.mzez))
  )

meristem <- meristem[order(meristem$log2FoldChange_mean, decreasing = T), ]
write.table(meristem, file = opt$output, row.names = F, quote = FALSE, sep = '\t')
