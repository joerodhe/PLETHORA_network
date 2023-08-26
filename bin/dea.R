suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(gridExtra))


vplot <- function(res, f1, f2, xl=8, yl=20){
  
  res <- na.omit(as.data.frame(res))
  
  res = mutate(res, sig=ifelse(res$padj < 0.05, "Sig", "Not DE"))
  res$sig[res$sig == "Sig" & res$log2FoldChange > log2(1.5)] <- "Up"
  res$sig[res$sig == "Sig" & res$log2FoldChange < -log2(1.5)] <- "Down"
  res$id <- rownames(res)
  
  plot <- ggplot(res, aes(x=log2FoldChange, y=-log10(padj))) + 
  geom_point(aes(col=sig), alpha = 0.4) +
  geom_text(data=res %>% filter(log2FoldChange > 3 & padj < 0.01),
             aes(label=id), check_overlap = T, na.rm = T, size = 2.8) +
  theme_minimal() +
  scale_color_brewer(palette="Accent") +
  geom_hline(yintercept=-log10(0.05),linetype = "longdash", col="red") +
  annotate(geom = "text", x=-xl+2, y=2, label="-log10(p = 0.05)", color="red") +
  geom_vline(xintercept=c(-log2(1.5), log2(1.5)), linetype = "dashed") +
  annotate(geom = "text", x=log2(1.5)+0.5, y=yl-2, label="log2(1.5)", angle=270) +
  scale_y_continuous(breaks = seq(0,yl,5), limits = c(0,yl)) +
  scale_x_continuous(breaks = seq(-xl,xl, 2), limits = c(-xl,xl)) +
  labs(title = paste("Differential Expression:", f1, "compared with", f2), 
       x = "Log2(Fold-Change)", 
       y = "-Log10(P-adjusted value)", color = "Expression")
  
  return(plot)
}

enrich_plot <- function(gsea.result, comparison){
  gsea.result <- as.data.frame(gsea.result)
  plot <- ggplot()
  if(nrow(gsea.result) > 0){
  if(nrow(gsea.result) > 20) max.n <- 20
  else max.n <- nrow(gsea.result)
  gsea.result <- head(gsea.result[order(gsea.result$p.adjust),],max.n)

  plot <- ggplot(gsea.result, aes(NES, fct_reorder(Description, NES))) +
    geom_segment(aes(xend=0, yend = Description)) +
    geom_point(aes(color=p.adjust, size = setSize)) +
    scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
    trans = "log10",
    guide=guide_colorbar(reverse=TRUE, order=1)) +
    scale_size_continuous(range=c(2, 10)) +
    theme_dose(12) +
    xlab("NES") +
    ylab(NULL) +
    ggtitle(paste0("BP GO terms enriched in DEGs (", comparison, ")")) +
        scale_y_discrete(labels = function(x) str_wrap(x, 45)) +
    theme(axis.text.y = element_text(size = 10))
  }
  return(plot)
}

option_list = list(
    make_option(c("-m", "--metadata"), type = "character", default = NULL,
              help = "Path to metadata file.", metavar = "character"),
    make_option(c("-c", "--counts"), type = "character", default = NULL,
              help = "Path to parsed featurecounts matrix.", metavar = "character"),
    make_option(c("-f", "--foldChange"), type = "double", default = NULL,
              help = "Fol-Change to define DEGs.", metavar = "character"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output path.", metavar = "character"),
    make_option(c("-g", "--term2gene"), type = "character", default = NULL,
              help = "Path to term2gene table file.", metavar = "character"),
    make_option(c("-d", "--term2name"), type = "character", default = NULL,
              help = "Path to term2name table file.", metavar = "character"),
    make_option(c("-u", "--organism"), type = "character", default = NULL,
              help = "Organism name.", metavar = "character")     
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

#opt$counts <- "output/Arabidopsis_thaliana/FEATURECOUNTS/TAIR10.gene.counts.tsv"
#opt$metadata <- "temp/Arabidopsis_thaliana/metadata.tsv"
#opt$foldChange <- 2
#opt$term2gene <- "input/Arabidopsis_thaliana/GENOME/TERM2GENE.tsv"
#opt$term2name <- "workflow/TERM2NAME.tsv"

count.matrix <- read.delim(opt$counts, check.names = FALSE)
metadata <- read.delim(opt$metadata)

TERM2GENE <- read.delim(opt$term2gene)
TERM2NAME <- read.delim(opt$term2name)

TERM2GENE <- TERM2GENE[TERM2GENE$go.id %in% TERM2NAME$go.id, ]

dds <- DESeqDataSetFromMatrix(countData=count.matrix, colData = metadata, design=~dex)

dds <- estimateSizeFactors(dds)
#plot(sizeFactors(dds), colSums(counts(dds)))
#abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))
#normlzd_dds <- as.data.frame(counts(dds, normalized=T))

cat('Running MZ-DZ comparison....\n')
contrast <- c("dex", "MZ", "DZ")
dds <- DESeq(dds)
res.md <- results(dds, name = "MZ vs DZ", 
               contrast = contrast,
               alpha=0.05, lfcThreshold = log2(as.numeric(opt$foldChange)))
summary(res.md)

write.table(res.md, file.path(opt$output, "dea.mz-dz.tsv"), sep="\t", quote = FALSE)

md.vplot <- vplot(res.md, "MZ", "DZ", xl=10, y= 25)
res.md <- res.md[(!is.na(res.md$log2FoldChange)) & (!is.na(res.md$padj)),]
res.md.gen.up <- rownames(res.md[res.md$padj < 0.05 & res.md$log2FoldChange > log2(opt$foldChange),])
res.md.lfc <- setNames(object = res.md$log2FoldChange, rownames(res.md))
res.md.lfc <- sort(res.md.lfc, decreasing = TRUE)

gsea.result = GSEA(res.md.lfc, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, eps = 0, nPermSimple = 10000, seed = TRUE)
md.enrichment <- enrich_plot(gsea.result, "MZ-DZ")


vplots <- grid.arrange(md.vplot, md.enrichment, 
    nrow = 1, top =paste("Differential expression analysis for", opt$organism))

ggsave(plot = vplots, width = 14, height = 6, units = "in",
    filename = file.path(opt$output, "dea.mz-dz.svg"))

cat('Running MZ-EZ comparison....\n')
contrast <- c("dex", "MZ", "EZ")
res.me <- results(dds, name = "MZ vs EZ", 
               contrast = contrast,
               alpha=0.05, lfcThreshold = log2(as.numeric(opt$foldChange)))
summary(res.me)

write.table(res.me, file.path(opt$output, "dea.mz-ez.tsv"), sep="\t", quote = FALSE)

me.vplot <- vplot(res.me, "MZ", "EZ", xl=10, y= 25)
res.me <- res.me[(!is.na(res.me$log2FoldChange)) & (!is.na(res.me$padj)),]
res.me.gen.up <- rownames(res.me[res.me$padj < 0.05 & res.me$log2FoldChange > log2(opt$foldChange),])
res.me.lfc <- setNames(object = res.me$log2FoldChange, rownames(res.me))
res.me.lfc <- sort(res.me.lfc, decreasing = TRUE)

gsea.result = GSEA(res.me.lfc, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, eps = 0, nPermSimple = 10000, seed = TRUE)
me.enrichment <- enrich_plot(gsea.result, "MZ-EZ")


vplots <- grid.arrange(me.vplot, me.enrichment, 
    nrow = 1, top = paste("Differential expression analysis for", opt$organism))

ggsave(plot = vplots, width = 14, height = 6, units = "in",
    filename = file.path(opt$output, "dea.mz-ez.svg"))

svg(file.path(opt$output, "venn.DEGs.svg"))
grid.newpage()
venn.plot <- draw.pairwise.venn(
  area1 = length(res.md.gen.up), 
  area2 = length(res.me.gen.up), 
  cross.area = length(intersect(res.md.gen.up, res.me.gen.up)),
  category = c("Upregulated DEGs in meristem\n(MZ-DZ)", 
               "Upregulated DEGs in meristem\n(MZ-EZ)"), 
  fill = c("skyblue", "pink1"), 
  lty = "blank", cat.cex = 0.9, scaled = F, cat.pos = 180
  )
grid.arrange(gTree(children = venn.plot), 
    top=paste('DEA: Meristem up-regulated DEGs in', opt$organism),
    bottom = paste('Union:', length(union(res.md.gen.up, res.me.gen.up))))
dev.off()

degs.up <- intersect(res.md.gen.up, res.me.gen.up)
writeLines(degs.up, file.path(opt$output, "upregulated.DEGs.txt"))
