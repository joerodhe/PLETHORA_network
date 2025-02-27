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
suppressPackageStartupMessages(library(ggrepel))


vplot <- function(res, name, plts, xl=8, yl=20){
  
  res <- na.omit(as.data.frame(res))
  
  res = mutate(res, sig=ifelse(res$padj < 0.05, "Sig", "Not DE"))
  res$sig[res$sig == "Sig" & res$log2FoldChange > log2(1)] <- "Up"
  res$sig[res$sig == "Sig" & res$log2FoldChange < -log2(1)] <- "Down"
  res$id <- rownames(res)
  rownames(plts) <- plts$id
  res$plt <- NA
  res$plt[res$id %in% plts$id] <- plts[res$id[res$id %in% plts$id], "name"]
  
  plot <- ggplot(res, aes(x=log2FoldChange, y=-log10(padj))) + 
    geom_point(aes(col=sig), alpha = 0.4) +
    geom_text_repel(#data=res %>% filter(log2FoldChange > 3 & padj < 0.01),
      data=res %>% filter(id %in% plts$id),
      aes(label=plt), 
      na.rm = T, 
      size = 2.8) +
    theme_minimal() +
    scale_color_brewer(palette="Dark2") +
    geom_hline(yintercept=-log10(0.05),linetype = "longdash", col="red") +
    annotate(geom = "text", x=-xl+2, y=2, label="-log10(p = 0.05)", color="red") +
    geom_vline(xintercept=c(-log2(1), log2(1)), linetype = "dashed") +
    annotate(geom = "text", x=log2(1)+0.5, y=yl-2, label="log2(1)", angle=270) +
    scale_y_continuous(breaks = seq(0,yl,5), limits = c(0,yl)) +
    scale_x_continuous(breaks = seq(-xl,xl, 2), limits = c(-xl,xl)) +
    labs(title = paste("Differential Expression:", name), 
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

dea <- function(dds, name, contrast, foldChange, output.path, organism, plts){
  cat("Running DEA", name, "\n")
  res <- results(
    dds, 
    name = name, 
    contrast = contrast,
    alpha=0.05, 
    lfcThreshold = foldChange)
    )
  summary(res)
  write.table(
    res, 
    file.path(output.path, paste("dea", paste(tolower(contrast[-1]), collapse = "-"), "tsv", sep = ".")),
    sep="\t", 
    quote = FALSE
    )
  res.vplot <- vplot(res, name, plts, xl=10, y= 25)
  res <- res[(!is.na(res$log2FoldChange)) & (!is.na(res$padj)),]
  res.gen.up <- rownames(res[res$padj < 0.05 & res$log2FoldChange > log2(foldChange),])
  res.gen.down <- rownames(res[res$padj < 0.05 & res$log2FoldChange < log2(foldChange),])
  res.lfc <- setNames(object = res$log2FoldChange, rownames(res))
  res.lfc <- sort(res.lfc, decreasing = TRUE)
  gsea.result <- GSEA(
    res.lfc, 
    TERM2GENE = TERM2GENE, 
    TERM2NAME = TERM2NAME, 
    eps = 0, 
    nPermSimple = 10000, 
    seed = TRUE
    )
  res.enrichment <- enrich_plot(gsea.result, name)
  vplots <- grid.arrange(
    res.vplot, 
    res.enrichment, 
    nrow = 1, 
    top =paste("Differential expression analysis for", organism)
    )
  ggsave(
    plot = vplots, 
    width = 14, 
    height = 6, 
    units = "in",
    filename = file.path(output.path, paste("dea", paste(tolower(contrast[-1]), collapse = "-"), "svg", sep = "."))
    )
    return(
      list(
        res.gen.up,
        res.gen.down
      )
    )
}

plot_venn <- function(output.path, g1, g2, zone, n1, n2, organism){
  svg(file.path(output.path, paste("venn", zone, "DEGs.svg", sep=".")))
  grid.newpage()
  venn.plot <- draw.pairwise.venn(
    area1 = length(g1), 
    area2 = length(g2), 
    cross.area = length(intersect(g1, g2)),
    category = c(
      paste("Upregulated DEGs in", zone, "\n", n1), 
      paste("Upregulated DEGs in", zone, "\n", n2)
      ), 
    fill = c("skyblue", "pink1"), 
    lty = "blank", cat.cex = 0.9, scaled = F, cat.pos = 180
    )
  grid.arrange(gTree(children = venn.plot), 
      top=paste("DEA:", zone, "up-regulated DEGs in", organism),
      bottom = paste('Union:', length(union(g1, g2))))
  dev.off()

  degs.up <- intersect(g1, g2)
  writeLines(degs.up, file.path(output.path, paste("upregulated", zone, "DEGs.txt", sep=".")))
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
              help = "Organism name.", metavar = "character")  ,
    make_option(c("-p", "--plts"), type = "character", default = NULL,
              help = "plt table", metavar = "character")   
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

#opt$counts <- "output/Arabidopsis_thaliana/FEATURECOUNTS/TAIR10.gene.counts.tsv"
#opt$metadata <- "temp/Arabidopsis_thaliana/metadata.tsv"
#opt$foldChange <- 1
#opt$term2gene <- "input/Arabidopsis_thaliana/GENOME/TERM2GENE.tsv"
#opt$term2name <- "workflow/TERM2NAME.tsv"
#opt$output <- "output/Arabidopsis_thaliana/temp"
#opt$organism <- "Arabidopsis_thaliana"
#opt$plts <- "others/plt.ids.txt"

count.matrix <- read.delim(opt$counts, check.names = FALSE)
metadata <- read.delim(opt$metadata)

TERM2GENE <- read.delim(opt$term2gene)
TERM2NAME <- read.delim(opt$term2name)

TERM2GENE <- TERM2GENE[TERM2GENE$go.id %in% TERM2NAME$go.id, ]

dds <- DESeqDataSetFromMatrix(countData=count.matrix, colData = metadata, design=~dex)

dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
#plot(sizeFactors(dds), colSums(counts(dds)))
#abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))
#normlzd_dds <- as.data.frame(counts(dds, normalized=T))

plts <- read.delim(opt$plts, head = FALSE, col.names = c("name", "id"))

contrasts <- list(
  c("dex", "MZ", "DZ"),
  c("dex", "MZ", "EZ"),
  c("dex", "DZ", "EZ")
)

dea.results <- list()

for(i in 1:length(contrasts)){
 dea.results[[i]] <- dea(
  dds, 
  paste(contrasts[[i]][-1], collapse = "_vs_"), 
  contrasts[[i]], 
  opt$foldChange, 
  opt$output, 
  opt$organism, 
  plts
  )
}

zones <- c("MZ", "DZ", "EZ")

for(zone in zones){
  c <- sapply(contrasts, function(x) zone %in% x)
  idx <- which(c)
  zone.genes <- lapply(idx, function(i) dea.results[[i]][[which(zone == contrasts[[i]][-1])]])
  c.names <- sapply(contrasts[c], function(x) paste(x[-1], collapse="-"))
  plot_venn(opt$output, zone.genes[[1]], zone.genes[[2]], zone, c.names[1], c.names[2], opt$organism)
}
