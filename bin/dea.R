###############################################################################
# Description:
#   This script reads an expression count matrix, metadata, and GO term 
#   annotations, then performs differential expression analyses (DEA) 
#   using DESeq2. It subsequently generates volcano plots, GSEA analyses, 
#   and Venn diagrams to visualize overlap among upregulated DEGs in multiple 
#   root zones (MZ, DZ, EZ). Additionally, it highlights PLT (PLT transcription 
#   factors) in the volcano plots for emphasis.
###############################################################################

###############################################################################
# Load packages needed for DEA, plots, and enrichment analyses, 
# suppressing startup messages for a cleaner output
###############################################################################
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

###############################################################################
# Function: vplot
# ---------------
# Description:
#   Creates a volcano plot with log2(FoldChange) on the X-axis and 
#   -log10(adjusted p-value) on the Y-axis. Color codes significant up/down 
#   genes and labels PLT genes with text repelling for clarity.
#
# Parameters:
#   res  : DESeq2 results object (can be coerced to data frame).
#   name : A label for the plot (e.g., "MZ vs DZ").
#   plts : Data frame of PLT genes with columns 'name' and 'id'.
#   xl, yl : Numeric bounds for X and Y axes.
#
# Returns:
#   A ggplot object (volcano plot).
###############################################################################
vplot <- function(res, name, plts, xl=8, yl=20){
  
  # Convert results to data frame and remove NA
  res <- na.omit(as.data.frame(res))
  
  # Identify significant genes based on padj < 0.05
  res = mutate(res, sig=ifelse(res$padj < 0.05, "Sig", "Not DE"))
  # Mark up/down with log2FoldChange > +1 or < -1
  res$sig[res$sig == "Sig" & res$log2FoldChange > log2(1)] <- "Up"
  res$sig[res$sig == "Sig" & res$log2FoldChange < -log2(1)] <- "Down"
  
  # Prepare labeling
  res$id <- rownames(res)
  rownames(plts) <- plts$id
  res$plt <- NA
  res$plt[res$id %in% plts$id] <- plts[res$id[res$id %in% plts$id], "name"]
  
  # Construct volcano plot
  plot <- ggplot(res, aes(x=log2FoldChange, y=-log10(padj))) + 
    geom_point(aes(col=sig), alpha = 0.4) +
    geom_text_repel(
      data=res %>% filter(id %in% plts$id),
      aes(label=plt), 
      na.rm = TRUE, 
      size = 2.8
    ) +
    theme_minimal() +
    scale_color_brewer(palette="Dark2") +
    # Horizontal line at p=0.05 threshold
    geom_hline(yintercept=-log10(0.05),linetype = "longdash", col="red") +
    annotate(geom = "text", x=-xl+2, y=2, label="-log10(p = 0.05)", color="red") +
    # Vertical lines at +-1 log2 fold-change
    geom_vline(xintercept=c(-log2(1), log2(1)), linetype = "dashed") +
    annotate(geom = "text", x=log2(1)+0.5, y=yl-2, label="log2(1)", angle=270) +
    scale_y_continuous(breaks = seq(0,yl,5), limits = c(0,yl)) +
    scale_x_continuous(breaks = seq(-xl,xl, 2), limits = c(-xl,xl)) +
    labs(
      title = paste("Differential Expression:", name), 
      x = "Log2(Fold-Change)", 
      y = "-Log10(P-adjusted value)", 
      color = "Expression"
    )
  
  return(plot)
}

###############################################################################
# Function: enrich_plot
# ---------------------
# Description:
#   Creates a dot-plot showing the normalized enrichment score (NES) for the
#   top GO terms from a GSEA result, colored by adjusted p-value and sized by 
#   the number of genes in each GO set.
#
# Parameters:
#   gsea.result : A GSEA object from clusterProfiler
#   comparison  : A string label used in the plot title
#
# Returns:
#   A ggplot object with top enriched GO terms.
###############################################################################
enrich_plot <- function(gsea.result, comparison){
  gsea.result <- as.data.frame(gsea.result)
  plot <- ggplot()
  
  if(nrow(gsea.result) > 0){
    # Plot up to 20 terms, ordered by p.adjust
    if(nrow(gsea.result) > 20) max.n <- 20
    else max.n <- nrow(gsea.result)
    
    gsea.result <- head(gsea.result[order(gsea.result$p.adjust),],max.n)
    
    plot <- ggplot(gsea.result, aes(NES, fct_reorder(Description, NES))) +
      geom_segment(aes(xend=0, yend = Description)) +
      geom_point(aes(color=p.adjust, size = setSize)) +
      scale_color_gradientn(
        colours=c("#f7ca64", "#46bac2", "#7e62a3"),
        trans = "log10",
        guide=guide_colorbar(reverse=TRUE, order=1)
      ) +
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

###############################################################################
# Function: dea
# -------------
# Description:
#   Performs DESeq2 differential expression analysis (DEA) for a given contrast,
#   extracts up/down genes above a foldChange threshold, runs GSEA, 
#   and creates summary plots (volcano + enrichment).
#
# Parameters:
#   dds         : A DESeqDataSet object, already run through DESeq().
#   name        : String label for naming results/plots.
#   contrast    : A DESeq2 contrast vector, e.g., c("dex", "MZ", "DZ").
#   foldChange  : Minimum fold-change threshold (passed as log2).
#   output.path : Where to save output files.
#   organism    : Organism name for labeling plots.
#   plts        : PLT gene table (for highlighting on the volcano plot).
#
# Returns:
#   A list of length 2:
#     1) Vector of significantly upregulated gene IDs
#     2) Vector of significantly downregulated gene IDs
###############################################################################
dea <- function(dds, name, contrast, foldChange, output.path, organism, plts){
  cat("Running DEA", name, "\n")
  res <- results(
    dds, 
    name = name, 
    contrast = contrast,
    alpha = 0.05, 
    lfcThreshold = foldChange
  )
  summary(res)
  
  # Save DEA results table
  write.table(
    res, 
    file.path(output.path, paste("dea", paste(tolower(contrast[-1]), collapse = "-"), "tsv", sep = ".")),
    sep="\t", 
    quote = FALSE
  )
  
  # Generate volcano plot
  res.vplot <- vplot(res, name, plts, xl=10, y= 25)
  
  # Identify up/down genes
  res <- res[(!is.na(res$log2FoldChange)) & (!is.na(res$padj)),]
  res.gen.up <- rownames(res[res$padj < 0.05 & res$log2FoldChange > log2(foldChange),])
  res.gen.down <- rownames(res[res$padj < 0.05 & res$log2FoldChange < log2(foldChange),])
  
  # Prepare list for GSEA
  res.lfc <- setNames(object = res$log2FoldChange, rownames(res))
  res.lfc <- sort(res.lfc, decreasing = TRUE)
  
  # Run GSEA (with user-provided TERM2GENE and TERM2NAME in the global environment)
  gsea.result <- GSEA(
    res.lfc, 
    TERM2GENE = TERM2GENE, 
    TERM2NAME = TERM2NAME, 
    eps = 0, 
    nPermSimple = 10000, 
    seed = TRUE
  )
  
  # Create GSEA enrichment dot-plot
  res.enrichment <- enrich_plot(gsea.result, name)
  
  # Combine volcano & enrichment plots side-by-side
  vplots <- grid.arrange(
    res.vplot, 
    res.enrichment, 
    nrow = 1, 
    top = paste("Differential expression analysis for", organism)
  )
  
  # Save combined plot
  ggsave(
    plot = vplots, 
    width = 14, 
    height = 6, 
    units = "in",
    filename = file.path(output.path, paste("dea", paste(tolower(contrast[-1]), collapse = "-"), "svg", sep = "."))
  )
  
  return(list(res.gen.up, res.gen.down))
}

###############################################################################
# Function: plot_venn
# -------------------
# Description:
#   Plots a pairwise Venn diagram for two sets of DEGs (upregulated genes in 
#   two different contrasts). Saves the Venn diagram as an SVG. Also prints 
#   and writes the intersection of the two sets to a text file.
#
# Parameters:
#   output.path : Path to save the Venn diagram and intersection file.
#   g1, g2      : Character vectors of gene IDs from two different contrasts.
#   zone        : Root zone name for labeling (e.g., "MZ").
#   n1, n2      : Labels for the two sets (often something like "MZ-DZ", "MZ-EZ").
#   organism    : Organism name, used in diagram title.
###############################################################################
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
    lty = "blank", cat.cex = 0.9, scaled = FALSE, cat.pos = 180
  )
  
  # Display the Venn diagram plus a title describing the union's size
  grid.arrange(
    gTree(children = venn.plot), 
    top = paste("DEA:", zone, "up-regulated DEGs in", organism),
    bottom = paste('Union:', length(union(g1, g2)))
  )
  dev.off()

  # Find intersection, write to file
  degs.up <- intersect(g1, g2)
  writeLines(degs.up, file.path(output.path, paste("upregulated", zone, "DEGs.txt", sep=".")))
}

###############################################################################
# Define and parse command-line arguments
###############################################################################
option_list = list(
    make_option(c("-m", "--metadata"), type = "character", default = NULL,
                help = "Path to metadata file.", metavar = "character"),
    make_option(c("-c", "--counts"), type = "character", default = NULL,
                help = "Path to parsed featurecounts matrix.", metavar = "character"),
    make_option(c("-f", "--foldChange"), type = "double", default = NULL,
                help = "Fold-change threshold (base = 2) to define DEGs.", metavar = "character"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Output path for results.", metavar = "character"),
    make_option(c("-g", "--term2gene"), type = "character", default = NULL,
                help = "Path to term2gene table file (GO annotations).", metavar = "character"),
    make_option(c("-d", "--term2name"), type = "character", default = NULL,
                help = "Path to term2name table file (GO descriptions).", metavar = "character"),
    make_option(c("-u", "--organism"), type = "character", default = NULL,
                help = "Organism name for labeling.", metavar = "character"),
    make_option(c("-p", "--plts"), type = "character", default = NULL,
                help = "Path to PLT gene table.", metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

###############################################################################
# (Optional overrides for local testing):
# opt$counts <- "output/Arabidopsis_thaliana/FEATURECOUNTS/TAIR10.gene.counts.tsv"
# opt$metadata <- "temp/Arabidopsis_thaliana/metadata.tsv"
# opt$foldChange <- 1
# opt$term2gene <- "input/Arabidopsis_thaliana/GENOME/TERM2GENE.tsv"
# opt$term2name <- "workflow/TERM2NAME.tsv"
# opt$output <- "output/Arabidopsis_thaliana/temp"
# opt$organism <- "Arabidopsis_thaliana"
# opt$plts <- "others/plt.ids.txt"
###############################################################################

###############################################################################
# Step 1: Read input data
###############################################################################
count.matrix <- read.delim(opt$counts, check.names = FALSE)
metadata <- read.delim(opt$metadata)
TERM2GENE <- read.delim(opt$term2gene)
TERM2NAME <- read.delim(opt$term2name)

# Filter TERM2GENE to only those GO terms found in TERM2NAME
TERM2GENE <- TERM2GENE[TERM2GENE$go.id %in% TERM2NAME$go.id, ]

# Create DESeq2 dataset and run DESeq
dds <- DESeqDataSetFromMatrix(countData = count.matrix, colData = metadata, design = ~dex)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)

# Read table of PLT genes
plts <- read.delim(opt$plts, head = FALSE, col.names = c("name", "id"))

###############################################################################
# Step 2: Define contrasts (root zones) and run DEA
###############################################################################
contrasts <- list(
  c("dex", "MZ", "DZ"),
  c("dex", "MZ", "EZ"),
  c("dex", "DZ", "EZ")
)

dea.results <- list()

# For each contrast (e.g., MZ vs DZ), run differential expression analysis
for(i in seq_along(contrasts)) {
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

###############################################################################
# Step 3: Build Venn diagrams for each zone (MZ, DZ, EZ)
###############################################################################
zones <- c("MZ", "DZ", "EZ")

for(zone in zones){
  # Identify which contrasts involve this zone (returns T/F for each contrast)
  c <- sapply(contrasts, function(x) zone %in% x)
  idx <- which(c)
  
  # zone.genes gathers upregulated genes for the relevant contrasts
  zone.genes <- lapply(idx, function(i) dea.results[[i]][[which(zone == contrasts[[i]][-1])]])
  
  # Extract short labels for these contrasts
  c.names <- sapply(contrasts[c], function(x) paste(x[-1], collapse="-"))
  
  # Plot and save Venn diagram for the intersection of upregulated genes in the zone
  plot_venn(opt$output, zone.genes[[1]], zone.genes[[2]], zone, c.names[1], c.names[2], opt$organism)
}
