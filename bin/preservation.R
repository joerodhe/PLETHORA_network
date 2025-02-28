###############################################################################
# Load external function for circle plots and suppress library startup messages
###############################################################################
source("/home/ibtuser/Documentos/PLTs/ara/others/preservation.circleplot.R")
suppressPackageStartupMessages(library(WGCNA))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpattern))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.At.tair.db))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(RColorBrewer))
options(stringsAsFactors = FALSE)

###############################################################################
# Enable multi-threading for WGCNA (up to 30 threads)
###############################################################################
enableWGCNAThreads(30)

###############################################################################
# Function: enrich_plot
# ---------------------
# Description:
#   Takes the results of an over-representation analysis (ORA) and creates a
#   ggplot dot-plot showing log2foldEnrichment vs. term (Description).
###############################################################################
enrich_plot <- function(ora.result, comparison){
  ora.result <- as.data.frame(ora.result)
  plot <- ggplot()
  if(nrow(ora.result) > 0){
    # Limit the maximum number of terms plotted
    if(nrow(ora.result) > 20) max.n <- 20
    else max.n <- nrow(ora.result)
    # Select top terms by adjusted p-value
    ora.result <- head(ora.result[order(ora.result$p.adjust),],max.n)

    # Build the dot-plot
    plot <- ggplot(ora.result, aes(log2foldEnrichment, fct_reorder(Description, log2foldEnrichment))) +
      geom_segment(aes(xend=0, yend = Description)) +
      geom_point(aes(color=p.adjust, size = Count)) +
      # Use a custom color gradient and color bar
      scale_color_gradientn(
        colours=brewer.pal(3, name = "Dark2"),
        guide=guide_colorbar(reverse=TRUE, order=1),
        name = "Adjusted P-value"
      ) +
      scale_size_continuous(range=c(2, 10)) +
      theme_dose(12) +
      xlab("log2foldEnrichment") +
      ylab(NULL) +
      ggtitle(paste0("")) +
      scale_y_discrete(labels = function(x) str_wrap(x, 45)) +
      theme(axis.text.y = element_text(size = 10))
  }
  return(plot)
}

###############################################################################
# Function: plot_pca
# ------------------
# Description:
#   Performs PCA on expression data, then calls 'plotpc' (custom) to visualize
#   principal components using metadata for labeling.
###############################################################################
plot_pca <- function(exp, metadata){
    # Perform PCA on the transposed matrix (samples as rows)
    pca <- prcomp(t(exp))
    # Extract explained variance per principal component
    pca.exp.var <- expVar(pca)
    # plotpc is presumably a user-defined function that produces a PCA plot
    plotpc.object <- plotpc(pca, pca.exp.var, metadata)
    return(plotpc.object)
}

###############################################################################
# Specify paths to expression matrices for different species
###############################################################################
expression.matrices <- c(
    "output/Arabidopsis_thaliana/FEATURECOUNTS/TAIR10.gene.counts.tsv",
    "output/Cucumis_sativus/FEATURECOUNTS/Gy14.gene.counts.tsv",
    "output/Glycine_max/FEATURECOUNTS/Wm82.gene.counts.tsv",
    "output/Oryza_sativa_Japonica_Group/FEATURECOUNTS/MSUv7.gene.counts.tsv",
    "output/Solanum_lycopersicum/FEATURECOUNTS/ITAG4.gene.counts.tsv",
    "output/Zea_mays/FEATURECOUNTS/Zm-B73-REFERENCE-NAM-5.0.55.gene.counts.tsv"
)

###############################################################################
# Associate each expression file with the organism name
###############################################################################
organisms <- c(
    "Arabidopsis_thaliana",
    "Cucumis_sativus",
    "Glycine_max",
    "Oryza_sativa_Japonica_Group",
    "Solanum_lycopersicum",
    "Zea_mays"
)

###############################################################################
# Define the reference organism for ortholog-based mappings
###############################################################################
organism.ref <- "Arabidopsis_thaliana"

# Read a table of orthologs (genes as rows, columns as organisms)
orthologs <- read.delim("output/interspecie/orthologs.bbh.all.tsv", row.names = 1)
nrow(orthologs)

###############################################################################
# Read a reference network file (e.g., from previous analyses) for Arabidopsis
###############################################################################
ref.net <- read.csv("others/prop_table.csv")
nrow(ref.net)

# Filter the reference network genes to only those present in the ortholog table
ref.net.genes <- ref.net$gen_id[which(ref.net$gen_id %in% orthologs[,organism.ref])]
length(ref.net.genes)

###############################################################################
# Read/merge expression data for all species via a custom function
###############################################################################
expression <- read_expression(expression.matrices, organisms, orthologs)
dim(expression)
expression[1:6,1:5]

###############################################################################
# Build a metadata data frame that splits sample row names into 'root_zone' and 'organism'
###############################################################################
metadata <- data.frame(
    root_zone = sapply(rownames(expression), function(ID) strsplit(ID, '_')[[1]][1]),
    organism = sapply(rownames(expression), function(ID) strsplit(ID, '_')[[1]][2]),
    stringsAsFactors = TRUE
)
head(metadata)

###############################################################################
# Perform rlog normalization (via a user-defined function rlog_norm) 
# and inspect the result
###############################################################################
expression.normalized <- rlog_norm(t(expression), metadata, normalization = "rlog")
dim(expression.normalized)
expression.normalized[1:5,1:5]

# Plot PCA on the rlog-transformed data
plot_pca(expression.normalized, metadata)
gsgave("/home/ibtuser/Documentos/PLTs/ara/others/pca.expression.rlog.orthologs.bbh.all.svg")

###############################################################################
# Remove batch effects by organism using limma::removeBatchEffect
###############################################################################
expression.corrected <- removeBatchEffect(expression.normalized, batch = metadata$organism)
dim(expression.corrected)
expression.corrected[1:5,1:5]

# Plot PCA after batch correction
plot_pca(expression.corrected, metadata)
gsgave("/home/ibtuser/Documentos/PLTs/ara/others/pca.expression.rlog.corrected.organism.orthologs.bbg.all.svg")

###############################################################################
# Prepare a module assignment, marking some nodes as reference network nodes
###############################################################################
nodes <- rownames(expression.corrected)  # A vector of all gene IDs
ref.net.nodes <- rownames(orthologs[orthologs[,organism.ref] %in% ref.net.genes,])
length(ref.net.nodes)
ref.net.nodes <- ref.net.nodes[ref.net.nodes %in% nodes]
length(ref.net.nodes)

# Arbitrarily assign modules (0..7) to all nodes, then assign module 8 to reference nodes
nodes.modules <- setNames(
    object = sample(c(0:7), length(nodes), replace = TRUE),
    nodes
)
nodes.modules[ref.net.nodes] <- 8
table(nodes.modules)

###############################################################################
# Split the samples into 18 groups and prepare multiExpr/moduleList for WGCNA
###############################################################################
samples_list <- split(rownames(metadata), cut(seq_along(rownames(metadata)),18,labels=FALSE))
length(samples_list)

# Create empty objects that will store expression data and modules for each subset
multiExpr <- list()
moduleList <- list()
sampleLabels = c()
idx <- 1

# Populate multiExpr/moduleList for each subset
for(i in samples_list){
    sample.expression <- t(expression.corrected[,i])
    multiExpr[[idx]] <- list(data = as.data.frame(sample.expression))
    names(multiExpr[[idx]]$data) <- rownames(expression.corrected)
    sampleLabel <- sapply(rownames(sample.expression), function(x) strsplit(x, "_SR")[[1]][1])
    sampleLabels <- c(sampleLabels, unique(sampleLabel))
    moduleList[[idx]] <- nodes.modules
    idx <- idx+1
}

length(multiExpr)
names(multiExpr) <- sampleLabels
names(moduleList) <- sampleLabels

# Quick checks: dimensions, module counts, and example data slices
lapply(multiExpr, lapply, dim)
lapply(moduleList, length)
lapply(moduleList, table)
multiExpr[[1]]$data[,1:5]
multiExpr[[2]]$data[,1:5]
multiExpr[[3]]$data[,1:5]

###############################################################################
# Check for any bad genes/samples across all sets before modulePreservation
###############################################################################
gsg <- goodSamplesGenesMS(multiExpr, minNSamples = 1, minNGenes = 1,  verbose = 3)
exprSize <- checkSets(multiExpr)
if (!gsg$allOK){
    # Print information about the removed genes:
    if (sum(!gsg$goodGenes) > 0)
        printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes],collapse = ", ")))
    for (set in 1:exprSize$nSets){
        if (sum(!gsg$goodSamples[[set]]))
            printFlush(paste("In set", setLabels[set], "removing samples",
                paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
        # Remove offending genes and samples
        multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes]
        moduleList[[set]] = moduleList[[set]][gsg$goodGenes]
    }
    # Update exprSize
    exprSize = checkSets(multiExpr)
}

###############################################################################
# Perform module preservation analysis, referencing the first set's network
###############################################################################
mp <- modulePreservation(
    multiExpr, 
    moduleList,
    referenceNetworks = 1,
    loadPermutedStatistics = FALSE,
    corFnc = "bicor",
    nPermutations = 100,
    verbose = 4,
    maxGoldModuleSize = 1000, 
    maxModuleSize = 1000,
    parallelCalculation = TRUE, 
    checkData = FALSE,
    savePermutedStatistics = FALSE
)

###############################################################################
# Example: Identify a "PLT" module (8) and check its preservation across sets
###############################################################################
plts.ids <- read.table("input/Arabidopsis_thaliana/PLTs/plt.ids.txt", col.names = c("name", "id"))
plts.ids$node <- sapply(plts.ids$id, function(plt) rownames(orthologs[orthologs[[organism.ref]] == plt, ]))

ref <- 1
preservation.results <- mp$preservation$Z[[1]][-ref]
names(preservation.results) <- sub(
    names(preservation.results),
    pattern = "inColumnsAlsoPresentIn.",
    replacement = "",
    fixed = TRUE
)

plt.module <- "8"

# Build a data frame (Zs) with statistics for each set
Zs <- do.call(
    rbind,
    lapply(
        preservation.results, function(result)
            data.frame(
                Zsummary = result[plt.module,]$Zsummary.pres,
                Zdensity = result[plt.module,]$Zdensity.pres,
                Zconn = result[plt.module,]$Zconnectivity.pres
            )
    )
)

# Extract root_zone, organism, and mark monocots vs. dicots
Zs$Root_zone <- sapply(rownames(Zs), function(x) strsplit(x, "_")[[1]][1])
Zs$Organism <- sapply(rownames(Zs), function(x) strsplit(x, "_")[[1]][2])
monocots <- c("Oryza", "Zea")
Zs$class <- ifelse(Zs$Organism %in% monocots, "Monocots", "Dicots")

###############################################################################
# Plot a bar chart of Zsummary across sets, highlighting thresholds
###############################################################################
ggplot(Zs, aes(x=Root_zone, y=Zsummary, fill=Organism, pattern = class)) +
    geom_bar_pattern(
        stat="identity", 
        position = position_dodge(preserve = "single"),
        pattern_fill = "black",
        pattern_density = 0.1,
        pattern_spacing = 0.025,
        pattern_angle = 45,
        pattern_key_scale_factor = 0.6,
        color = "black"
    ) +
    scale_fill_brewer(palette = "Dark2") +
    scale_pattern_manual(values = c(Monocots = "circle", Dicots = "none")) +
    geom_hline(yintercept = 10, linetype = "dashed", color="red") +
    geom_hline(yintercept = 2, linetype = "dashed", color="orange") +
    labs(title = "Zsummary. Reference: Arabidopsis MZ network", pattern = "") +
    guides(
        pattern = guide_legend(override.aes = list(fill = "white")),
        fill = guide_legend(override.aes = list(pattern = "none"))
    ) +
    theme_classic()

# Save the bar plot
ggsave("others/zsummary.barplot.pdf")

###############################################################################
# KME scatterplots (example code, partial)
###############################################################################
par(mar = c(3.3, 3.3, 4, 0.5))
par(mgp = c(1.9, 0.6, 0))
ref = 1
ind = 5
for (set in 1:nSets){
    if (set!=ref){
        verboseScatterplot(KMEpathway[, ref], KMEpathway[, set],
            xlab = spaste("KME ", colnames(KMEpathway)[set]),
            ylab = "KME MZ_Arabidopsis",  abline = TRUE)
        # main = spaste(LETTERS[ind], ". KME in ", organisms[order2[set]], "\nvs. ",
        # organisms[ref], "\n"),
        # cex.main = 1.4, cex.lab = 1.2, cex.axis = 1.2
        ind = ind + 1;
    } else {
        plot(c(0,1), type ="n", axes = FALSE, xlab = "", ylab ="")
    }
}
# Close the plotting device
dev.off()

###############################################################################
# GO enrichment of the reference network genes (example)
###############################################################################
ref.net.genes.ora <- enrichGO(
    gene = ref.net.genes, 
    keyType = 'TAIR', 
    OrgDb = org.At.tair.db, 
    ont = "All", 
    pAdjustMethod = "BH", 
    pvalueCutoff = 0.05, 
    qvalueCutoff = 0.05
)

# Calculate log2(Fold-Enrichment)
ref.net.genes.ora <- mutate(
    ref.net.genes.ora, 
    log2foldEnrichment = log2(
        (as.numeric(sub("/\\d+", "", GeneRatio)) / as.numeric(sub("\\d+/", "", GeneRatio))) /
        (as.numeric(sub("/\\d+", "", BgRatio))   / as.numeric(sub("\\d+/", "", BgRatio)))
    )
)

ref.net.genes.ora <- as.data.frame(ref.net.genes.ora)
head(ref.net.genes.ora)
ref.net.genes.ora <- ref.net.genes.ora[order(ref.net.genes.ora$log2foldEnrichment, decreasing = TRUE), ]
head(ref.net.genes.ora[,c("Description", "Count"), drop = FALSE], 10)

# Extract genes for the top GO terms
go.genes <- strsplit(ref.net.genes.ora$geneID[1], "/")[[1]]
go.genes <- lapply(ref.net.genes.ora$geneID[1:10], function(x) strsplit(x, "/")[[1]])
go.genes <- unique(unlist(go.genes))
length(go.genes)

###############################################################################
# Example: Searching "RNA polymerase" terms in a user-defined object 'terms2'
###############################################################################
candidates = grep("RNA polymerase", terms2, fixed = TRUE)
length(candidates)
go.terms <- terms2[candidates]

writeLines(paste(names(go.terms), go.terms), "go.terms.txt")

go.ids <- names(terms2[candidates])

# Build a list of GO offspring across CC, BP, MF
go.offspring <- c(
    as.list(GOCCOFFSPRING),
    as.list(GOBPOFFSPRING),
    as.list(GOMFOFFSPRING)
)
length(go.offspring)
head(go.offspring)

go.offspring.ids <- names(go.offspring)

# Build a parent-child table for the selected GO IDs
go.relations <- do.call(rbind,
    lapply(go.ids, function(go.id){
        child.terms <- go.offspring[go.offspring.ids == go.id]
        x <- data.frame(
            parent = rep(go.id, length(child.terms)),
            childterm = child.terms[[1]]
        )
        return(x)
    })
)

go.ids.offspring <- unique(c(
    unique(go.relations$parent),
    as.character(na.omit(go.relations$childterm))
))
length(go.ids.offspring)

# Read general GO annotations (TERM2GENE format)
go.anns <- read.delim("input/Arabidopsis_thaliana/GENOME/TERM2GENE.tsv")
nrow(go.anns)
head(go.anns)

# Gather all annotated genes for those parent+child GO IDs
go.genes <- data.frame()
for(i in 1:length(go.ids.offspring)){
    go.genes <- rbind(
        go.genes,
        go.anns[which(go.anns$go.id == go.ids.offspring[i]),]
    )
}
nrow(go.genes)

###############################################################################
# Intersect reference network genes with these GO terms
###############################################################################
ref.net.genes.go <- ref.net.genes[ref.net.genes %in% go.genes]
length(ref.net.genes.go)
writeLines(ref.net.genes.go)
ref.net.nodes.go <- rownames(orthologs[orthologs[,organism.ref] %in% ref.net.genes.go,])
length(ref.net.nodes.go)

# Only keep nodes present in the expression data
ref.net.nodes.go <- ref.net.nodes.go[ref.net.nodes.go %in% colnames(multiExpr[[1]]$data)]
length(ref.net.nodes.go)

###############################################################################
# Aggregate expression within each set, then build correlation for ref.net.nodes
###############################################################################
expr <- list()
for(set in 1:length(multiExpr)){
    expr[[set]] <- list(data = t(apply(multiExpr[[set]]$data, 1, tapply, colnames(multiExpr[[set]]$data), mean, na.rm = TRUE)))
}

###############################################################################
# Function: build_cor
# -------------------
# Description:
#   Calculates bicor on a matrix, replaces diagonals and upper triangle with NA,
#   then melts to a data frame.
###############################################################################
build_cor <- function(x){
    x <- bicor(x, nThreads = 16)
    diag(x) <- NA
    x[upper.tri(x)] <- NA
    x.melt <- melt(x, na.rm = TRUE)
    return(x.melt)
}

###############################################################################
# Combine correlation data from each set into a single data frame
###############################################################################
x <-do.call(cbind, lapply(expr, function(expr.dataset){
    x <- build_cor(expr.dataset$data[,ref.net.nodes])
    rownames(x) <- paste(x[,1], x[,2], sep = "_")
    x <- x[ , 3, drop = FALSE]
    return(x)
}))

colnames(x) <- names(multiExpr)

###############################################################################
# Example: separate columns into reference vs. non-reference sets
###############################################################################
col.ref <- c(1,4,7,10,13,16)
col.rest <- setdiff(1:18, col.ref)

# Calculate variance of edges across the reference columns
edges.var <- apply(x[,col.ref],1,var)

# Identify edges with the smallest variance
edges.less.var <- sort(edges.var, decreasing = FALSE)[1:100]
nodes.less.var <- unique(unlist(strsplit(names(edges.less.var), "_")))
length(nodes.less.var)

###############################################################################
# Map low-variance node IDs back to reference organism gene IDs
###############################################################################
genes.less.var <- orthologs[rownames(orthologs) %in% nodes.less.var,1]
writeLines(genes.less.var)

###############################################################################
# Enrichment analysis for these genes, restricting GO to BP, with ref.net as universe
###############################################################################
genes.less.var.ora <- enrichGO(
    gene = genes.less.var, 
    keyType = 'TAIR', 
    OrgDb = org.At.tair.db, 
    ont = "BP", 
    universe = ref.net$gen_id,
    pAdjustMethod = "BH", 
    pvalueCutoff = 0.05, 
    qvalueCutoff = 0.02
)

genes.less.var.ora <- mutate(
    genes.less.var.ora, 
    log2foldEnrichment = log2(
      (as.numeric(sub("/\\d+", "", GeneRatio)) / as.numeric(sub("\\d+/", "", GeneRatio))) /
      (as.numeric(sub("/\\d+", "", BgRatio))   / as.numeric(sub("\\d+/", "", BgRatio)))
    )
)

enrich_plot(genes.less.var.ora, "Down-ora")
ggsave("tmo.ora.pdf", width = 7, height = 4)

###############################################################################
# Build adjacency-like matrices for the selected nodes and visualize distributions
###############################################################################
pathwayAdjs <- list()
pathwayAdjs.pval <- list()
for (set in 1:length(multiExpr)){
    printFlush(paste("Working on set", names(multiExpr)[set]))
    bc <- bicorAndPvalue(expr[[set]]$data[, nodes.less.var], use = "p")
    # Weighted adjacency = correlation * sign
    pathwayAdjs[[set]] <- abs(bc$bicor) * sign(bc$bicor)
    pathwayAdjs.pval[[set]] <- bc$p
}

# Compute node degree (sum of absolute adjacency) in the reference set
ref.nodes.degree <- apply(abs(pathwayAdjs[[1]]), 2, sum)-1
ref.nodes.degree <- order(ref.nodes.degree, decreasing = TRUE)

###############################################################################
# Plot distributions of edge weights for each set
###############################################################################
sizeGrWindow(8,24)
par(mfrow =c(6,3))
par(mar = c(2.5, 1, 2.5, 1))
for (set in 1:length(pathwayAdjs)){
    adjs <- pathwayAdjs[[set]]
    adjs[upper.tri(adjs)] <- NA
    diag(adjs) <- NA
    edges <- as.numeric(na.omit(c(adjs)))
    den <- density(edges)
    edges.limit <- max(abs(edges))
    cat("Edge threshold for network", names(multiExpr)[set], "is", edges.limit, "\n")
    plot(den, frame = FALSE, col = "blue", main = names(multiExpr)[set])
}

###############################################################################
# Circle plot visualizations (stored in a PDF)
###############################################################################
pdf(file = "/home/ibtuser/Documentos/PLTs/ara/others/ribosome.biogenesis.1.pdf", wi=10, h=24)
par(mfrow =c(6,3))
par(mar = c(0.3, 1, 1.5, 1))
for (set in 1:length(multiExpr)){
    circlePlot(
        pathwayAdjs[[set]], 
        pathwayAdjs.pval[[set]],
        colnames(pathwayAdjs[[set]]), 
        ref.nodes.degree, 
        filterPval = TRUE,
        main = names(multiExpr)[set],
        variableLabelAngle = TRUE, 
        min.cex.labels = 0.7,
        max.cex.labels = 1, 
        radii = c(0.6,0.5),
        center = c(0.1, 0.001), 
        variable.cex.labels = TRUE
    )
}
dev.off()

###############################################################################
# (Commented-out code: alternative or exploratory analyses)
###############################################################################
# multiMEs = multiSetMEs(multiExpr, universalColors = nodes.modules)
#
# mz.nodes.plt.regulated <- mz.nodes.plt.regulated[mz.nodes.plt.regulated %in% names(moduleList[[set]])]
# ref.net.nodes <- ref.net.nodes[ref.net.nodes %in% names(moduleList[[set]])]
# nGenes <- length(ref.net.nodes)
#
# nSets <- length(multiExpr)
#
# KMEpathway = matrix(0, nGenes, nSets)
# for (set in 1:nSets){
#     KMEpathway[, set] <- cor(
#         multiExpr[[set]]$data[,ref.net.nodes], 
#         multiMEs[[set]]$data[, 3], 
#         use = "p")
# }
#
# row.names(KMEpathway) <- ref.net.nodes
# colnames(KMEpathway) <- names(multiExpr)
#
# w <- rcorr(as.matrix(t(x[,col.ref])), type="pearson")
# w.pval <- w$P
# diag(w.pval) <- NA
# w.pval[upper.tri(w.pval)] <- NA
# edges.preservation.pval <- melt(w.pval[1:100,1:100], na.rm = TRUE)
# edges.preserved <- edges.preservation.pval[edges.preservation.pval$value < 0.05,]
#
# pvals <- c()
# for(i in 1:nrow(x)){
#     w <- wilcox.test(as.numeric(x[i,col.ref]), as.numeric(x[i,col.rest]))
#     pvals <- c(pvals, w$p.value)
# }
#
# names(pvals) <- rownames(x)
#
# pvals.adjusted <- p.adjust(pvals, method = "BH")
#
# edges.diff <- names(pvals[pvals < 0.05])
# nodes.diff <- strsplit(edges.diff, "_")
