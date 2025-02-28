suppressPackageStartupMessages(library(WGCNA))
suppressPackageStartupMessages(library(flashClust))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(e1071))
suppressPackageStartupMessages(library(mixtools))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(showtext))
suppressPackageStartupMessages(library(parallel))

# Set global options
options(stringsAsFactors = FALSE)
showtext_auto()
set.seed(3006)

################################################################################
# Function: checkExpression
# -------------------------
# Description:
#   Ensures that each gene and sample meets the criteria for having enough valid
#   expression data. Removes offending genes or samples based on the WGCNA 
#   function goodSamplesGenes.
#
# Params:
#   expression.matrix: Expression matrix (samples x genes)
#
# Returns:
#   The same matrix (possibly with some rows/columns removed if they fail checks).
################################################################################
checkExpression <- function(expression.matrix){
  # Check which genes and samples are suitable
  gsg = goodSamplesGenes(expression.matrix, verbose = 3)
  
  # If there are genes or samples failing the checks, remove them
  if (!gsg$allOK){
    # Optionally, print the gene names removed
    if (sum(!gsg$goodGenes)>0)
      dynamicTreeCut::printFlush(
        paste("Removing genes:",
              paste(names(expression.matrix)[!gsg$goodGenes], collapse = ", "))
      );
    
    # Optionally, print the sample names removed
    if (sum(!gsg$goodSamples)>0)
      dynamicTreeCut::printFlush(
        paste("Removing samples:",
              paste(rownames(expression.matrix)[!gsg$goodSamples], collapse = ", "))
      );
    
    # Actually remove the offending genes/samples from the data
    expression.matrix = expression.matrix[gsg$goodSamples, gsg$goodGenes]
  }
  return(expression.matrix)
}

################################################################################
# Function: module.enrichment
# ---------------------------
# Description:
#   Performs an Over-Representation Analysis (ORA) for the genes belonging
#   to a specific module, using user-provided TERM2GENE and TERM2NAME tables.
#   Plots an enrichment dot-plot if there are enriched terms.
#
# Params:
#   module      : The module color or identifier.
#   TERM2GENE   : Data frame mapping each GO term to its associated genes.
#   TERM2NAME   : Data frame mapping GO terms to their descriptions.
#   background  : Vector of background genes for ORA.
#   output.path : Output folder where results and plots are stored.
#
# Returns:
#   A ggplot object containing the enrichment plot (or an empty ggplot if none).
################################################################################
module.enrichment <- function(module, TERM2GENE, TERM2NAME, background, output.path){

  cat("Performing Overrepresentantion functional enrichment over", module, "module...\n")

  # Default empty plot
  plot <- ggplot()
  
  # Check if the module node file is present
  if(file.exists(file.path(output.path, paste0(module, '.nodes.tsv')))){

    # Read the list of genes (nodes) for the current module
    module.nodes <- read.delim(file.path(output.path, paste0(module, '.nodes.tsv')))
    query <- module.nodes[,1]

    # Filter the module genes to ensure they are also in the background
    query <- query[which(query %in% background)]

    cat("Number of genes in background:", length(background), "\n")
    cat("Number of genes in query (after matching with background):", length(query), "\n")

    # Perform enrichment only if there are genes overlapping with background
    if (length(query) > 0) {
      # Enricher function from clusterProfiler
      ora.results <- enricher(
        gene = query,
        universe = background,
        TERM2GENE = TERM2GENE,
        TERM2NAME = TERM2NAME
      )
      
      cat("Number of terms enriched:", nrow(ora.results), "\n")

      # Only plot if there are enriched terms
      if (!is.null(ora.results)) if (nrow(ora.results) > 0) {

        # Compute log2(Fold-Enrichment) for each term
        ora.results <- mutate(
          ora.results, 
          log2foldEnrichment = log2(
            (as.numeric(sub("/\\d+", "", GeneRatio)) / as.numeric(sub("\\d+/", "", GeneRatio))) /
            (as.numeric(sub("/\\d+", "", BgRatio))   / as.numeric(sub("\\d+/", "", BgRatio)))
          )
        )

        # Create a dot-plot showing enriched terms
        plot <- ggplot(
          ora.results, 
          showCategory = 10,
          aes(log2foldEnrichment, fct_reorder(Description, log2foldEnrichment))
        ) +
          geom_segment(aes(xend=0, yend = Description)) +
          geom_point(aes(color=p.adjust, size = Count)) +
          scale_color_gradientn(
            colours=c("#f7ca64", "#46bac2", "#7e62a3"),
            trans = "log10",
            guide=guide_colorbar(reverse=TRUE, order=1)
          ) +
          scale_size_continuous(range=c(2, 10)) +
          theme_dose(12) +
          xlab("Log2(Fold-Enrichment)") +
          ylab(NULL) +
          ggtitle(paste("GO enrichment for", module, "module")) +
          scale_y_discrete(labels = function(x) str_wrap(x, 45)) +
          theme(axis.text.y = element_text(size = 10))
      }
    }
  }
  return(plot)
}

################################################################################
# Function: decompose_gmm
# -----------------------
# Description:
#   Decomposes a numeric vector (edge weights) into two Gaussian distributions
#   using the mixtools package (normalmixEM). Also identifies and handles outliers.
#
# Params:
#   edge.weight: A named numeric vector of edge weights from TOM or adjacency.
#
# Returns:
#   A list containing the mixture means, lambdas, sigmas, and the 
#   probabilities data frame with posterior probabilities of belonging to 
#   each component.
################################################################################
decompose_gmm <- function(edge.weight){
  # Identify outliers using boxplot statistics
  stats <- boxplot.stats(edge.weight)
  outliers <- which(edge.weight %in% stats$out)
  noutliers <- length(outliers)
  cat('Number of outliers:', noutliers, '\n')
  if( noutliers == 0) outliers <- length(edge.weight)+1

  # Decompose Gaussian distribution into 2 components
  maxit = 10000
  cat('Performing GMM with', maxit, 'max iterations...\n')
  gaussian.dist <- normalmixEM(edge.weight[-outliers], k=2, maxit=maxit)
  cat('Done...\n')

  # Build a posterior probability data frame for each observation
  probs <- as.data.frame(
    gaussian.dist$posterior,
    row.names = names(edge.weight)[-outliers]
  )

  # Reintroduce outliers with artificially forced probabilities of 0/1
  if (noutliers > 0) {
    cat('Adding outliers...\n')
    probs.outliers <- do.call(rbind, lapply(edge.weight[outliers], function(outlier){
      return(data.frame(
        'comp.1' = ifelse(outlier<=stats$stats[1], 1, 0),
        'comp.2' = ifelse(outlier>=stats$stats[5], 1, 0)
      ))
    }))
    rownames(probs.outliers) <- names(edge.weight)[outliers]
    probs <- rbind(probs, probs.outliers)
    edge.weight <- c(edge.weight[-outliers], edge.weight[outliers])
  }

  # Combine edge weight and the posterior probabilities
  cat('Building dataframe...\n')
  decomposed.dist <- cbind(edge.weight, probs)

  return(list(
    mu=gaussian.dist$mu, 
    lambda=gaussian.dist$lambda,
    sigma=gaussian.dist$sigma, 
    probabilities=decomposed.dist
  ))
}

################################################################################
# Function: plot_group_proportion
# -------------------------------
# Description:
#   Plots the proportion (count) of each identified group (e.g., Connected vs 
#   NotConnected edges) to get a quick bar chart of how many edges belong to each 
#   mixture group.
#
# Params:
#   dist.dat     : Data frame from the mixture decomposition containing group info.
#   observation.type : A label for the y-axis (e.g. "edges", "nodes").
#   plot.title   : Plot title.
#
# Returns:
#   A ggplot object (bar plot).
################################################################################
plot_group_proportion <- function(dist.dat, observation.type, plot.title){
  plot <- ggplot(dist.dat, aes(x=factor(group))) +
    geom_bar() +
    xlab(element_blank()) +
    ylab(paste("Number of", observation.type)) +
    labs(title=plot.title) +
    theme_classic()
  return(plot)
}

################################################################################
# Function: plot_dist
# -------------------
# Description:
#   Plots a histogram of the edge weights and overlays the two fitted Gaussian
#   curves from the mixture model.
#
# Params:
#   dist.dat        : A list from decompose_gmm containing mixture parameters 
#                     and probabilities.
#   plot.title      : A title for the histogram.
#   connected.group : Which component is considered "connected" (comp.1 or comp.2).
#
# Returns:
#   A ggplot object showing the distribution of edge weights plus the 
#   component curves in different colors.
################################################################################
plot_dist <- function(dist.dat, plot.title, connected.group){

  # Convert "comp.1"/"comp.2" labels into numeric indices to pick the right params
  connected.group <- ifelse(connected.group == 'comp.1', 1, 2)
  unconnected.group <- ifelse(connected.group == 1, 2, 1)

  # Build a histogram of edge weights
  plot <- data.frame(x = dist.dat$probabilities$edge.weight) %>%
    ggplot() +
    geom_histogram(
      aes(x, after_stat(density)), 
      binwidth = 0.05, colour = "black", fill = "white"
    ) +
    # Add the first (connected) curve in red
    stat_function(
      geom = "line", 
      fun = plot_mix_comps,
      args = list(dist.dat$mu[connected.group], 
                  dist.dat$sigma[connected.group], 
                  lam = dist.dat$lambda[connected.group]),
      colour = "red", 
      lwd = 1
    ) +
    # Add the second (unconnected) curve in blue
    stat_function(
      geom = "line", 
      fun = plot_mix_comps,
      args = list(dist.dat$mu[unconnected.group], 
                  dist.dat$sigma[unconnected.group], 
                  lam = dist.dat$lambda[unconnected.group]),
      colour = "blue", 
      lwd = 1
    ) +
    ylab("Density") + 
    xlab("Topological Overlap Value") +
    labs(title=plot.title) +
    theme_classic()

  return(plot)
}

################################################################################
# Function: plot_mix_comps
# ------------------------
# Description:
#   A helper function for stat_function() calls in plot_dist(), which draws the
#   Gaussian components. 
#
# Params:
#   x     : x data for stat_function.
#   mu    : Mean of the Gaussian component.
#   sigma : Standard deviation of the Gaussian component.
#   lam   : Mixture proportion λ.
#
# Returns:
#   The probability density for the given x, scaled by lam.
################################################################################
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

################################################################################
# Function: module.export
# -----------------------
# Description:
#   Exports edges/nodes for a given module as TSV files (for Cytoscape), 
#   decomposes its TOM edges via a GMM approach to determine "connected" edges,
#   identifies hub genes based on kME, and plots distribution figures.
#
# Params:
#   module        : Module color or label (e.g., "turquoise").
#   dynamicColors : A named vector of module assignments for each gene.
#   gene.names    : Vector of gene names in the data.
#   TOM           : The topological overlap matrix for all genes.
#   temp.dir      : Path for temporary output files.
#   out.dir       : Path for final output files.
#   org.name      : Name of the organism (for figure titles).
#   datKME        : Data frame of signed kME values for each gene x module.
#   meristem.nodes: A vector of special genes (e.g., from a background list).
#
# Returns:
#   (No direct return). Writes multiple files and figures as a side effect.
################################################################################
module.export <- function(module,
                          dynamicColors, 
                          gene.names, 
                          TOM, 
                          temp.dir, 
                          out.dir, 
                          org.name, 
                          datKME,
                          meristem.nodes){
  # Print status
  cat('Decomposing', module, 'edge weights...\n')

  # Identify genes belonging to this module
  inModule <- is.finite(match(dynamicColors,module))
  modGenes <- gene.names[inModule]

  # Subset TOM to only those genes
  modTOM <- TOM[inModule,inModule]
  dimnames(modTOM) <- list(modGenes,modGenes)

  # Export module network to Cytoscape edge and node files in a temporary folder
  cyt <- exportNetworkToCytoscape(
    modTOM, 
    edgeFile = file.path(temp.dir, paste0(module, '.edges.tsv')), 
    nodeFile = file.path(temp.dir, paste0(module, '.nodes.tsv')), 
    weighted = TRUE, 
    threshold = 0, 
    nodeNames = modGenes
  )

  # Read edges from the just-created file
  mod.edges <- read.delim(file.path(temp.dir, paste0(module, '.edges.tsv')))

  # Create a named vector of edge weights
  mod.weights <- setNames(
    object = mod.edges$weight,
    rownames(mod.edges)
  )

  # Decompose edge weight distribution with GMM
  TOM.GMM <- decompose_gmm(mod.weights)

  # Identify which Gaussian component is "connected" or "not connected"
  cat('Spliting groups based on probability threshold...\n')
  weight.comp.1 <- TOM.GMM$probabilities$edge.weight[TOM.GMM$probabilities$comp.1 >= 0.5]
  weight.comp.1.median <- median(weight.comp.1)
  cat('Median edge weight of distribution 1 in', module, 'module is', weight.comp.1.median, '\n')
  
  weight.comp.2 <- TOM.GMM$probabilities$edge.weight[TOM.GMM$probabilities$comp.2 >= 0.5]
  weight.comp.2.median <- median(weight.comp.2)
  cat('Median edge weight of distribution 2 in', module, 'module is', weight.comp.2.median, '\n')

  # Decide which component is the "connected" group based on medians
  connected.group <- ifelse(
    weight.comp.1.median > weight.comp.2.median,
    'comp.1',
    'comp.2'
  )
  if(is.na(connected.group)) connected.group <- 'comp.2'
  cat('Connected group for', module, 'is', connected.group, '\n')

  # Label edges as "Connected" or "NotConnected" in the probability table
  TOM.GMM$probabilities$group <- ifelse(
    TOM.GMM$probabilities[[connected.group]] >= 0.5,
    "Connected",
    "NotConnected"
  )

  cat('Number of connected nodes in', 
      module, 'module is', nrow(TOM.GMM$probabilities[which(TOM.GMM$probabilities$group == 'Connected'), ]),
      '\n')

  # Write out the GMM results as a TSV
  cat('Writing GMM results...\n')
  write.table(
    as.data.frame(TOM.GMM$probabilities), 
    file = file.path(temp.dir, paste0(module, '.gmm.tsv')),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE
  )

  # Determine the edge weight threshold for "connected"
  edge.threshold <- min(TOM.GMM$probabilities$edge.weight[which(TOM.GMM$probabilities$group == 'Connected')])

  # Plot distribution & group proportions
  cat('Plotting...\n')
  gmm.propotion.plot <- plot_group_proportion(
    TOM.GMM$probabilities,
    observation.type="edges",
    plot.title=paste("Proportion of edges\nMin W:", round(edge.threshold, 2))
  )

  gmm.distribution.plot <- plot_dist(
    TOM.GMM, 
    plot.title='Topological Overlap Matrix gaussian decomposition',
    connected.group=connected.group
  )

  # Combine both plots into a single figure
  gmm.plot <- grid.arrange(
    gmm.propotion.plot,
    gmm.distribution.plot,
    widths=c(0.3, 0.7),
    nrow = 1,
    top = paste("Gaussian Model Mixture of", org.name, "WGNCA", module, "module edge weights")
  )
  
  # Save the combined figure
  ggsave(
    plot = gmm.plot, 
    width = 10, 
    height = 4, 
    file = file.path(out.dir, paste0(module, '.gmm.svg'))
  )

  # Identify hub genes as those above the 90th percentile of kME
  datKME <- datKME[modGenes, module]
  hub.threshold <- quantile(datKME, probs = 0.9)
  cat('Hub threshold for', module, 'is', hub.threshold, '(corresponding to 10% of higher kME)\n')

  nodes.metadata <- data.frame(
    kME = datKME,
    hub = ifelse(datKME > hub.threshold, TRUE, FALSE),
    MZ.expressed = ifelse(modGenes %in% meristem.nodes, TRUE, FALSE),
    row.names = modGenes
  )

  # Export the final network (edges/nodes) above the determined threshold
  cyt = exportNetworkToCytoscape(
    modTOM, 
    edgeFile = file.path(out.dir, paste0(module, '.edges.tsv')), 
    nodeFile = file.path(out.dir, paste0(module, '.nodes.tsv')), 
    weighted = TRUE,
    nodeAttr = nodes.metadata, 
    threshold = edge.threshold, 
    nodeNames = modGenes
  )  
}

################################################################################
# Define command-line options
################################################################################
option_list = list(
  make_option(c("-c", "--expression"), type = "character", default = NULL,
              help = "Path to parsed featurecounts matrix converted to rlog values.", metavar = "character"),
  make_option(c("-s", "--softPower"), type = "integer", default = 6,
              help = "Force to use a fixed softPower.", metavar = "integer"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output path.", metavar = "character"),
  make_option(c("-g", "--term2gene"), type = "character", default = NULL,
              help = "Path to term2gene table file.", metavar = "character"),
  make_option(c("-d", "--term2name"), type = "character", default = NULL,
              help = "Path to term2name table file.", metavar = "character"),
  make_option(c("-b", "--background"), type = "character", default = NULL,
              help = "Background genes for Overrepresentation analysis (meristem DEGs)", metavar = "character"),   
  make_option(c("-t", "--threads"), type = "integer", default = 1,
              help = "Max number of cores to be used in TOM matrix calculation.", metavar = "integer"),
  make_option(c("-m", "--mergeThreshold"), type = "double", default = 1,
              help = "Threshold to merge WGCNA modules (eigengenes).", metavar = "double"),
  make_option(c("-f", "--forceTOMCalc"), type = "character", default = NULL,
              help = "Force re-calculation of TOM matrix.", metavar = "character"),   
  make_option(c("-u", "--temp"), type = "character", default = NULL,
              help = "Temp path", metavar = "character")       
)

# Parse arguments
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# (Below block of opt assignments is for testing/debugging; remove or comment in production)
opt$expression <- "output/Arabidopsis_thaliana/FEATURECOUNTS/TAIR10.gene.rpkm.tsv"
opt$softPower <- 12
opt$term2gene <- "input/Arabidopsis_thaliana/GENOME/TERM2GENE.tsv"
opt$term2name <- "workflow/TERM2NAME.tsv"
opt$output <- "output/Arabidopsis_thaliana/temp2"
opt$background <- "output/Arabidopsis_thaliana/DEA/upregulated.MZ.DEGs.txt"
opt$forceTOMCalc <- "False"
opt$threads <- 3
opt$temp <- 'temp/Arabidopsis_thaliana/WGCNA2'
opt$org.name <- 'Arabidopsis thaliana'
opt$gmm.threshold <- 0.5

cat("Building networks from rlog matrix:", opt$expression, "\n")
cat("Force re-building networks:", opt$forceTOMCalc, "\n")

# Enable multi-threading for WGCNA
enableWGCNAThreads(opt$threads)

# Define path to store intermediate R workspace
temp.WGCNA <- file.path(opt$temp, 'WGCNA.RData')

# If we already have previous data and user wants to reuse it, load it
if(file.exists(temp.WGCNA) & opt$forceTOMCalc == "False"){
  cat("Loading data of previous session...\n")
  load(temp.WGCNA)
} else {
  # Otherwise, read the expression matrix and prepare it
  expression.matrix <- read.delim(opt$expression)
  # Transpose to get samples as rows and genes as columns
  expression.matrix <- t(expression.matrix)
  
  # Clean expression matrix by removing bad samples or genes
  expression.matrix <- checkExpression(expression.matrix)
  gene.names <- colnames(expression.matrix)

  # Check for sample outliers with hierarchical clustering
  sampleTree = hclust(dist(expression.matrix), method = "average")

  # Choose a set of candidate powers
  cat("Analyzing powers...\n")
  powers = c(c(1:20), seq(from = 22, to=30, by=2))
  sft = pickSoftThreshold(expression.matrix, powerVector = powers)

  cat("Calculating TOM similarity...\n")
  # Build adjacency (co-expression network) using the selected softPower
  adjacency = adjacency(expression.matrix, power = opt$softPower)
  # Convert adjacency to topological overlap (TOM)
  TOM = TOMsimilarity(adjacency, verbose = 3)
  rm(adjacency)

  # 1 - TOM dissimilarity
  dissTOM = 1 - TOM
  cat("Running hierarchical clustering over TOM dissimilarity matrix...\n")
  geneTree = hclust(as.dist(dissTOM), method = "average")

  # Set minimum module size
  minModuleSize = 30
  cat("Performing module identification through dynamic tree cut...\n")
  dynamicMods = cutreeDynamic(
    dendro = geneTree, 
    distM = dissTOM,
    deepSplit = 2, 
    pamRespectsDendro = FALSE,
    minClusterSize = minModuleSize
  )
  rm(dissTOM)

  # Convert numeric labels into colors
  dynamicColors = labels2colors(dynamicMods)
  cat("Calculating eigengenes...\n")
  MEList = moduleEigengenes(expression.matrix, colors = dynamicColors)
  MEs = MEList$eigengenes
  rm(MEList)

  # Dissimilarity of module eigengenes
  MEDiss = 1 - cor(MEs)
  rm(MEs)
  
  cat("Performing hierarchical clusting over eigengenes...\n")
  METree = hclust(as.dist(MEDiss), method = "average")

  cat("Merging modules...\n")
  merge = mergeCloseModules(expression.matrix, dynamicColors, 
                            cutHeight = opt$mergeThreshold, verbose = 3)
  mergedColors = merge$colors
  mergedMEs = merge$newMEs
  rm(merge)

  # Save all variables to resume session in the future
  gc()
  cat("Saving session...\n")
  cat('Writing file', temp.WGCNA, '\n')
  save.image(file = temp.WGCNA)
}

# (Re-)Calculate module eigengenes if needed
MEList = moduleEigengenes(expression.matrix, colors = dynamicColors)
datME = MEList$eigengenes

# Calculate kME for each gene
datKME <- signedKME(expression.matrix, datME, outputColumnName='')
datKME.file <- file.path(opt$temp, 'kme.tsv')
cat('Saving KME table file in:', datKME.file, '\n')
write.table(datKME, file = datKME.file, sep = '\t', quote = FALSE)

# Read the background gene list
background <- readLines(opt$background)

# Export each module with dedicated analysis
x <- mclapply(
  unique(dynamicColors),
  module.export,
  dynamicColors = dynamicColors,
  gene.names = gene.names,
  TOM = TOM,
  temp.dir = opt$temp,
  out.dir = opt$output,
  org.name = opt$org.name,
  datKME = datKME,
  meristem.nodes = background,
  mc.cores=opt$threads
)

# Example call for a specific module
module.export(
  "turquoise",
  dynamicColors = dynamicColors,
  gene.names = gene.names,
  TOM = TOM,
  temp.dir = opt$temp,
  out.dir = opt$output,
  org.name = opt$org.name,
  datKME = datKME,
  meristem.nodes = background
)

# Plot sample tree
svg(file.path(opt$output, "sampleTree.svg"))
par(cex = 1.3);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", 
     sub="", xlab="", cex.lab = 1.0,cex.axis = 1.5, cex.main = 1)
dev.off()

# Plot the scale-free topology analysis results
svg(file.path(opt$output, "powers.svg"))
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2", type="n",
     main = paste("Scale independence. Selected threshold:", opt$softPower)
)
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", 
     type="n",
     main = paste("Mean connectivity")
) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# Plot hierarchical clustering of module eigengenes
svg(file.path(opt$output, "eigengenes.svg"))
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
# Red line for the merge threshold
abline(h=opt$mergeThreshold, col = "red")
dev.off()

# Plot dendrogram of genes with dynamic module colors
svg(file.path(opt$output, "dendro.svg"))
plotDendroAndColors(
  geneTree, 
  cbind(dynamicColors, mergedColors),
  c("Dynamic Tree Cut", "Merged dynamic"),
  dendroLabels = FALSE, 
  hang = 0.03,
  addGuide = TRUE, 
  guideHang = 0.05
)
dev.off()

# Read TERM2GENE and TERM2NAME to perform GO-based enrichment analyses
TERM2GENE <- read.delim(opt$term2gene)
TERM2NAME <- read.delim(opt$term2name)
TERM2GENE <- TERM2GENE[TERM2GENE$go.id %in% TERM2NAME$go.id, ]

# Filter background to keep only those in TERM2GENE
background <- background[which(background %in% TERM2GENE$gen.id)]

# For each module, run ORA and save the plot
for (color in unique(dynamicColors)) {
  module.ora <- module.enrichment(
    color,
    TERM2GENE = TERM2GENE, 
    TERM2NAME = TERM2NAME, 
    background = background,
    output.path = opt$output
  )
  ggsave(filename = file.path(opt$output, paste0(color,".ora.svg")), plot = module.ora)
}

# Below is an example of referencing a stored WGCNA session and 
# applying further subsetting or adjacency-based analysis for a reference set of genes

# Load a previously computed WGCNA session (example)
load('/home/ibtuser/Documentos/PLTs/ara/temp/Arabidopsis_thaliana/WGCNA/WGCNA.RData')

# Example reference network data
ref.net <- read.csv("output/Arabidopsis_thaliana/NETWORK/prop_table.csv")
ref.genes <- ref.net$gen_id
ref.genes <- ref.genes[ref.genes %in% gene.names]

# Ensure rownames and colnames of TOM match the gene list
rownames(TOM) <- gene.names
colnames(TOM) <- gene.names

# Subset the TOM to only reference genes
TOM <- TOM[ref.genes, ref.genes]

# Also subset dynamicColors to only reference genes
names(dynamicColors) <- gene.names
dynamicColors <- dynamicColors[ref.genes]

# Compute dissimilarity and do hierarchical clustering with the reference subset
diss <- 1 - TOM
hier <- hclust(as.dist(diss), method="average")
diag(diss) <- NA

# Example adjacency approach with the reference subset
ref.genes.expression <- expression.matrix[,ref.genes]
adjacency <- adjacency(ref.genes.expression, power = 8)
adjacency.diss <- 1 - adjacency
adjacency.hier <- hclust(as.dist(adjacency.diss), method = "average")
dynamicMods <- cutreeDynamic(
  dendro = adjacency.hier, 
  distM = adjacency.diss,
  deepSplit = 2, 
  pamRespectsDendro = FALSE,
  minClusterSize = 30
)
adjacency.modules <- labels2colors(dynamicMods)
diag(adjacency.diss) <- NA

# Visualize adjacency-based clustering in a heatmap
png('output/Arabidopsis_thaliana/WGCNA/dendro.adjacency.refnet.png')
TOMplot(adjacency.diss^4, adjacency.hier, adjacency.modules,
        main = "Adjacency heatmap, genes in PLETHORA reference network")
dev.off()
