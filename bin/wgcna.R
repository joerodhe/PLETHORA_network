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
options(stringsAsFactors = FALSE)
showtext_auto()
set.seed(3006)

checkExpression <- function(expression.matrix){
  gsg = goodSamplesGenes(expression.matrix, verbose = 3)
  if (!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
      dynamicTreeCut::printFlush(
        paste("Removing genes:", paste(names(expression.matrix)[!gsg$goodGenes], collapse = ", "))
        );
    if (sum(!gsg$goodSamples)>0)
      dynamicTreeCut::printFlush(
        paste("Removing samples:", paste(rownames(expression.matrix)[!gsg$goodSamples], collapse = ", "))
        );
    # Remove the offending genes and samples from the data:
    expression.matrix = expression.matrix[gsg$goodSamples, gsg$goodGenes]
  }
  return(expression.matrix)
}

module.enrichment <- function(module, TERM2GENE, TERM2NAME, background, output.path){

  cat("Performing Overrepresentantion functional enrichment over", module, "module...\n")

  plot <- ggplot()
  
  if(file.exists(file.path(output.path, paste0(module, '.nodes.tsv')))){

  module.nodes <- read.delim(file.path(output.path, paste0(module, '.nodes.tsv')))
  query <- module.nodes[,1]

  #gene.list.file <- file.path(output.path, paste0(module, ".genes.txt"))
  #cat("Writing gene module list (", length(query),") in file:", gene.list.file, "\n")
  #writeLines(query, con = gene.list.file)

  query <- query[which(query %in% background)]

  cat("Number of genes in background:", length(background), "\n")
  cat("Number of genes in query (after matching with background):", length(query), "\n")

  
  if (length(query) > 0) {

    ora.results <- enricher(
      gene = query,
      universe = background,
      TERM2GENE = TERM2GENE,
      TERM2NAME = TERM2NAME
    )

    cat("Number of terms enriched:", nrow(ora.results), "\n")

    if (!is.null(ora.results)) if (nrow(ora.results) > 0) {

      ora.results <- mutate(ora.results, 
                      log2foldEnrichment = log2((as.numeric(sub("/\\d+", "", GeneRatio)) / as.numeric(sub("\\d+/", "", GeneRatio))) / (as.numeric(sub("/\\d+", "", BgRatio)) / as.numeric(sub("\\d+/", "", BgRatio)))))

      plot <- ggplot(ora.results, showCategory = 10, aes(log2foldEnrichment, fct_reorder(Description, log2foldEnrichment))) +
      geom_segment(aes(xend=0, yend = Description)) +
      geom_point(aes(color=p.adjust, size = Count)) +
      scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
      trans = "log10",
      guide=guide_colorbar(reverse=TRUE, order=1)) +
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

decompose_gmm <- function(edge.weight){
  # identify outliers
  stats <- boxplot.stats(edge.weight)
  outliers <- which(edge.weight %in% stats$out)
  noutliers <- length(outliers)
  cat('Number of outliers:', noutliers, '\n')
  if( noutliers == 0) outliers <- length(edge.weight)+1
  # decompose gaussian distribution
  maxit = 10000
  cat('Performing GMM with', maxit, 'max iterations...\n')
  gaussian.dist <- normalmixEM(edge.weight[-outliers], k=2, maxit=maxit)
  # extract probabilities for each observation to belong to one or another group
  # of the decomposed gaussian distribution
  cat('Done...\n')
  probs <- as.data.frame(
    gaussian.dist$posterior,
    row.names = names(edge.weight)[-outliers]
    )
  if (noutliers > 0) {
    # add probabilities of outliers (0 or 1 probs)
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
  # add raw (weight) data
  cat('Building dataframe...\n')
  decomposed.dist <- cbind(edge.weight, probs)
  

  return(list(
    mu=gaussian.dist$mu, 
    lambda=gaussian.dist$lambda,
    sigma=gaussian.dist$sigma, 
    probabilities=decomposed.dist))
}

plot_group_proportion <- function(dist.dat, observation.type, plot.title){
  plot <- ggplot(dist.dat, aes(x=factor(group))) +
    geom_bar() +
    xlab(element_blank()) +
    ylab(paste("Number of", observation.type)) +
    labs(title=plot.title) +
    theme_classic() #+
    #theme(axis.title=element_text(size=18), title = element_text(size=22))

  return(plot)
}

plot_dist <- function(dist.dat, plot.title, connected.group){

  connected.group <- ifelse(connected.group == 'comp.1', 1, 2)
  unconnected.group <- ifelse(connected.group == 'comp.1', 2, 1)

  plot <- data.frame(x = dist.dat$probabilities$edge.weight) %>%
    ggplot() +
    geom_histogram(aes(x, after_stat(density)), binwidth = 0.05, colour = "black", 
                  fill = "white") +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(
                    dist.dat$mu[connected.group], 
                    dist.dat$sigma[connected.group], 
                    lam = dist.dat$lambda[connected.group]),
                  colour = "red", lwd = 1) +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(
                    dist.dat$mu[unconnected.group], 
                    dist.dat$sigma[unconnected.group], 
                    lam = dist.dat$lambda[unconnected.group]),
                  colour = "blue", lwd = 1) +
    ylab("Density") + 
    xlab("Topological Overlap Value") +
    labs(title=plot.title) +
    theme_classic() #+
    #theme(axis.title=element_text(size=18), title = element_text(size=22))

  return(plot)
}

plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

module.export <- function(module,
                          dynamicColors, 
                          gene.names, 
                          TOM, 
                          temp.dir, 
                          out.dir, 
                          org.name, 
                          datKME,
                          meristem.nodes){
  # TOM histogram
  cat('Decomposing', module, 'edge weights...\n')

  inModule <- is.finite(match(dynamicColors,module))
  modGenes <- gene.names[inModule]
  modTOM <- TOM[inModule,inModule]
  dimnames(modTOM) <- list(modGenes,modGenes)
  cyt <- exportNetworkToCytoscape(
    modTOM, 
    edgeFile = file.path(temp.dir, paste0(module, '.edges.tsv')), 
    nodeFile = file.path(temp.dir, paste0(module, '.nodes.tsv')), 
    weighted = TRUE, 
    threshold = 0, 
    nodeNames = modGenes
    )

  mod.edges <- read.delim(file.path(temp.dir, paste0(module, '.edges.tsv')))
  mod.weights <- setNames(
    object = mod.edges$weight,
    rownames(mod.edges)
  )

  TOM.GMM <- decompose_gmm(mod.weights)

  # add group based on probability threshold
  cat('Spliting groups based on probability threshold...\n')
  weight.comp.1 <- TOM.GMM$probabilities$edge.weight[TOM.GMM$probabilities$comp.1 >= 0.5]
  
  weight.comp.1.median <- median(weight.comp.1)
  cat('Median edge weight of distribution 1 in', module, 'module is', weight.comp.1.median, '\n')
  weight.comp.2 <- TOM.GMM$probabilities$edge.weight[TOM.GMM$probabilities$comp.2 >= 0.5]
  weight.comp.2.median <- median(weight.comp.2)
  cat('Median edge weight of distribution 2 in', module, 'module is', weight.comp.2.median, '\n')
  connected.group <- ifelse(
    weight.comp.1.median > weight.comp.2.median,
    'comp.1',
    'comp.2'
  )
  if(is.na(connected.group)) connected.group <- 'comp.2'
  cat('Connected group for', module, 'is', connected.group, '\n')
  TOM.GMM$probabilities$group <- ifelse(
    TOM.GMM$probabilities[[connected.group]] >= 0.5,
    "Connected",
    "NotConnected"
  )

  cat('Number of connected nodes in', 
    module, 'module is', nrow(TOM.GMM$probabilities[which(TOM.GMM$probabilities$group == 'Connected'), ]),
    '\n')

  cat('Writing GMM results...\n')

  write.table(
    as.data.frame(TOM.GMM$probabilities), 
    file = file.path(temp.dir, paste0(module, '.gmm.tsv')),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE)

  edge.threshold <- min(TOM.GMM$probabilities$edge.weight[which(TOM.GMM$probabilities$group == 'Connected')])

  cat('Plotting...\n')
  gmm.propotion.plot <- plot_group_proportion(TOM.GMM$probabilities,
    observation.type="edges",
    plot.title=paste("Proportion of edges\nMin W:", round(edge.threshold, 2))
    )

  gmm.distribution.plot <- plot_dist(TOM.GMM, 
    plot.title='Topological Overlap Matrix gaussian decomposition',
    connected.group=connected.group
    )

  gmm.plot <- grid.arrange(gmm.propotion.plot, gmm.distribution.plot, widths=c(0.3, 0.7),
    nrow = 1, top = paste("Gaussian Model Mixture of", org.name, "WGNCA", module, "module edge weights"))
  
  x <- ggsave(
    plot = gmm.plot, 
    width = 10, height = 4, 
    file = file.path(out.dir, paste0(module, '.gmm.svg')))

  datKME <- datKME[modGenes, module]

  hub.threshold <- quantile(datKME, probs = 0.9)


  cat('Hub threshold for', module, 'is', hub.threshold, '(corresponding to 10% of higher kME)\n')

  nodes.metadata <- data.frame(
    kME = datKME,
    hub = ifelse(datKME > hub.threshold, TRUE, FALSE),
    MZ.expressed = ifelse(modGenes %in% meristem.nodes, TRUE, FALSE),
    row.names = modGenes
  )

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
              help = "opt$mergeThreshold threshold to merge WGCNA modules (eigengenes).", metavar = "double"),
    make_option(c("-f", "--forceTOMCalc"), type = "character", default = NULL,
              help = "Force re-calculation of TOM matrix.", metavar = "character"),   
    make_option(c("-u", "--temp"), type = "character", default = NULL,
              help = "Temp path", metavar = "character")       
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

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

enableWGCNAThreads(opt$threads)

temp.WGCNA <- file.path(opt$temp, 'WGCNA.RData')

if(file.exists(temp.WGCNA) & opt$forceTOMCalc == "False"){
  cat("Loading data of previous session...\n")
  load(temp.WGCNA)
} else {
  expression.matrix <- read.delim(opt$expression)
  expression.matrix <- t(expression.matrix)
  # check for genes and samples with too many missing values
  expression.matrix <- checkExpression(expression.matrix)
  gene.names <- colnames(expression.matrix)
  # check sample outliers
  sampleTree = hclust(dist(expression.matrix), method = "average")
  #Choosing a soft-threshold to fit a scale-free topology to the network
  cat("Analyzing powers...\n")
  powers = c(c(1:20), seq(from = 22, to=30, by=2))
  sft = pickSoftThreshold(expression.matrix, powerVector = powers)
  cat("Calculating TOM similarity...\n")
  adjacency = adjacency(expression.matrix, power = opt$softPower)
  TOM = TOMsimilarity(adjacency, verbose = 3)
  rm(adjacency)
  dissTOM = 1-TOM
  cat("Running hierarchical clustering over TOM dissimilarity matrix...\n")
  # Call the hierarchical clustering function
  geneTree = hclust(as.dist(dissTOM), method = "average")
  # We like large modules, so we set the minimum module size relatively high:
  minModuleSize = 30
  cat("Performing module identification through dynamic tree cut...\n")
  # Module identification using dynamic tree cut
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                  deepSplit = 2, pamRespectsDendro = FALSE,
                  minClusterSize = minModuleSize)
  rm(dissTOM)
  # Convert numeric lables into colors
  dynamicColors = labels2colors(dynamicMods)
  cat("Calculating eigengenes...\n")
  # Calculate eigengenes
  MEList = moduleEigengenes(expression.matrix, colors = dynamicColors)
  MEs = MEList$eigengenes
  rm(MEList)
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs)
  rm(MEs)
  cat("Performing hierarchical clusting over eigengenes...\n")
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average")
  cat("Merging modules...\n")
  # Call an automatic merging function
  merge = mergeCloseModules(expression.matrix, dynamicColors, cutHeight = opt$mergeThreshold, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs
  rm(merge)
  # save session
  gc()
  cat("Saving session...\n")
  cat('Writing file', temp.WGCNA, '\n')
  save.image(file = temp.WGCNA)
}

MEList = moduleEigengenes(expression.matrix, colors = dynamicColors)
datME = MEList$eigengenes

datKME <- signedKME(expression.matrix, datME, outputColumnName='')
datKME.file <- file.path(opt$temp, 'kme.tsv')

cat('Saving KME table file in:', datKME.file, '\n')
write.table(datKME, file = datKME.file, sep = '\t', quote = FALSE)

background <- readLines(opt$background)

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

# plot sample tree
svg(file.path(opt$output, "sampleTree.svg"))
par(cex = 1.3);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.0,cex.axis = 1.5, cex.main = 1)
dev.off()

# Plot powers 
svg(file.path(opt$output, "powers.svg"))
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence. Selected threshold:", opt$softPower));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    labels=powers,cex=cex1,col="red");
#Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
    xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
    main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# Plot the eigengenes hierarchical clustering
svg(file.path(opt$output, "eigengenes.svg"))
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
# Plot the cut line into the dendrogram
abline(h=opt$mergeThreshold, col = "red")
dev.off()

# Plot dendogram
svg(file.path(opt$output, "dendro.svg"))
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

TERM2GENE <- read.delim(opt$term2gene)
TERM2NAME <- read.delim(opt$term2name)
TERM2GENE <- TERM2GENE[TERM2GENE$go.id %in% TERM2NAME$go.id, ]

background <- background[which(background %in% TERM2GENE$gen.id)]

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

load('/home/ibtuser/Documentos/PLTs/ara/temp/Arabidopsis_thaliana/WGCNA/WGCNA.RData')
ref.net <- read.csv("output/Arabidopsis_thaliana/NETWORK/prop_table.csv")
ref.genes <- ref.net$gen_id
ref.genes <- ref.genes[ref.genes %in% gene.names]
rownames(TOM) <- gene.names
colnames(TOM) <- gene.names
TOM <- TOM[ref.genes, ref.genes]
names(dynamicColors) <- gene.names
dynamicColors <- dynamicColors[ref.genes]
#color1=dynamicColors
#restGenes= (color1 != "grey")
diss <- 1-TOM
hier <- hclust(as.dist(diss), method="average")
diag(diss) <- NA

ref.genes.expression <- expression.matrix[,ref.genes]
adjacency <- adjacency(ref.genes.expression, power = 8)
adjacency.diss <- 1-adjacency
adjacency.hier <- hclust(as.dist(adjacency.diss), method = "average")
dynamicMods <- cutreeDynamic(dendro = adjacency.hier, distM = adjacency.diss,
                  deepSplit = 2, pamRespectsDendro = FALSE,
                  minClusterSize = 30)
adjacency.modules <- labels2colors(dynamicMods)
diag(adjacency.diss) <- NA


png('output/Arabidopsis_thaliana/WGCNA/dendro.adjacency.refnet.png')
#sizeGrWindow(7,7)
TOMplot(adjacency.diss^4, adjacency.hier, adjacency.modules, main = "Adjacency heatmap, genes in PLETHORA reference network" )
dev.off()
