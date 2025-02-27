suppressPackageStartupMessages(library(optparse))        # For parsing command-line options
suppressPackageStartupMessages(library(ggplot2))         # For plotting
suppressPackageStartupMessages(library(factoextra))      # For visualizing PCA results
suppressPackageStartupMessages(library(DOSE))            # For enrichment analysis themes
suppressPackageStartupMessages(library(forcats))         # For factor reordering in plots
suppressPackageStartupMessages(library(stringr))         # For string manipulation
suppressPackageStartupMessages(library(rjson))           # For handling JSON files
suppressPackageStartupMessages(library(gridExtra))       # For arranging multiple plots together
suppressPackageStartupMessages(library(clusterProfiler)) # For gene set enrichment analysis (GSEA)


# Function to plot PCA results.
# pca_object: PCA output from prcomp()
# exp_var: Dataframe with explained variance information.
plotPCA <- function(pca_object, exp_var){
  # Convert PCA scores to a data frame.
  PC <- as.data.frame(pca_object$x)
  # Assign sample labels (assumes three zones with three samples each).
  PC$Sample <- rep(c('MZ', 'EZ', 'DZ'), each = 3)
  # Create a scatter plot for PC1 vs PC2 colored by sample zone.
  plot <- ggplot(PC, aes(x = PC1, y = PC2, color = Sample)) +
    geom_point() +
    theme_classic() +
    labs(title = "Principal Component Analysis", 
         x = paste0("PC1: ", round(exp_var$exp_var[1]), "% Variance"), 
         y = paste0("PC2: ", round(exp_var$exp_var[2]), "% Variance"), 
         color = "Sample zone")
  return(plot)
}

# Function to calculate the explained variance for each principal component.
# pca_object: PCA output from prcomp()
expVar <- function(pca_object){
  # Variance of each PC is the square of standard deviations.
  exp_var <- pca_object$sdev^2
  # Calculate percentage of total variance explained by each PC.
  exp_var <- (exp_var / sum(exp_var)) * 100 
  # Create a sequence for the number of components.
  nPC <- 1:length(pca_object$sdev)
  # Combine into a data frame.
  exp_var <- data.frame(exp_var, nPC)
  # Calculate the cumulative explained variance.
  exp_var$var_cum <- cumsum(exp_var$exp_var)
  return(exp_var)
}

# Function to plot the explained variance for each principal component.
# exp_var: Dataframe returned by expVar()
# org: Organism name for the plot title.
plotExpVar <- function(exp_var, org){
  # Select a set of colors excluding gray shades.
  color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  color <- sample(color, length(exp_var$nPC), replace = T)
  # Create a bar plot of the explained variance per PC and add cumulative variance line.
  plot <- ggplot(exp_var) +
    geom_bar(aes(x = nPC, y = exp_var), stat = "identity") +
    geom_line(aes(x = nPC, y = var_cum)) +
    geom_point(aes(x = nPC, y = var_cum)) +
    scale_x_continuous(breaks = 1:length(exp_var$nPC)) +
    # Add horizontal lines indicating 50%, 75%, and 90% thresholds.
    geom_hline(yintercept = 90, color = 'Green') +
    annotate(geom = "text", x = length(exp_var$nPC), y = 93, label = "90%", size = 3) +
    geom_hline(yintercept = 75, color = 'Blue') +
    annotate(geom = "text", x = length(exp_var$nPC), y = 78, label = "75%", size = 3) +
    geom_hline(yintercept = 50, color = 'Red') +
    annotate(geom = "text", x = length(exp_var$nPC), y = 53, label = "50%", size = 3) +
    labs(title = paste("Explained variance by each PC for", org, "transcriptome"), 
         x = "Principal Component (PC)", 
         y = "Explained variance (%)") +
    # Add text labels for each bar.
    geom_text(aes(label = paste0(round(exp_var, 0), "%"), x = nPC, y = exp_var), vjust = -0.5, size = 3) +
    theme_classic()
  
  return(plot)
}

# Function to compute the centroid (mean coordinates) for a specific sample type in PCA space.
# sample.type: The sample group label (e.g., "MZ", "EZ", "DZ").
# metadata: Metadata table containing sample annotations.
# pca.components: Matrix of PCA scores.
# components: The specific principal components to consider (e.g., "PC1", "PC2").
get_centroid <- function(sample.type, metadata, pca.components, components){
    # Select sample IDs matching the sample type.
    sampleIDs <- rownames(metadata)[metadata$dex == sample.type]
    # Subset the PCA scores for the selected samples and desired components.
    pca.components <- pca.components[rownames(pca.components) %in% sampleIDs, components]
    # Calculate the mean value for each specified component.
    centroid.coords <- sapply(components, function(component) mean(pca.components[ ,component]))
    return(centroid.coords)
}

# Function to calculate the Euclidean distance between two numeric vectors.
CalculateEuclideanDistance <- function(vect1, vect2) sqrt(sum((vect1 - vect2)^2))

# Function to create a biplot for a specific sample type showing the most contributing genes.
# sample.type: The sample group to highlight.
# gen.list: List of gene names selected by their contribution.
# metadata: Metadata table with sample information.
# pca.exp.var: Dataframe with explained variance of the PCA.
plot_biplot <- function(sample.type, gen.list, metadata, pca.exp.var){
  biplot <- fviz_pca_biplot(
        pca,                                 # PCA object from prcomp()
        geom.ind = "point",                  # Plot individuals as points
        fill.ind = metadata$dex,             # Color individuals by sample group
        pointsize = 5,
        pointshape = 21, 
        repel = TRUE,                        # Avoid overlapping labels
        select.ind = list(name = rownames(metadata)[metadata$dex == sample.type]), 
        select.var = list(name = gen.list),  # Only show selected genes/variables
        col.var = "contrib"                  # Color variables by their contribution
        ) + 
        scale_color_gradient(low = "blue", high = "red") +
        labs(
            title = paste("Variables most contributing to the variance of", sample.type), 
            x = paste0("PC1 ", "(", round(pca.exp.var$exp_var[1], 2), "%)"),
            y = paste0("PC2 ", "(", round(pca.exp.var$exp_var[2], 2), "%)"),
            color = "Contribution", 
            fill = "Zone"
        )  +
        theme(
            legend.position = "bottom", 
            legend.text = element_text(size = 6)
        ) + 
        coord_fixed(ratio = 3/2)
  return(biplot)
}

# Function to perform GO over-representation analysis using Panther DB.
# genes: List of genes of interest.
# ref_genes: List of reference genes.
# zone: The sample zone for which enrichment is performed.
# taxon.id: Taxon identifier.
pca_go <- function(genes, ref_genes, zone, taxon.id){
  # Collapse gene lists into comma-separated strings.
  genes <- paste(genes, collapse = ",")
  ref_genes <- paste(ref_genes, collapse = ",")
  
  # Create a temporary file to store the JSON output from Panther.
  ora.json <- tempfile(fileext = ".json")
  
  # Construct and execute the curl command to query Panther DB for GO enrichment.
  system(paste0('curl -X POST "http://pantherdb.org/services/oai/pantherdb/enrich/overrep?geneInputList=',
   genes, '&organism=', taxon.id, '&refInputList=', ref_genes, '&refOrganism=', taxon.id, '&annotDataSet=GO%3A0008150&enrichmentTestType=FISHER&correction=FDR" -H "accept: application/json" > ', ora.json))
  
  # Read and parse the JSON result.
  ora.result <- fromJSON(file = ora.json)$results$result
  
  # Convert the list of results to a data frame.
  ora.result.df <- as.data.frame(do.call(rbind, lapply(ora.result, function(result) result)))
  ora.result.df$FDR <- as.numeric(ora.result.df$fdr)
  ora.result.df$Count <- as.numeric(ora.result.df$number_in_list)
  # Filter for significant results (FDR <= 0.5).
  ora.result.df <- ora.result.df[which(ora.result.df$FDR <= 0.5), ]
  # Limit the number of results to 20 (or less if fewer results are available).
  if(nrow(ora.result.df) > 20) max_result <- 20
  else max_result <- nrow(ora.result.df)
  ora.result.df <- head(ora.result.df[order(ora.result.df$FDR), ], max_result)
  # Calculate log2 fold enrichment.
  ora.result.df$log2foldEnrichment <- log2(as.numeric(ora.result.df$fold_enrichment))
  # Extract a more descriptive term by splitting the term string.
  ora.result.df$Description <- sapply(ora.result.df$term, function(term) str_split(term, ", ", n = 2, simplify = TRUE)[2])
  
  # Plot the enrichment results.
  plot <- ggplot(ora.result.df, aes(log2foldEnrichment, fct_reorder(Description, log2foldEnrichment))) +
    geom_segment(aes(xend = 0, yend = Description)) +
    geom_point(aes(color = FDR, size = Count)) +
    scale_color_gradientn(colours = c("#f7ca64", "#46bac2", "#7e62a3"),
                          trans = "log10",
                          guide = guide_colorbar(reverse = TRUE, order = 1)) +
    scale_size_continuous(range = c(2, 10)) +
    theme_dose(12) +
    xlab("Log2(Fold-Enrichment)") +
    ylab(NULL) +
    ggtitle(paste("BP GO enriched for", zone)) +
    scale_y_discrete(labels = function(x) str_wrap(x, 45)) +
    theme(axis.text.y = element_text(size = 10))
    
  return(plot)
}

# Function to perform Gene Set Enrichment Analysis (GSEA) on gene distances.
# gene.distance: Named vector of gene distances (used as ranking metric).
# TERM2GENE: Table mapping GO terms to genes.
# TERM2NAME: Table mapping GO term IDs to descriptive names.
pca_gsea <- function(gene.distance, TERM2GENE, TERM2NAME){
  # Sort the gene distances in decreasing order.
  gene.distance <- sort(gene.distance, decreasing = TRUE)
  # Remove duplicated gene entries.
  gene.distance <- gene.distance[!duplicated(gene.distance)]
  
  # Run the GSEA analysis.
  gsea.result <- GSEA(gene.distance, 
                      TERM2GENE = TERM2GENE, 
                      TERM2NAME = TERM2NAME, 
                      seed = TRUE, 
                      scoreType = "pos", 
                      eps = 0, 
                      nPermSimple = 10000)
  gsea.result <- as.data.frame(gsea.result)
  plot <- ggplot()
  if(nrow(gsea.result) > 0){
    # Limit the number of GSEA results to 25.
    if(nrow(gsea.result) > 25) max.n <- 25
    else max.n <- nrow(gsea.result)
    gsea.result <- head(gsea.result[order(gsea.result$p.adjust), ], max.n)
    
    # Plot the GSEA results.
    plot <- ggplot(gsea.result, aes(NES, fct_reorder(Description, NES))) +
      geom_segment(aes(xend = 0, yend = Description)) +
      geom_point(aes(color = p.adjust, size = setSize)) +
      scale_color_gradientn(colours = c("#f7ca64", "#46bac2", "#7e62a3"),
                            trans = "log10",
                            guide = guide_colorbar(reverse = TRUE, order = 1)) +
      scale_size_continuous(range = c(2, 10)) +
      theme_dose(12) +
      xlab("NES") +
      ylab(NULL) +
      ggtitle(paste("BP GO terms enriched")) +
      scale_y_discrete(labels = function(x) str_wrap(x, 45)) +
      theme(axis.text.y = element_text(size = 10))
  }
  return(plot)
}

### Parse Command-Line Options ###
# Define the list of options the script accepts.
option_list = list(
    make_option(c("-m", "--metadata"), type = "character", default = NULL,
              help = "Path to metadata file.", metavar = "character"),
    make_option(c("-n", "--rlog"), type = "character", default = NULL,
              help = "Path to rlog normalized counts matrix.", metavar = "character"),
    make_option(c("-o", "--organism"), type = "character", default = NULL,
              help = "Organism name.", metavar = "character"),
    make_option(c("-t", "--taxon"), type = "character", default = NULL,
              help = "Taxon ID to perform GO enrichment through Panther.db", metavar = "character"),
    make_option(c("-p", "--plots"), type = "character", default = NULL,
              help = "Path to save plots.", metavar = "character"),
    make_option(c("-x", "--tempDir"), type = "character", default = NULL,
              help = "Temporary directory.", metavar = "character"),
    make_option(c("-g", "--term2gene"), type = "character", default = NULL,
              help = "Path to term2gene table file.", metavar = "character"),
    make_option(c("-d", "--term2name"), type = "character", default = NULL,
              help = "Path to term2name table file.", metavar = "character")          
)

# Create an option parser object and parse the provided arguments.
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# Print out the options being used.
cat('Reading metadata file:', opt$metadata, '\n')
cat('Reading rlog matrix:', opt$rlog, '\n')
cat('Working on organism:', opt$organism, '\n')
cat('Taxon ID:', opt$taxon, '\n')
cat('Saving plot in dir:', opt$plots, '\n')
cat('Reading term2gene file:', opt$term2gene, '\n')
cat('Reading term2name file:', opt$term2name, '\n\n')

### Load Input Data ###
# Read gene expression matrix and metadata.
gene.rlog <- read.delim(opt$rlog, check.names = FALSE)
metadata <- read.delim(opt$metadata)
# Extract unique sample types from the metadata.
sample.types <- unique(metadata$dex)
# Read GO term mapping files.
TERM2GENE <- read.delim(opt$term2gene)
TERM2NAME <- read.delim(opt$term2name)
# Filter TERM2GENE to include only GO IDs that are in TERM2NAME.
TERM2GENE <- TERM2GENE[TERM2GENE$go.id %in% TERM2NAME$go.id, ]

cat("Number of genes with GO terms annotated:", length(unique(TERM2GENE$gen.id)), '\n')

### Principal Component Analysis ###
# Perform PCA on the transposed rlog counts (genes as columns, samples as rows).
pca <- prcomp(t(gene.rlog))
# Save the PCA object for future use.
saveRDS(pca, file = file.path(opt$tempDir, "pca.RDS"))
# Calculate explained variance for each PC.
pca.exp.var <- expVar(pca)

# Generate PCA plot and explained variance plot.
pca.plot <- plotPCA(pca, pca.exp.var)
pca.exp.var.cum <- plotExpVar(pca.exp.var, opt$organism)

# Save the plots as SVG files.
ggsave(plot = pca.plot, file = file.path(opt$plots, "pca.svg"))
ggsave(pca.exp.var.cum, file = file.path(opt$plots, "var.svg"))

# Extract details about variable contributions for PC1 and PC2.
axes <- c(1,2)
pca.var.details <- facto_summarize(pca, element = "var", result = c("coord", "contrib", "cos2"), axes = axes)

# Get PCA coordinates for individuals.
pca.ind.details <- get_pca_ind(pca)
pca.ind.details <- data.frame(pca.ind.details$coord[, axes, drop = FALSE], stringsAsFactors = TRUE)

# Calculate a scaling ratio to adjust variable coordinates for the biplot.
r <- min(
   (max(pca.ind.details[,"Dim.1"]) - min(pca.ind.details[,"Dim.1"]) / (max(pca.var.details[,"Dim.1"]) - min(pca.var.details[,"Dim.1"]))),
   (max(pca.ind.details[,"Dim.2"]) - min(pca.ind.details[,"Dim.2"]) / (max(pca.var.details[,"Dim.2"]) - min(pca.var.details[,"Dim.2"])))
  )
pca.var.details[, c("Dim.1", "Dim.2")] <- pca.var.details[, c("Dim.1", "Dim.2")] * (r * 0.7)

# Compute centroids (mean positions) for each sample type in PCA space.
centroids <- lapply(sample.types, get_centroid, metadata = metadata, pca.components = pca$x, components = c("PC1", "PC2"))
names(centroids) <- sample.types

# Identify top contributing variables for each sample type.
var.zone <- c()
for(sample.type in sample.types){
    centroid.coords <- centroids[sample.type]
    pca.var.details[[sample.type]] <- apply(pca.var.details, 1, function(row){
        1 / CalculateEuclideanDistance(as.numeric(row[2:3]), centroid.coords[[1]])
    })
    var.zone <- c(var.zone, 
      head(rownames(pca.var.details[order(pca.var.details[[sample.type]], decreasing = TRUE), ]), 10)
      )
}

# Create an overall biplot using the selected top contributing variables.
biplot <- fviz_pca_biplot(
    pca,
    geom.ind = "point",
    fill.ind = metadata$dex, 
    pointsize = "cos2",
    pointshape = 21, 
    repel = TRUE, 
    select.var = list(name = var.zone),
    col.var = "contrib"
    ) + 
    scale_color_gradient(low = "blue", high = "red") +
    labs(
        title = paste("Variables most contributing to the variance of", opt$organism, "samples"), 
        x = paste0("PC1 ", "(", round(pca.exp.var$exp_var[1], 2), "%)"),
        y = paste0("PC2 ", "(", round(pca.exp.var$exp_var[2], 2), "%)"),
        color = "Contribution", fill = "Zone"
    ) 

# Save the overall biplot.
ggsave(plot = biplot, file = file.path(opt$plots, "biplot.svg"))

### Detailed Biplots and Enrichment Analysis for Each Sample Type ###
for(sample.type in sample.types){
  # Select top 10 genes with highest contribution for the sample type.
  sample.genes <- rownames(head(pca.var.details[order(pca.var.details[[sample.type]], decreasing = TRUE), ], 10))
  
  # Generate a biplot for the specific sample type.
  biplot <- plot_biplot(
    sample.type = sample.type, 
    gen.list = sample.genes, 
    metadata = metadata,
    pca.exp.var = pca.exp.var
  )
  
  # Perform GSEA for the selected genes using their contribution scores.
  gsea.plot <- pca_gsea(
    setNames(object = pca.var.details[[sample.type]], rownames(pca.var.details)),
    TERM2GENE = TERM2GENE,
    TERM2NAME = TERM2NAME
  )
  
  # Combine the biplot and GSEA plot into a single plot grid.
  biplot.enriched <- grid.arrange(biplot, gsea.plot,
    nrow = 1, top = paste("Principal Component Analysis for", opt$organism))
  
  # Save the combined plot for each sample type.
  ggsave(plot = biplot.enriched, width = 14, height = 6, units = "in",
      file = file.path(opt$plots, paste(sample.type, "biplot", "svg", sep = ".")))
}
