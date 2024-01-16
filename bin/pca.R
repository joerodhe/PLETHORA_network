suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(clusterProfiler))

plotPCA <- function(pca_object, exp_var){
  PC <- as.data.frame(pca_object$x)
  PC$Sample <- rep(c('MZ', 'EZ', 'DZ'), each = 3)
  plot <- ggplot(PC, aes(x = PC1, y = PC2, color = Sample)) +
    geom_point() +
    theme_classic() +
    labs(title = "Principal Component Analysis", 
         x = paste0("PC1: ", round(exp_var$exp_var[1]), "% Variance"), 
         y = paste0("PC2: ", round(exp_var$exp_var[2]), "% Variance"), color = "Sample zone")
  return(plot)
}

expVar <- function(pca_object){
  exp_var <- pca_object$sdev^2
  exp_var <- (exp_var / sum(exp_var)) * 100 # Proporcion de la varianza explicada
  nPC <- 1:length(pca_object$sdev) # Numero de componentes
  exp_var <- data.frame(exp_var, nPC)
  exp_var$var_cum <- cumsum(exp_var$exp_var)
  return(exp_var)
}

plotExpVar <- function(exp_var, org){
  color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  color <- sample(color, length(exp_var$nPC), replace = T)
  plot <- ggplot(exp_var) +
    geom_bar(aes(x = nPC, y = exp_var), stat = "identity") +
    geom_line(aes(x = nPC, y = var_cum)) +
    geom_point(aes(x = nPC, y = var_cum)) +
    scale_x_continuous(breaks = 1:length(exp_var$nPC)) +
    geom_hline(yintercept = 90, color = 'Green') +
    annotate(geom = "text", x=length(exp_var$nPC), y=93, label="90%", size=3) +
    geom_hline(yintercept = 75, color = 'Blue') +
    annotate(geom = "text", x=length(exp_var$nPC), y=78, label="75%", size=3) +
    geom_hline(yintercept = 50, color = 'Red') +
    annotate(geom = "text", x=length(exp_var$nPC), y=53, label="50%", size=3) +
    labs(title = paste("Explained variance by each PC for", org,"transcriptome"), x = "Principal Component (PC)", y = "Explained variance (%)") +
    geom_text(aes(label=paste0(round(exp_var, 0), "%"), x=nPC, y=exp_var), vjust=-0.5, size=3) +
    theme_classic()
  
  return(plot)
  }

get_centroid <- function(sample.type, metadata, pca.components, components){
    sampleIDs <- rownames(metadata)[metadata$dex == sample.type]
    pca.components <- pca.components[rownames(pca.components) %in% sampleIDs, components]
    centroid.coords <- sapply(components, function(component) mean(pca.components[ ,component]))
    #centroid.angle <- angle(x_col=centroid.coords[1], y_col=centroid.coords[2])
    return(centroid.coords)
}

CalculateEuclideanDistance <- function(vect1, vect2) sqrt(sum((vect1 - vect2)^2))

plot_biplot <- function(sample.type, gen.list, metadata, pca.exp.var){
  biplot <- fviz_pca_biplot(
        pca,
        geom.ind = "point",
        fill.ind = metadata$dex, 
        pointsize = 5,
        pointshape = 21, 
        repel = TRUE,
        select.ind = list(name = rownames(metadata)[metadata$dex == sample.type]), 
        select.var = list(name = gen.list),
        col.var = "contrib"
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
            legend.position="bottom", 
            legend.text = element_text( size=6)
        ) + 
        coord_fixed(ratio = 3/2)
  return(biplot)
}

pca_go <- function(genes, ref_genes, zone, taxon.id){

  genes <- paste(genes, collapse = ",")
  ref_genes <- paste(ref_genes, collapse = ",")

  ora.json <- tempfile(fileext = ".json")

  system(paste0('curl -X POST "http://pantherdb.org/services/oai/pantherdb/enrich/overrep?geneInputList=',
   genes,'&organism=',taxon.id, '&refInputList=',ref_genes,'&refOrganism=',taxon.id,'&annotDataSet=GO%3A0008150&enrichmentTestType=FISHER&correction=FDR" -H "accept: application/json" > ', ora.json))

  ora.result <- fromJSON(file = ora.json)$results$result

  ora.result.df <- as.data.frame(do.call(rbind, lapply(ora.result, function(result) result)))
  ora.result.df$FDR <- as.numeric(ora.result.df$fdr)
  ora.result.df$Count <- as.numeric(ora.result.df$number_in_list)
  ora.result.df <- ora.result.df[which(ora.result.df$FDR <= 0.5), ]
  if(nrow(ora.result.df) > 20) max_result <- 20
  else max_result <- nrow(ora.result.df)
  ora.result.df <- head(ora.result.df[order(ora.result.df$FDR), ], max_result)
  ora.result.df$log2foldEnrichment <- log2(as.numeric(ora.result.df$fold_enrichment))
  ora.result.df$Description <- sapply(ora.result.df$term, function(term) str_split(term, ", ", n=2, simplify = TRUE)[2])

  plot <- ggplot(ora.result.df, aes(log2foldEnrichment, fct_reorder(Description, log2foldEnrichment))) +
    geom_segment(aes(xend=0, yend = Description)) +
    geom_point(aes(color=FDR, size = Count)) +
    scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
    trans = "log10",
    guide=guide_colorbar(reverse=TRUE, order=1)) +
    scale_size_continuous(range=c(2, 10)) +
    theme_dose(12) +
    xlab("Log2(Fold-Enrichment)") +
    ylab(NULL) +
    ggtitle(paste("BP GO enriched for", zone)) +
        scale_y_discrete(labels = function(x) str_wrap(x, 45))
    theme(axis.text.y = element_text(size = 10))
    
    return(plot)
}

pca_gsea <- function(gene.distance, TERM2GENE, TERM2NAME){
  gene.distance <- sort(gene.distance, decreasing = TRUE)
  gene.distance <- gene.distance[!duplicated(gene.distance)]

  gsea.result <- GSEA(gene.distance, 
                      TERM2GENE = TERM2GENE, 
                      TERM2NAME = TERM2NAME, 
                      seed = TRUE, 
                      scoreType = "pos", 
                      eps=0, nPermSimple=10000)
  gsea.result <- as.data.frame(gsea.result)
  plot <- ggplot()
  if(nrow(gsea.result) > 0){
  if(nrow(gsea.result) > 25) max.n <- 25
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
    ggtitle(paste("BP GO terms enriched")) +
        scale_y_discrete(labels = function(x) str_wrap(x, 45)) +
    theme(axis.text.y = element_text(size = 10))
  }
  return(plot)
}

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
              help = "temporary directory.", metavar = "character"),
    make_option(c("-g", "--term2gene"), type = "character", default = NULL,
              help = "Path to term2gene table file.", metavar = "character"),
    make_option(c("-d", "--term2name"), type = "character", default = NULL,
              help = "Path to term2name table file.", metavar = "character")          
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

cat('Reading metadata file:', opt$metadata, '\n')
cat('Reading rlog matrix:', opt$rlog, '\n')
cat('Working on organism:', opt$organism, '\n')
cat('Taxon ID:', opt$taxon, '\n')
cat('Saving plot in dir:', opt$plots, '\n')
cat('Reading term2gene file:', opt$term2gene, '\n')
cat('Reading term2name file:', opt$term2name, '\n\n')


gene.rlog <- read.delim(opt$rlog, check.names = FALSE)
metadata <- read.delim(opt$metadata)
sample.types <- unique(metadata$dex)
TERM2GENE <- read.delim(opt$term2gene)
TERM2NAME <- read.delim(opt$term2name)

TERM2GENE <- TERM2GENE[TERM2GENE$go.id %in% TERM2NAME$go.id, ]

cat("Number of genes with GO terms annotated:", length(unique(TERM2GENE$gen.id)), '\n')

pca <- prcomp(t(gene.rlog))
saveRDS(pca, file = file.path(opt$tempDir, "pca.RDS"))
pca.exp.var <- expVar(pca)

pca.plot <- plotPCA(pca, pca.exp.var)
pca.exp.var.cum <- plotExpVar(pca.exp.var, opt$organism)

ggsave(plot = pca.plot, file = file.path(opt$plots, "pca.svg"))
ggsave(pca.exp.var.cum, file = file.path(opt$plots, "var.svg"))

axes <- c(1,2)
pca.var.details <- facto_summarize(pca, element = "var", result = c("coord", "contrib", "cos2"), axes = axes)

pca.ind.details <- get_pca_ind(pca)
pca.ind.details <- data.frame(pca.ind.details$coord[, axes, drop=FALSE], stringsAsFactors = TRUE)

r <- min(
   (max(pca.ind.details[,"Dim.1"])-min(pca.ind.details[,"Dim.1"])/(max(pca.var.details[,"Dim.1"])-min(pca.var.details[,"Dim.1"]))),
   (max(pca.ind.details[,"Dim.2"])-min(pca.ind.details[,"Dim.2"])/(max(pca.var.details[,"Dim.2"])-min(pca.var.details[,"Dim.2"])))
  )

pca.var.details[,c("Dim.1", "Dim.2")] <- pca.var.details[,c("Dim.1", "Dim.2")] * (r*0.7)


centroids <- lapply(sample.types, get_centroid, metadata = metadata, pca.components = pca$x, components = c("PC1", "PC2"))
names(centroids) <- sample.types

var.zone <- c()
for(sample.type in sample.types){
    centroid.coords <- centroids[sample.type]
    pca.var.details[[sample.type]] <- apply(pca.var.details, 1, function(row){
        1/CalculateEuclideanDistance(as.numeric(row[2:3]), centroid.coords[[1]])
    })
    var.zone <- c(var.zone, 
      head(rownames(pca.var.details[order(pca.var.details[[sample.type]], decreasing = TRUE), ]), 10)
      )
}

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

ggsave(plot = biplot, file = file.path(opt$plots, "biplot.svg"))


for(sample.type in sample.types){

  sample.genes <- rownames(head(pca.var.details[order(pca.var.details[[sample.type]], decreasing=TRUE),], 10))

  biplot <- plot_biplot(
    sample.type = sample.type, 
    gen.list = sample.genes, 
    metadata = metadata,
    pca.exp.var = pca.exp.var
  )

  gsea.plot <- pca_gsea(
    setNames(object = pca.var.details[[sample.type]], rownames(pca.var.details)),
    TERM2GENE = TERM2GENE,
    TERM2NAME = TERM2NAME
  )

  biplot.enriched <- grid.arrange(biplot, gsea.plot,
    nrow = 1, top = paste("Principal Component Analysis for", opt$organism))

    ggsave(plot = biplot.enriched, width = 14, height = 6, units = "in",
      file = file.path(opt$plots, paste(sample.type, "biplot", "svg", sep = ".")))
}
