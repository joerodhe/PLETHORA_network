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
enableWGCNAThreads(30)


enrich_plot <- function(ora.result, comparison){
  ora.result <- as.data.frame(ora.result)
  plot <- ggplot()
  if(nrow(ora.result) > 0){
  if(nrow(ora.result) > 20) max.n <- 20
  else max.n <- nrow(ora.result)
  ora.result <- head(ora.result[order(ora.result$p.adjust),],max.n)

  plot <- ggplot(ora.result, aes(log2foldEnrichment, fct_reorder(Description, log2foldEnrichment))) +
    geom_segment(aes(xend=0, yend = Description)) +
    geom_point(aes(color=p.adjust, size = Count)) +
    scale_color_gradientn(colours=brewer.pal(3, name = "Dark2"),
        #trans = "log10"
        guide=guide_colorbar(reverse=TRUE, order=1), name = "Adjusted P-value"
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

plot_pca <- function(exp, metadata){

    pca <- prcomp(t(exp))
    pca.exp.var <- expVar(pca)

    plotpc.object <- plotpc(pca, pca.exp.var, metadata)

    return(plotpc.object)
}

# Sets that we work with
expression.matrices <- c(
    "output/Arabidopsis_thaliana/FEATURECOUNTS/TAIR10.gene.counts.tsv",
    "output/Cucumis_sativus/FEATURECOUNTS/Gy14.gene.counts.tsv",
    "output/Glycine_max/FEATURECOUNTS/Wm82.gene.counts.tsv",
    "output/Oryza_sativa_Japonica_Group/FEATURECOUNTS/MSUv7.gene.counts.tsv",
    "output/Solanum_lycopersicum/FEATURECOUNTS/ITAG4.gene.counts.tsv",
    "output/Zea_mays/FEATURECOUNTS/Zm-B73-REFERENCE-NAM-5.0.55.gene.counts.tsv"
)

organisms <- c(
    "Arabidopsis_thaliana",
    "Cucumis_sativus",
    "Glycine_max",
    "Oryza_sativa_Japonica_Group",
    "Solanum_lycopersicum",
    "Zea_mays"
)

organism.ref <- "Arabidopsis_thaliana"
orthologs <- read.delim("output/interspecie/orthologs.bbh.all.tsv", row.names = 1)
nrow(orthologs)

ref.net <- read.csv("others/prop_table.csv")
nrow(ref.net)
ref.net.genes <- ref.net$gen_id[which(ref.net$gen_id %in% orthologs[,organism.ref])]
length(ref.net.genes)


expression <- read_expression(expression.matrices, organisms, orthologs)
dim(expression)
expression[1:6,1:5]

metadata <- data.frame(
    root_zone = sapply(rownames(expression), function(ID) strsplit(ID, '_')[[1]][1]),
    organism = sapply(rownames(expression), function(ID) strsplit(ID, '_')[[1]][2]),
    stringsAsFactors = TRUE
    )

head(metadata)

expression.normalized <- rlog_norm(t(expression), metadata, normalization = "rlog")
dim(expression.normalized)
expression.normalized[1:5,1:5]


plot_pca(expression.normalized, metadata)
gsgave("/home/ibtuser/Documentos/PLTs/ara/others/pca.expression.rlog.orthologs.bbh.all.svg")


expression.corrected <- removeBatchEffect(expression.normalized, batch = metadata$organism)
dim(expression.corrected)
expression.corrected[1:5,1:5]

plot_pca(expression.corrected, metadata)
gsgave("/home/ibtuser/Documentos/PLTs/ara/others/pca.expression.rlog.corrected.organism.orthologs.bbg.all.svg")

nodes <- rownames(expression.corrected)
ref.net.nodes <- rownames(orthologs[orthologs[,organism.ref] %in% ref.net.genes,])
length(ref.net.nodes)
ref.net.nodes <- ref.net.nodes[ref.net.nodes %in% nodes]
length(ref.net.nodes)
nodes.modules <- setNames(
    object = sample(c(0:7), length(nodes), replace = TRUE),
    nodes
)
nodes.modules[ref.net.nodes] <- 8
table(nodes.modules)
samples_list <- split(rownames(metadata), cut(seq_along(rownames(metadata)),18,labels=FALSE))
length(samples_list)

# Object that will contain the expression data
multiExpr <- list()
moduleList <- list()
sampleLabels = c()
idx <- 1
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

lapply(multiExpr, lapply, dim)
lapply(moduleList, length)
lapply(moduleList, table)

multiExpr[[1]]$data[,1:5]
multiExpr[[2]]$data[,1:5]
multiExpr[[3]]$data[,1:5]


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
        # Remove the offending genes and samples
        multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes]
        moduleList[[set]] = moduleList[[set]][gsg$goodGenes]
    }
    # Update exprSize
    exprSize = checkSets(multiExpr)
}

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

Zs <- do.call(
    rbind,
    lapply(
        preservation.results, function(result)
            data.frame(
                Zsummary = result[plt.module,]$Zsummary.pres,
                Zdensity = result[plt.module,]$Zdensity.pres,
                Zconn = result[plt.module,]$Zconnectivity.pres
            )
    ))

Zs$Root_zone <- sapply(rownames(Zs), function(x) strsplit(x, "_")[[1]][1])
Zs$Organism <- sapply(rownames(Zs), function(x) strsplit(x, "_")[[1]][2])

monocots <- c("Oryza", "Zea")

Zs$class <- ifelse(Zs$Organism %in% monocots, "Monocots", "Dicots")

ggplot(Zs, aes(x=Root_zone, y=Zsummary, fill=Organism, pattern = class)) +
    geom_bar_pattern(
        stat="identity", 
        position = position_dodge(preserve = "single"),
        pattern_fill = "black",
        pattern_density = 0.1,
        pattern_spacing = 0.025,
        pattern_angle = 45,
        pattern_key_scale_factor = 0.6,
        color = "black") +
    scale_fill_brewer(palette = "Dark2") +
    scale_pattern_manual(values = c(Monocots = "circle", Dicots = "none")) +
    geom_hline(yintercept = 10, linetype = "dashed", color="red") +
    geom_hline(yintercept = 2, linetype = "dashed", color="orange") +
    labs(title = "Zsummary. Reference: Arabidopsis MZ network", pattern = "") +
    guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none"))) +
    theme_classic()

ggsave("others/zsummary.barplot.pdf")

# KME scatterplots
par(mar = c(3.3, 3.3, 4, 0.5))
par(mgp = c(1.9, 0.6, 0))
ref = 1
ind = 5
for (set in 1:nSets){
if (set!=ref){
verboseScatterplot(KMEpathway[, ref], KMEpathway[, set],
xlab = spaste("KME ", colnames(KMEpathway)[set]),
ylab = "KME MZ_Arabidopsis",  abline = TRUE)
#main = spaste(LETTERS[ind], ". KME in ", organisms[order2[set]], "\nvs. ",
#organisms[ref], "\n"),
#cex.main = 1.4, cex.lab = 1.2, cex.axis = 1.2
ind = ind + 1;
} else {
plot(c(0,1), type ="n", axes = FALSE, xlab = "", ylab ="")
}
}
# If plotting into a file, close it
dev.off()

ref.net.genes.ora <- enrichGO(
    gene = ref.net.genes, 
    keyType = 'TAIR', 
    OrgDb = org.At.tair.db, 
    ont = "All", 
    pAdjustMethod = "BH", 
    pvalueCutoff = 0.05, 
    qvalueCutoff = 0.05
    )

ref.net.genes.ora <- mutate(
    ref.net.genes.ora, 
    log2foldEnrichment = log2((as.numeric(sub("/\\d+", "", GeneRatio)) / as.numeric(sub("\\d+/", "", GeneRatio))) / (as.numeric(sub("/\\d+", "", BgRatio)) / as.numeric(sub("\\d+/", "", BgRatio)))))

ref.net.genes.ora <- as.data.frame(ref.net.genes.ora)
head(ref.net.genes.ora)
ref.net.genes.ora <- ref.net.genes.ora[order(ref.net.genes.ora$log2foldEnrichment, decreasing = TRUE), ]
head(ref.net.genes.ora[,c("Description", "Count"), drop = FALSE], 10)

go.genes <- strsplit(ref.net.genes.ora$geneID[1], "/")[[1]]
go.genes <- lapply(ref.net.genes.ora$geneID[1:10], function(x) strsplit(x, "/")[[1]])
go.genes <- unique(unlist(go.genes))
length(go.genes)


candidates = grep("RNA polymerase", terms2, fixed = TRUE)
length(candidates)
go.terms <- terms2[candidates]

writeLines(paste(names(go.terms), go.terms), "/home/ibtuser/Documentos/PLTs/ara/others/go.terms.txt")

go.ids <- names(terms2[candidates])


go.offspring <- c(
    as.list(GOCCOFFSPRING),
    as.list(GOBPOFFSPRING),
    as.list(GOMFOFFSPRING)
)
length(go.offspring)
head(go.offspring)

go.offspring.ids <- names(go.offspring)

go.relations <- do.call(rbind,
    lapply(go.ids, function(go.id){
        child.terms <- go.offspring[go.offspring.ids == go.id]
        x <- data.frame(
            parent = rep(go.id, length(child.terms)),
            childterm = child.terms[[1]]
        )
        return(x)
        }
    )
)

go.ids.offspring <- unique(c(
    unique(go.relations$parent),
    as.character(na.omit(go.relations$childterm))
))
length(go.ids.offspring)

go.anns <- read.delim("input/Arabidopsis_thaliana/GENOME/TERM2GENE.tsv")
nrow(go.anns)
head(go.anns)

go.genes <- data.frame()
for(i in 1:length(go.ids.offspring)){
    go.genes <- rbind(
        go.genes,
        go.anns[which(go.anns$go.id == go.ids.offspring[i]),]
    )
}
nrow(go.genes)


ref.net.genes.go <- ref.net.genes[ref.net.genes %in% go.genes]
length(ref.net.genes.go)
writeLines(ref.net.genes.go)
ref.net.nodes.go <- rownames(orthologs[orthologs[,organism.ref] %in% ref.net.genes.go,])
length(ref.net.nodes.go)

ref.net.nodes.go <- ref.net.nodes.go[ref.net.nodes.go %in% colnames(multiExpr[[1]]$data)]
length(ref.net.nodes.go)

expr <- list()
for(set in 1:length(multiExpr)){
    expr[[set]] <- list(data = t(apply(multiExpr[[set]]$data, 1, tapply, colnames(multiExpr[[set]]$data), mean, na.rm = TRUE)))
}

build_cor <- function(x){
    x <- bicor(x, nThreads = 16)
    diag(x) <- NA
    x[upper.tri(x)] <- NA
    x.melt <- melt(x, na.rm = TRUE)
    return(x.melt)
}

x <-do.call(cbind, lapply(expr, function(expr.dataset){
        x <- build_cor(expr.dataset$data[,ref.net.nodes])
        rownames(x) <- paste(x[,1], x[,2], sep = "_")
        x <- x[ , 3, drop = FALSE]
        return(x)
    }))

colnames(x) <- names(multiExpr)

col.ref <- c(1,4,7,10,13,16)
col.rest <- setdiff(1:18, col.ref)

edges.var <- apply(x[,col.ref],1,var)

edges.less.var <- sort(edges.var, decreasing = FALSE)[1:100]
nodes.less.var <- unique(unlist(strsplit(names(edges.less.var), "_")))
length(nodes.less.var)

genes.less.var <- orthologs[rownames(orthologs) %in% nodes.less.var,1]

writeLines(genes.less.var)

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
    log2foldEnrichment = log2((as.numeric(sub("/\\d+", "", GeneRatio)) / as.numeric(sub("\\d+/", "", GeneRatio))) / (as.numeric(sub("/\\d+", "", BgRatio)) / as.numeric(sub("\\d+/", "", BgRatio)))))

enrich_plot(genes.less.var.ora, "ok")

ggsave("okkkkk.pdf", width = 7, height = 4)

pathwayAdjs <- list()
pathwayAdjs.pval <- list()
for (set in 1:length(multiExpr)){
    printFlush(paste("Working on set", names(multiExpr)[set]))
    bc <- bicorAndPvalue(expr[[set]]$data[, nodes.less.var], use = "p")
    pathwayAdjs[[set]] <- abs(bc$bicor) * sign(bc$bicor)
    pathwayAdjs.pval[[set]] <- bc$p
}

ref.nodes.degree <- apply(abs(pathwayAdjs[[1]]), 2, sum)-1
ref.nodes.degree <- order(ref.nodes.degree, decreasing = TRUE)


sizeGrWindow(8,24)
par(mfrow =c(6,3))
#par(mfrow =c(4,2))
par(mar = c(2.5, 1, 2.5, 1))
for (set in 1:length(pathwayAdjs)){
    adjs <- pathwayAdjs[[set]]
    adjs[upper.tri(adjs)] <- NA
    diag(adjs) <- NA
    edges <- as.numeric(na.omit(c(adjs)))
    den <- density(edges)
    edges.limit <- max(abs(edges))
    cat("Edge threshold for network", names(multiExpr)[set], "is", edges.limit, "\n")
    plot(den, frame = FALSE, col = "blue",main = names(multiExpr)[set])
    #hist(edges, main = names(multiExpr)[set])
    #abline(v=c(-edges.limit, edges.limit), col = c("blue", "red"), lty=2)
}


#sizeGrWindow(8,24)
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
        variable.cex.labels = TRUE)
}
dev.off()

#multiMEs = multiSetMEs(multiExpr, universalColors = nodes.modules)
#
#mz.nodes.plt.regulated <- mz.nodes.plt.regulated[mz.nodes.plt.regulated %in% names(moduleList[[set]])]
#ref.net.nodes <- ref.net.nodes[ref.net.nodes %in% names(moduleList[[set]])]
#nGenes <- length(ref.net.nodes)
#
#nSets <- length(multiExpr)
#
#KMEpathway = matrix(0, nGenes, nSets)
#for (set in 1:nSets){
#    KMEpathway[, set] <- cor(
#        multiExpr[[set]]$data[,ref.net.nodes], 
#        multiMEs[[set]]$data[, 3], 
#        use = "p")
#}
#
#row.names(KMEpathway) <- ref.net.nodes
#colnames(KMEpathway) <- names(multiExpr)
#

#w <- rcorr(as.matrix(t(x[,col.ref])), type="pearson")
#w.pval <- w$P
#diag(w.pval) <- NA
#w.pval[upper.tri(w.pval)] <- NA
#edges.preservation.pval <- melt(w.pval[1:100,1:100], na.rm = TRUE)
#edges.preserved <- edges.preservation.pval[edges.preservation.pval$value < 0.05,]
#
#
#
#pvals <- c()
#for(i in 1:nrow(x)){
#    w <- wilcox.test(as.numeric(x[i,col.ref]), as.numeric(x[i,col.rest]))
#    pvals <- c(pvals, w$p.value)
#}
#
#names(pvals) <- rownames(x)
#
#pvals.adjusted <- p.adjust(pvals, method = "BH")
#
#edges.diff <- names(pvals[pvals < 0.05])
#nodes.diff <- strsplit(edges.diff, "_")
#length(unique(unlist(nodes.diff)))

## log2(counts + 0.1) All mz expressed and all PLT TFBS (42 genes)
#                Zsummary   Zdensity     Zconn
#DZ_Arabidopsis  8.623043  5.1459132 12.100173
#EZ_Arabidopsis 15.712200 12.5364511 18.887949
#MZ_Cucumis     11.249632 10.0757283 12.423536
#DZ_Cucumis     10.435592  9.3429399 11.528245
#EZ_Cucumis     12.167795 10.4930727 13.842517
#MZ_Glycine      8.203014  6.7627370  9.643292
#DZ_Glycine      3.109191  1.8908642  4.327519
#EZ_Glycine      7.642925  6.5778551  8.707995
#MZ_Oryza        7.444263  6.5707264  8.317800
#DZ_Oryza        4.764301  2.5090091  7.019594
#EZ_Oryza        3.799626  2.4228707  5.176381
#MZ_Solanum     19.433505 15.7646656 23.102344
#DZ_Solanum     10.891124  9.4468852 12.335363
#EZ_Solanum      4.135676  2.4922416  5.779110
#MZ_Zea          2.192255  0.8991507  3.485359
#DZ_Zea          4.190139  3.2928087  5.087470
#EZ_Zea          5.017914  3.7953345  6.240493

# log(counts + 0.1) All mz expressed TFBS in A thaliana (248 genes)
#                Zsummary   Zdensity     Zconn
#DZ_Arabidopsis  5.976255  3.3454466  8.607064
#EZ_Arabidopsis 13.919538 11.3380683 16.501008
#MZ_Cucumis      9.917341  8.3661675 11.468514
#DZ_Cucumis      5.060283  2.6801194  7.440447
#EZ_Cucumis      9.008231  7.7993891 10.217073
#MZ_Glycine      5.908547  3.5108668  8.306227
#DZ_Glycine      4.243834  3.2196507  5.268016
#EZ_Glycine      4.667122  3.5036650  5.830579
#MZ_Oryza        7.844304  6.8117208  8.876888
#DZ_Oryza        4.313809  2.9151574  5.712460
#EZ_Oryza        4.091242  2.6437304  5.538754
#MZ_Solanum     16.712589 14.7390001 18.686178
#DZ_Solanum      9.081776  7.8830167 10.280536
#EZ_Solanum      6.939055  5.9864880  7.891622
#MZ_Zea          2.195394  0.9364931  3.454294
#DZ_Zea          4.235797  2.7686975  5.702897
#EZ_Zea          3.749348  2.2238026  5.274894

# log(counts + 0.1) All mz expressed Manual net as ref net (90 genes)
#                Zsummary   Zdensity     Zconn
#DZ_Arabidopsis  9.410903  6.7242692 12.097537
#EZ_Arabidopsis 18.610396 14.9210954 22.299697
#MZ_Cucumis      9.662665  6.9866600 12.338670
#DZ_Cucumis     10.927949 10.0126763 11.843223
#EZ_Cucumis     14.361193 12.8186509 15.903736
#MZ_Glycine      8.238428  6.8499378  9.626918
#DZ_Glycine      5.200005  4.0844553  6.315555
#EZ_Glycine      4.926837  3.5728254  6.280848
#MZ_Oryza        9.434156  7.5385278 11.329784
#DZ_Oryza        7.348191  6.1404646  8.555918
#EZ_Oryza        1.448590  0.4728998  2.424280
#MZ_Solanum     19.815989 16.7264863 22.905492
#DZ_Solanum     10.176244  9.4420001 10.910488
#EZ_Solanum      5.880304  4.3020399  7.458568
#MZ_Zea          3.144651  1.0663481  5.222954
#DZ_Zea          2.572438  1.5463779  3.598498
#EZ_Zea          6.511062  4.6307016  8.391423

# counts Ara ref net (277 genes)
#                Zsummary  Zdensity     Zconn
#DZ_Arabidopsis 17.742067 12.183208 23.300926
#EZ_Arabidopsis 43.255858 35.884971 50.626746
#MZ_Cucumis     12.927381  9.670542 16.184220
#DZ_Cucumis      9.645396  7.226562 12.064229
#EZ_Cucumis     10.098909  8.253608 11.944210
#MZ_Glycine      6.440306  4.331749  8.548863
#DZ_Glycine     11.772110  8.379467 15.164753
#EZ_Glycine      6.499274  4.156377  8.842171
#MZ_Oryza       11.855119  8.991717 14.718521
#DZ_Oryza        4.818942  2.098337  7.539547
#EZ_Oryza        3.357615  1.146744  5.568487
#MZ_Solanum     21.458302 17.595319 25.321286
#DZ_Solanum      9.935214  7.431170 12.439257
#EZ_Solanum      6.895197  4.446847  9.343546
#MZ_Zea          5.058125  2.101936  8.014314
#DZ_Zea          4.309930  3.435078  5.184782
#EZ_Zea          6.175596  3.913725  8.437467

## counts  All mz expressed and all PLT TFBS (42 genes)
#                Zsummary  Zdensity     Zconn
#DZ_Arabidopsis 19.977804 13.533004 26.422604
#EZ_Arabidopsis 37.889987 31.927226 43.852748
#MZ_Cucumis      6.645665  4.288800  9.002531
#DZ_Cucumis     14.102583 11.588691 16.616475
#EZ_Cucumis      9.382936  7.800779 10.965093
#MZ_Glycine      3.343264  1.361540  5.324988
#DZ_Glycine      8.218757  5.835408 10.602105
#EZ_Glycine      7.557719  4.874728 10.240709
#MZ_Oryza       12.197060  9.126606 15.267513
#DZ_Oryza        6.978280  4.549592  9.406967
#EZ_Oryza        3.861792  2.178007  5.545578
#MZ_Solanum     32.444020 26.694587 38.193453
#DZ_Solanum      8.793021  5.676448 11.909594
#EZ_Solanum     10.316235  7.125908 13.506562
#MZ_Zea          4.535489  1.775825  7.295153
#DZ_Zea          3.727648  1.372713  6.082583
#EZ_Zea          6.127573  3.616684  8.638461