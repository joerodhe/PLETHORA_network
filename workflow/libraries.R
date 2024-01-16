update.packages(ask = FALSE)

if (!require("ggplot2", quietly = TRUE))
    install.packages("ggplot2")

if (!require("igraph", quietly = TRUE))
    install.packages("igraph")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!require("edgeR", quietly = TRUE))
    BiocManager::install("edgeR")

if (!require("DESeq2", quietly = TRUE))
    BiocManager::install("DESeq2")

if (!require("optparse", quietly = TRUE))
    install.packages("optparse")

if (!require("svglite", quietly = TRUE))
    install.packages("svglite")

if (!require("factoextra", quietly = TRUE))
    install.packages("factoextra")

if (!require("clusterProfiler", quietly = TRUE))
    BiocManager::install("clusterProfiler")

if (!require("WGCNA", quietly = TRUE))
    BiocManager::install("WGCNA")

if (!require("flashClust", quietly = TRUE))
    install.packages("flashClust")

if (!require("rbiopai", quietly = TRUE))
    install.packages("rbioapi")

if (!require("ape", quietly = TRUE))
    install.packages("ape")

if (!require("VennDiagram", quietly = TRUE))
    install.packages("VennDiagram")

if (!require("e1071", quietly = TRUE))
    install.packages('e1071')

if (!require("mixtools", quietly = TRUE))
    install.packages('mixtools')

if (!require("showtext", quietly = TRUE))
    install.packages('showtext')

if (!require("Seurat", quietly = TRUE))
    install.packages('Seurat')

if (!require("qlcMatrix", quietly = TRUE))
    install.packages('qlcMatrix')

if (!require("SingleCellExperiment", quietly = TRUE))
    BiocManager::install("SingleCellExperiment")

if (!require("ggpattern", quietly = TRUE))
    BiocManager::install("ggpattern")

if (!require("org.At.tair.db", quietly = TRUE))
    BiocManager::install("org.At.tair.db")



install.packages(c("doParallel", "foreach", "ape", "Rdpack", "benchmarkme", "devtools", "GenomicFeatures", "rtracklayer"))


devtools::install_github("drostlab/metablastr")
devtools::install_github("drostlab/rdiamond")
devtools::install_github("drostlab/orthologr")
