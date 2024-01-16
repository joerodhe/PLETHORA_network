
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