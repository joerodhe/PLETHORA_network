suppressPackageStartupMessages(library(ACTIONet))
suppressPackageStartupMessages(library(NetLibR))
suppressPackageStartupMessages(library(SCINET))
suppressPackageStartupMessages(library(qlcMatrix))
suppressPackageStartupMessages(library(ACTIONetExperiment))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(reshape2))

construct.archetype.signature.profile <- function(sce, ACTIONet.out, reduction_slot = "S_r") {

    # eigengene x archetypes
    reduced.archetype.profile = (t(reducedDims(sce)[[reduction_slot]]) %*% ACTIONet.out$reconstruct.out$C_stacked)

    V <- rowData(sce)$rotation
    rownames(V) <- rownames(rowData(sce))

    perm = order(sapply(colnames(V), function(str) as.numeric(stringr::str_replace(str, "PC", ""))))
    V = V[, perm]

    # gene x archetypes
    archetype.signature.profile = V %*% reduced.archetype.profile
    rownames(archetype.signature.profile) = rownames(sce)

    return(archetype.signature.profile)
}

build_sc_net <- function(ref.net, 
    A, 
    annotation.list, 
    cell.type, 
    genes, 
    marker.genes,
    threads = 8){
    cat("Working on building networks for", cell.type, "cells...\n")
    cat("Subsampling gene activity scores...\n")
    cols <- which(annotation.list == cell.type)
    cat("Number of cells:", length(cols), "\n")
    subsampled.activity.scores <- subsample_gene_activities(
        A = A, 
        rows = 1:nrow(A), 
        samples = cols, 
        thread_no = threads, 
        total_subsamples = 100, 
        cells_per_subsample = 500, 
        seed = 0
    )
    rownames(subsampled.activity.scores) <- genes
    cat("Plotting marker genes activity...\n")
    marker.rows <- match(marker.genes, genes)
    marker.activity.scores <- t(subsampled.activity.scores[marker.rows, ])
    colnames(marker.activity.scores) <- marker.genes
    marker.activity.scores <- reshape2::melt(marker.activity.scores)
    colnames(marker.activity.scores) <- c('id', 'gene', 'activity')

    marker.genes.activity.plot <- ggplot(
        marker.activity.scores, 
        aes(x=activity, fill=gene)
        ) +
        geom_density()
    cat("Building paired datasets object...\n")
    paired.datasets <- pair.datasets(
        ref.net, 
        subsampled.activity.scores
        )
    ref.net.adj <- as(get.adjacency(
        paired.datasets$net), 'TsparseMatrix'
        )
    cat("Building SC networks...\n")
    nets.ref.subsampled.aggregated <- construct_cell_networks_summary(
        net = ref.net.adj, 
        gene_activities = paired.datasets$activity.scores, 
        thread_no = threads, 
        total_subsamples = 100, 
        cells_per_subsample = 500, 
        seed = 0
    )
    cat("Aggregating networks...\n")
    ## Aggregate networks and compute stats
    Mu.ref.subsampled.aggregated <- Reduce(
        "+", nets.ref.subsampled.aggregated) / length(nets.ref.subsampled.aggregated
        )
    net.withRef.subsampled.aggregated.graph <- graph_from_adjacency_matrix(
        Mu.ref.subsampled.aggregated, 
        mode = "undirected", 
        weighted = TRUE
        )
    net.withRef.subsampled.aggregated.graph <- set.vertex.attribute(
        graph = net.withRef.subsampled.aggregated.graph,
        name = "name",
        value = paired.datasets$genes
    )
    return(list(
        net.withRef.subsampled.aggregated.graph, marker.genes.activity.plot
        ))
}

sce <- readRDS('others/sce.counts.RDS')

sce.reduced <- reduce.sce(
  sce, 
  reduced_dim = 50,
  return_V = TRUE
  )

saveRDS(sce.reduced, file = "others/sce.reduced.RDS")

ACTIONet.out <- run.ACTIONet(
  sce.reduced, 
  k_max = 20, 
  thread_no = 100
  )

ACTIONet.out$signature.profile <- construct.archetype.signature.profile(sce.reduced, ACTIONet.out)

saveRDS(ACTIONet.out, file = "others/actinet.out.RDS")

ACTIONet.out <- readRDS("others/actinet.out.RDS")

W <- ACTIONet.out$signature.profile
H <- ACTIONet.out$reconstruct.out$H_stacked
A <-  W %*% H

sce.reduced <- readRDS("others/sce.reduced.RDS")
genes <- rownames(sce.reduced)
cell.types <- as.character(sce.reduced$final.anno)
if("Unknown" %in% cell.types) cell.types <- cell.types[cell.types != "Unknown"]

rm(sce.reduced)
rm(ACTIONet.out)
rm(W)
rm(H)
gc()

marker.genes <- readLines("others/plt_list.txt", warn = FALSE)

ref.net <- read.csv("output/Arabidopsis_thaliana/NETWORK/network_table.csv")
ref.net.graph <- as.matrix(ref.net[,c(1,2)])
ref.net.graph <- graph_from_edgelist(
    ref.net.graph, 
    directed = FALSE
    )

out.path <- "output/Arabidopsis_thaliana/NETWORK"

for(cell.type in unique(cell.types)){

    sc.results <- build_sc_net(
        ref.net = ref.net.graph, 
        A = A, 
        annotation.list = cell.types, 
        cell.type = cell.type, 
        genes = genes,
        marker.genes = marker.genes
        )

    sc.net <- sc.results[[1]]
    saveRDS(sc.net, file.path(out.path, paste0(cell.type, ".sc.igraph.RDS")))

    markers.plot <- sc.results[[2]]
    ggsave(
        markers.plot,
        filename = file.path(out.path, paste0(cell.type, ".marker.genes.svg"))
        )
    
    sc.net.edges <- as.data.frame(cbind(
        get.edgelist(sc.net, names = TRUE), 
        E(sc.net)$weight
        ))
    colnames(sc.net.edges) <- c("source", "target", "weight")
    write.csv(
        sc.net.edges, 
        file.path(out.path, paste0(cell.type, ".sc.csv")), 
        quote = FALSE, 
        row.names = FALSE
        )
    ref.net[[cell.type]] <- apply(ref.net, 1, function(edge){
        sc.net.edges$weight[
            sc.net.edges$source == edge['source'] & sc.net.edges$target == edge["target"]
            ][1]
    })
    rm(sc.net)
    rm(sc.results)
    rm(markers.plot)
    gc()
}

write.csv(
    ref.net, 
    "output/Arabidopsis_thaliana/NETWORK/interaction.table.sc.csv",
    quote = FALSE,
    row.names = FALSE, na = ""
    )

# TO PRINT CELL TYPE ANNOTATIONS

#lapply(colnames(colData(sce.reduced)), function(colname){
#    column <- sce.reduced[[colname]]
#    column.class <- class(column)
#    if(column.class == "character" | column.class == "factor"){
#        tcolumn <- as.data.frame(table(column))
#        colnames(tcolumn) <- c(colname, "Freq")
#        if(nrow(tcolumn) < 20 & nrow(tcolumn) > 4) tcolumn
#    }
#})


#meristem.cols <- which(sce.reduced$time.anno == "Distal Columella")
#activity.scores <- compute_gene_activities(A = A, samples = meristem.cols, thread_no = 8)

#rownames(activity.scores) <- genes
#saveRDS(activity.scores, file = "others/activity.scores.RDS")
#marker.rows <- match(marker.genes, genes)
#
#marker.activity.scores <- t(activity.scores[marker.rows, ])
#colnames(marker.activity.scores) <- marker.genes
#marker.activity.scores <- reshape2::melt(marker.activity.scores)
#colnames(marker.activity.scores) <- c('id', 'gene', 'activity')
#
#ggplot(marker.activity.scores, aes(x=activity, fill=gene)) +
#  geom_density()

#paired.datasets <- pair.datasets(ref.net, activity.scores)
#
#saveRDS(paired.datasets, file = "others/paired.datasets.RDS")
#
#paired.datasets <- readRDS("others/paired.datasets.RDS")
#
#ref.net.adj <- as(get.adjacency(paired.datasets$net), 'TsparseMatrix')
#
## Construct Meristem specific reference-based networks
#meristem.nets.ref <- construct_cell_networks(
#    net = ref.net.adj, 
#    gene_activities = paired.datasets$activity.scores, 
#    thread_no = 4
#    )
#
### Aggregate networks and compute stats
#Mu.ref <- Reduce("+", meristem.nets.ref) / length(meristem.nets.ref)
#meristem.net.withRef.graph <- graph_from_adjacency_matrix(
#    Mu.ref, 
#    mode = "undirected", 
#    weighted = TRUE
#    )
#
#topo.spec.withRef <- topo.spec(
#    meristem.net.withRef.graph, 
#    sample_no = 100
#    )
#print(paired.datasets$genes[order(topo.spec.withRef, decreasing = TRUE)[1:20]])
#
#subsampled.activity.scores <- subsample_gene_activities(
#    A = A, 
#    rows = 1:nrow(A), 
#    samples = meristem.cols, 
#    thread_no = 10, 
#    total_subsamples = 30, 
#    cells_per_subsample = 10, 
#    seed = 0
#    )
#
#rownames(subsampled.activity.scores) <- genes
#
#paired.datasets.subsampled <- pair.datasets(ref.net, subsampled.activity.scores)
#
## Construct meristem specific reference-based networks for the first 30 meristem samples
#meristem.nets.ref.subsampled <- construct_cell_networks(
#    net = ref.net.adj, 
#    gene_activities = paired.datasets.subsampled$activity.scores, 
#    thread_no = 8
#    )
#
#
### Aggregate networks and compute stats
#Mu.ref.subsampled <- Reduce(
#    "+", meristem.nets.ref.subsampled) / length(meristem.nets.ref.subsampled
#    )
#
#meristem.net.withRef.subsampled.graph <- graph_from_adjacency_matrix(
#    Mu.ref.subsampled, 
#    mode = "undirected", 
#    weighted = TRUE
#    )
#
#topo.spec.withRef.subsampled <- topo.spec(
#    meristem.net.withRef.subsampled.graph, 
#    sample_no = 100
#    )
#
#print(paired.datasets$genes[order(topo.spec.withRef.subsampled, decreasing = TRUE)[1:20]])
#
## Construct meristem specific reference-based networks for the first 30 meristem samples
#meristem.nets.ref.subsampled.aggregated <- construct_cell_networks_summary(
#    net = ref.net.adj, 
#    gene_activities = paired.datasets.subsampled$activity.scores, 
#    thread_no = 8, 
#    total_subsamples = 30, 
#    cells_per_subsample = 10, 
#    seed = 0
#    )
#
#
### Aggregate networks and compute stats
#Mu.ref.subsampled.aggregated <- Reduce(
#    "+", meristem.nets.ref.subsampled.aggregated) / length(meristem.nets.ref.subsampled.aggregated
#    )
#
#meristem.net.withRef.subsampled.aggregated.graph <- graph_from_adjacency_matrix(
#    Mu.ref.subsampled.aggregated, 
#    mode = "undirected", 
#    weighted = TRUE
#    )
#
#topo.spec.withRef.subsampled.aggregated <- topo.spec(
#    meristem.net.withRef.subsampled.aggregated.graph, 
#    sample_no = 100
#    )
#
#print(paired.datasets$genes[order(topo.spec.withRef.subsampled.aggregated, decreasing = TRUE)[1:20]])