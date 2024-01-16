suppressPackageStartupMessages(library(ACTIONet))
suppressPackageStartupMessages(library(NetLibR))
suppressPackageStartupMessages(library(SCINET))
suppressPackageStartupMessages(library(qlcMatrix))
suppressPackageStartupMessages(library(ACTIONetExperiment))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(mixtools))
suppressPackageStartupMessages(library(RColorBrewer))

build_sc_network <- function(ref.net, 
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
    net.withRef.subsampled.aggregated.graph <- NA
    marker.genes.activity.plot <- NA
    #marker.genes.activity.plot <- ggplot(
    #    marker.activity.scores, 
    #    aes(x=activity, fill=gene)
    #    ) +
    #    geom_density()
    #cat("Building paired datasets object...\n")
    #paired.datasets <- pair.datasets(
    #    ref.net, 
    #    subsampled.activity.scores
    #    )
    #ref.net.adj <- as(get.adjacency(
    #    paired.datasets$net), 'TsparseMatrix'
    #    )
    #cat("Building SC networks...\n")
    #nets.ref.subsampled.aggregated <- construct_cell_networks_summary(
    #    net = ref.net.adj, 
    #    gene_activities = paired.datasets$activity.scores, 
    #    thread_no = threads, 
    #    total_subsamples = 100, 
    #    cells_per_subsample = 500, 
    #    seed = 0
    #)
    #cat("Aggregating networks...\n")
    ### Aggregate networks and compute stats
    #Mu.ref.subsampled.aggregated <- Reduce(
    #    "+", nets.ref.subsampled.aggregated) / length(nets.ref.subsampled.aggregated
    #    )
    #net.withRef.subsampled.aggregated.graph <- graph_from_adjacency_matrix(
    #    Mu.ref.subsampled.aggregated, 
    #    mode = "undirected", 
    #    weighted = TRUE
    #    )
    #net.withRef.subsampled.aggregated.graph <- set.vertex.attribute(
    #    graph = net.withRef.subsampled.aggregated.graph,
    #    name = "name",
    #    value = paired.datasets$genes
    #)
    return(list(
        net.withRef.subsampled.aggregated.graph, marker.genes.activity.plot, marker.activity.scores
        ))
}

xxx <- function(n, center, radius, increase, maxdist){
    fin <- TRUE
    xm <- matrix(ncol = 2, nrow = n)
    startp <- 0
    i <- 1
    z <- 1
    print(n)
    while(fin){
        nn <- calcular_numero_de_puntos_en_circunferencia(radius, maxdist)
        if(nn == 0) nn <- 1
        cat("Number of nodes in circunference", z, ":", nn, "\n")
        if(n-nn < 1){
            cat("Last circunference...\n")
            fin <- FALSE
            nn <- n
        }else{
            n <- n-nn
            }
        if(nn == 1 & fin){
            cat("Positionig as center(", center ,") as radius is", radius, "and maxdist is", maxdist, "in index:", i, "\n")
            xm[i,] <- center
            nn <- 1
        }else{
            f <- i+nn-1
            print(i:f)
            xm[i:f,] <- find_coords(center = center, radius = radius, n.nodes = nn, startp = startp)
        }
        radius <- radius + increase   
        startp <- startp + round(runif(1, min = 0, max = 360))
        i <- i+nn  
        z <- z+1 
    }
    return(xm)
}

calcular_numero_de_puntos_en_circunferencia <- function(radio, distancia_minima_entre_puntos) {
  # Calcula el diámetro de los puntos teniendo en cuenta la distancia mínima entre ellos
  diametro_punto <- distancia_minima_entre_puntos
  
  # Calcula el área de un círculo con el diámetro de los puntos
  area_circulo_punto <- pi * (diametro_punto / 2)^2
  
  # Calcula el área de la circunferencia objetivo
  area_circunferencia_objetivo <- pi * radio^2
  
  # Calcula el número de puntos que pueden ser acomodados dentro de la circunferencia
  num_puntos <- floor(area_circunferencia_objetivo / area_circulo_punto)
  
  return(num_puntos)
}


find_coords <- function(center = c(0,0), radius_var = FALSE, radius, n.nodes, startp) {
  
  # Calculate the radius based on the size of the square
  #n.nodes <- nrow(coords)
  cat("Number of nodes to place:", n.nodes, "\n")
  if(n.nodes > 10 & radius_var)
    radius <- runif(n = n.nodes, min = 0.1, max = radius)
  radius <- radius / sqrt(3)
  
  angle <- 360 / n.nodes

  # Calculate angles for each vertex
  angles <- seq((0+startp), (360+startp)-angle, by = angle)

  # Calculate the coordinates of each vertex
  x_coordinates <- radius * cos(pi * angles / 180)
  y_coordinates <- radius * sin(pi * angles / 180)

  # Shift the coordinates to the center
  x_coordinates <- x_coordinates + center[1]
  y_coordinates <- y_coordinates + center[2]
  


  # Combine x and y coordinates into a matrix
  coordinates <- matrix(c(x_coordinates, y_coordinates), ncol = 2)
  
  return(coordinates)
}

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

gmm <- function(edge.weight, maxit = 1000){
    stats <- boxplot.stats(edge.weight)
    outliers <- which(edge.weight %in% stats$out)
    noutliers <- length(outliers)
    cat('Number of outliers:', noutliers, '\n')
    if( noutliers == 0) outliers <- length(edge.weight)+1
    # decompose gaussian distribution
    
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

plot_dist <- function(dist.dat, plot.title, connected.group, median){

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
                  colour = "#2b9cff", lwd = 1) +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(
                    dist.dat$mu[unconnected.group], 
                    dist.dat$sigma[unconnected.group], 
                    lam = dist.dat$lambda[unconnected.group]),
                  colour = "#ff0000", lwd = 1) +
    ylab("Density") + 
    xlab("Meristem Single-Cell network edge weight") +
    labs(title=plot.title, subtitle = paste("t=", median)) +
    geom_vline(xintercept = median, linetype = "dotted", color = "green") +
    theme_classic() #+
    #theme(axis.title=element_text(size=18), title = element_text(size=22))

  return(plot)
}

plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}


normzero <- function(x){
    res <- rep(NA, length(x))
    inds <- which(!is.na(x))
    x <- as.numeric(na.omit(x))
    res[inds] <- (x-min(x))/(max(x)-min(x))
    return(res)
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
cell.types <- as.character(sce.reduced$)
if("Unknown" %in% cell.types) cell.types <- cell.types[cell.types != "Unknown"]

rm(sce.reduced)
rm(ACTIONet.out)
rm(W)
rm(H)
gc()

marker.genes <- readLines("others/plt_list.txt", warn = FALSE)

ref.net <- read.csv("output/Arabidopsis_thaliana/NETWORK/interaction.table.sc.csv")
ref.net <- ref.net[,1:3]
head(ref.net)
ref.net.graph <- as.matrix(ref.net[,c(1,2)])
ref.net.graph <- graph_from_edgelist(
    ref.net.graph, 
    directed = FALSE
    )
ref.net.graph
plot(ref.net.graph, layout = create_layout(ref.net.graph, layout = 'igraph', algorithm = 'kk'))


bc <- closeness(ref.net.graph)
vertices <- V(ref.net.graph)$name
plt <- vertices %in% marker.genes
table(plt)
plt.ids <- read.delim("others/plt.ids.txt", header = FALSE, col.names = c("name", "id"))
plt.tar <- rep(NA, length(vertices))
names(plt.tar) <- vertices
plt.tar[plt.ids$id] <- plt.ids$name

ref.nodes <- read.csv("others/prop_table.csv")
nrow(ref.nodes)
ref.nodes <- ref.nodes[!(ref.nodes$gen_id %in% plt.ids$id),]
head(ref.nodes,10)

ref.nodes$nReg <- sapply(ref.nodes$gen_id, function(target) nrow(ref.net[ref.net$target == target,]))
ref.nodes$regs <- sapply(ref.nodes$gen_id, function(target) 
    ref.net$source[ref.net$target == target]
)
ref.nodes$regs.name <- sapply(ref.nodes$regs, function(regs) 
    sort(plt.ids$name[plt.ids$id %in% regs])
)

plt.groups <- unique(ref.nodes$regs.name)
backbone.nodes <- plt.groups[order(sapply(plt.groups, length))]

ref.nodes$backbone.node <- sapply(ref.nodes$regs.name, function(plt.group){
    b <- TRUE
    i <- 1
    while(b){
        iden <- identical(plt.group, backbone.nodes[[i]])
        if(iden) b <- FALSE
        else i <- i+1
    }
    return(i)
})

backbone.ref.net <- match(ref.net$target, ref.nodes$gen_id)
backbone.ref.net <- ref.nodes$backbone.node[backbone.ref.net]
backbone.ref.net <- cbind(ref.net$source, backbone.ref.net)
nrow(backbone.ref.net)
backbone.ref.net <- na.omit(backbone.ref.net)
backbone.nodes.hub <- backbone.nodes[sapply(backbone.nodes, length) == 1]
head(backbone.ref.net)
nrow(backbone.ref.net)

backbone.ref.net.graph <- graph_from_edgelist(backbone.ref.net, directed = FALSE)

plot(backbone.ref.net.graph)

backbone.nodes.degree <- sapply(1:length(backbone.nodes), function(backbone.node.ind){
    backbone.node.n.targets <- nrow(ref.nodes[ref.nodes$backbone.node == backbone.node.ind, ])
})

backbone.nodes.degree <- normzero(backbone.nodes.degree)

ref.nodes$regs.ind <- sapply(ref.nodes$regs, function(regs){
    plt.ids$ind[plt.ids$id %in% regs]
})

plt.ids$ind <- sapply(plt.ids$id, function(plt) which(vertices == plt))

backbone.graph <- make_full_graph(length(backbone.nodes))
vertices <- V(ref.net.graph)$name
backbone.nodes.dist <- rescale(backbone.nodes.degree, to = c(0.05, 1))

backbone.edges <- as_edgelist(backbone.graph)
backbone.edges.bool <- apply(backbone.edges, 1, function(x){
    all(backbone.nodes[[x[1]]] %in% backbone.nodes[[x[2]]]) & length(backbone.nodes[[x[1]]]) == 1
})
table(backbone.edges.bool)
head(backbone.edges.bool,11)
backbone.edges[1:11,]
backbone.nodes[1:11]
backbone.edges <- backbone.edges[backbone.edges.bool,]
head(backbone.edges)
backbone.graph <- graph_from_edgelist(backbone.edges)

E(backbone.graph)$weight <- apply(backbone.edges, 1, function(nodes) backbone.nodes.degree[nodes[1]] + backbone.nodes.degree[nodes[2]]) 

smoothing((max(backbone.nodes.degree+1) - backbone.nodes.degree))
ggraph(backbone.graph, layout="stress", weights = smoothing((max(E(backbone.graph)$weight)+1) - E(backbone.graph)$weight, strength=0.8)) +
#ggraph(backbone.graph, layout="centrality", centrality = smoothing((max(backbone.nodes.degree+1) - backbone.nodes.degree), strength = 0.5)^10) +
#ggraph(backbone.graph, layout="kk") +
    #draw_circle(use = "focus",max.circle = 3) +
    geom_edge_link(width=0.2,colour="grey")+
    geom_node_point(aes(size = backbone.nodes.degree), shape = 21)+
    geom_node_text(aes(label = V(backbone.graph)),  colour = 'black', size=4,
                  show.legend = FALSE, family = "serif") +
    scale_size_continuous(range = c(10,20)) +
    #ylim(-100, 100) +
    #xlim(-100, 100) +
    theme_graph()

set.seed(15)
backbone.layout <- layout_with_stress(backbone.graph, weights = smoothing((max(E(backbone.graph)$weight)+1) - E(backbone.graph)$weight, strength=0.8))
backbone.layout <- layout_with_dh(backbone.graph)
backbone.layout <- layout_with_centrality(backbone.graph, smoothing((max(backbone.nodes.degree+1) - backbone.nodes.degree), strength = 0.5)^2)
backbone.layout <- apply(backbone.layout, 2, rescale, to = c(-10,10))

###############
layout <- matrix(ncol = 2, nrow = length(vertices))
for(i in 1:length(backbone.nodes)){
    cat("Working on backbone node:", i, "\n")
    plttt <- FALSE
    if(length(backbone.nodes[[i]]) == 1){
        pltt <- plt.ids$id[plt.ids$name == backbone.nodes[[i]]]
        cat("Adding",  pltt,"coordinates...\n")
        backbone.node.ind <- which(vertices == pltt)
        layout[backbone.node.ind,] <- backbone.layout[i,]
        plttt <- TRUE
    }
    group.targets <- ref.nodes$gen_id[ref.nodes$backbone.node == i]
    group.targets.inds <- which(vertices %in% group.targets)
    layout[group.targets.inds, ] <- xxx(
        n = nrow(layout[group.targets.inds, ,drop=FALSE]),
        center = backbone.layout[i,],
        #radius = ifelse(plttt, 0.33,0.18), increase = 0.2, maxdist = 0.35
        #radius = ifelse(plttt, 1.2,0.5), increase = 0.7, maxdist = 1.5
        radius = ifelse(plttt, 1,0.1), increase = 0.4, maxdist = 1
    )
}

layout[10,] <- c(10.5,-20)
layout[10,] <- c(7,6)

plt.ids$grp <- sapply(plt.ids$name, function(plt) which(plt == backbone.nodes))
plt.ids$grp[5] <- 38

V(ref.net.graph)$grp <- sapply(vertices, function(node){
    if(node %in% ref.nodes$gen_id){
        plts <- ref.nodes$regs.name[ref.nodes$gen_id == node][[1]]
        p <- sample(plts, 1)
    }else{
        p <- plt.ids$name[plt.ids$id == node]
    }
    return(p)
})

ggraph(ref.net.graph, layout="manual", x=layout[,1], y=layout[,2]) +
    geom_mark_hull(
            aes(x, y, fill = grp, color = grp), size=0,
            concavity = 4,
            expand = unit(2, "mm"),
            alpha = 0.1, show.legend = TRUE
        ) +
    geom_edge_link(aes(col = gmmresult[[2]], alpha = gmmresult[[2]]), width=0.2)+
    geom_node_point(aes(fill = factor(plt.tar), size = factor(plt), alpha = factor(plt)), shape = 21, color = "black")+
    scale_fill_manual(
        breaks = c ("PLT1", "PLT2", "PLT3", "PLT4", "PLT5", "PLT7"), 
        values = brewer.pal(6, name = "Dark2"),
        name = "", na.value = "#343434") + 
    scale_color_manual(
        breaks = c ("PLT1", "PLT2", "PLT3", "PLT4", "PLT5", "PLT7"), 
        values = c("#1b9e772d", "#d95f021d", "#7570b329", "#e7298b34", "#67a61e23", "#e6a90246"),
        name = "", na.value = "#343434") + 
    scale_edge_color_gradient2(
        low = "#ff0000",
        mid = "#4d80e5c9",
        high = "#2b9cff",
        midpoint = gmmresult[[3]], name="Meristem SC\nedge weight"
    ) +
    scale_edge_alpha(range = c(0.2,1), guide = "none") +
    scale_size_discrete(range = c(1,5), guide = "none") +
    scale_alpha_discrete(range = c(0.9,1), guide = "none") +
    guides(fill = guide_legend(override.aes = list(size = 5))) +
    theme_graph(fg_text_colour = 'white', base_family = 'Helvetica') 

ggsave("ok4.pdf")

out.path <- "output/Arabidopsis_thaliana/SC"

cell.types <- sce.reduced$time.anno
cell.types <- sce.reduced$final.anno
table(cell.types)

cell.type <- "Columella"

markers.activity <- data.frame()

for(cell.type in ct){

    sc.results <- build_sc_network(
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
    
    markers.exp <- sc.results[[3]]

    markers.exp[['Cell.type']] <- cell.type

    head(markers.exp)

    nrow(sc.results[[3]])

    mex <- rbind(markers.exp, markers.exp1)

    markers.activity <- rbind(markers.activity, markers.exp)

    nrow(markers.activity)

    markers.activity %>% mutate(Cell.type = fct_reorder(Cell.type, activity)) %>%
    ggplot(
        aes(x=activity, y=Cell.type, fill=gene)) +
        geom_density_ridges(alpha = 0.85, size=0.2) +
        scale_fill_brewer(palette = "Dark2", labels = plt.ids$name, name = "") +
        theme_ridges()


    ggsave("cell.types.activity.pdf")

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

ref.net <- read.csv("output/Arabidopsis_thaliana/NETWORK/interaction.table.sc.csv")

ref.net.weights <- ref.net[,-c(1:3)]
dim(ref.net.weights)
ref.net.weights[1:6,1:4]

ref.net.weights.scaled <- apply(ref.net.weights, 2, normzero)
dim(ref.net.weights.scale)
ref.net.weights.scaled[1:6,1:4]

ref.net.scaled <- cbind(ref.net[,1:3], ref.net.weights.scaled)
write.csv(ref.net.scales, "interaction_table.centered.csv", quote=FALSE, row.names = FALSE, na = "")


decompose_gmm <- function(dat, network.name, strength = 0.25){

    cat("Decomposing network:", network.name, "\n")

    dat.smoothing.centered.na <- rep(NA, length(dat))

    dat.na <- which(is.na(dat))
    dat.non.na <- which(!is.na(dat))

    dat <- as.numeric(na.omit(dat))

    cat("Using streng =", strength, "for smoothing\n")

    dat.smoothing <- smoothing(dat, strength = strength)

    dat.smoothing.centered <- normzero(dat.smoothing)

    dat.smoothing.centered.na[dat.na] <- 0
    dat.smoothing.centered.na[dat.non.na] <- dat.smoothing.centered

    x <- gmm(dat.smoothing.centered)
    xm1 <- median(x$probabilities$edge.weight[x$probabilities$comp.1 >= 0.5])
    xm2 <- median(x$probabilities$edge.weight[x$probabilities$comp.2 >= 0.5])

    comp <- ifelse(xm1>xm2, "comp.1", "comp.2")

    if(xm1>xm2) xmed <- quantile(x$probabilities$edge.weight[x$probabilities$comp.1 >= 0.5])[1]
    else xmed <- quantile(x$probabilities$edge.weight[x$probabilities$comp.2 >= 0.5])[1]

    gmm.plot <- plot_dist(x, 
        plot.title="",
        connected.group=comp,
        median = round(xmed,3)
        )
    
    return(list(gmm.plot, dat.smoothing.centered.na, xmed))
}

strength.min <- 0
strength.max <- 3.5

ref.net.weights.scaled.var <- apply(ref.net.weights.scaled, 2, var, na.rm = TRUE)

scaling.factor <- (strength.max - strength.min) / max(ref.net.weights.scaled.var) - min(ref.net.weights.scaled.var)

ref.net.weights.scaled.var <- (ref.net.weights.scaled.var - min(ref.net.weights.scaled.var)) * scaling.factor + strength.min

for(set in 1:ncol(ref.net.weights)){

    network.name <- names(ref.net.weights)[set]
    gmmresult1 <- decompose_gmm(ref.net.weights[,set], network.name, strength = ref.net.weights.scaled.var[set])

    gmmresult1[[1]]

    gmmresult[[2]]

    gmmresult[[3]]

    ggsave(gmmplot, filename = file.path(out.path, paste0(network.name, ".gmm.svg")))

}

ggsave("gmm.ok.pdf")

plot(ref.net.graph.)


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