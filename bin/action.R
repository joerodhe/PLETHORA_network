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

###############################################################################
# Function: build_sc_network
# --------------------------
# Description:
#   Constructs a single-cell network using a reference network (ref.net),
#   gene activity matrix (A), and annotation of cell types. Subsamples gene
#   activity scores for a given cell.type, pairs them with the reference network,
#   and computes an aggregated network. Additionally, plots marker gene
#   activities (density plot).
#
# Params:
#   ref.net       : igraph object of the reference network.
#   A             : Gene activity matrix (W %*% H from ACTIONet).
#   annotation.list: A vector specifying the cell type of each column in A.
#   cell.type     : The cell type for which the single-cell network is built.
#   genes         : A vector of gene names corresponding to rows of A.
#   marker.genes  : A subset of genes of interest (e.g. markers).
#   threads       : Number of threads to use for parallel computations.
#
# Returns:
#   A list containing:
#     1) An igraph object of the aggregated SC network,
#     2) A ggplot object (density plot of marker gene activities),
#     3) A data frame of marker gene activity scores.
###############################################################################
build_sc_network <- function(ref.net, 
                             A, 
                             annotation.list, 
                             cell.type, 
                             genes, 
                             marker.genes,
                             threads = 8){
  cat("Working on building networks for", cell.type, "cells...\n")
  cat("Subsampling gene activity scores...\n")
  
  # Identify columns (cells) that match the target cell type
  cols <- which(annotation.list == cell.type)
  cat("Number of cells:", length(cols), "\n")
  
  # Subsample gene activities for the selected cells
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
  # Filter to marker genes
  marker.rows <- match(marker.genes, genes)
  marker.activity.scores <- t(subsampled.activity.scores[marker.rows, ])
  colnames(marker.activity.scores) <- marker.genes
  
  # Reshape for ggplot
  marker.activity.scores <- reshape2::melt(marker.activity.scores)
  colnames(marker.activity.scores) <- c('id', 'gene', 'activity')
  
  net.withRef.subsampled.aggregated.graph <- NA
  marker.genes.activity.plot <- NA
  
  # Plot density of marker gene activities
  marker.genes.activity.plot <- ggplot(
    marker.activity.scores, 
    aes(x=activity, fill=gene)
  ) +
    geom_density()
  
  cat("Building paired datasets object...\n")
  # Pair the subsampled activity scores with the reference network
  paired.datasets <- pair.datasets(
    ref.net, 
    subsampled.activity.scores
  )
  
  # Convert igraph network to TsparseMatrix adjacency
  ref.net.adj <- as(get.adjacency(
    paired.datasets$net), 'TsparseMatrix'
  )
  
  cat("Building SC networks...\n")
  # Build single-cell networks for each subsample
  nets.ref.subsampled.aggregated <- construct_cell_networks_summary(
    net = ref.net.adj, 
    gene_activities = paired.datasets$activity.scores, 
    thread_no = threads, 
    total_subsamples = 100, 
    cells_per_subsample = 500, 
    seed = 0
  )
  
  cat("Aggregating networks...\n")
  # Aggregate across subsamples
  Mu.ref.subsampled.aggregated <- Reduce(
    "+", nets.ref.subsampled.aggregated) / length(nets.ref.subsampled.aggregated
    )
  
  # Convert to igraph
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
  
  # Return the aggregated network, marker plot, and activity scores
  return(list(
    net.withRef.subsampled.aggregated.graph,
    marker.genes.activity.plot,
    marker.activity.scores
  ))
}

###############################################################################
# Function: arr_points (arbitrary naming)
# --------------------------------
# Description:
#   Places n points around concentric circles (one circle at a time).
#   Continues adding circles (increasing radius each time) until all points
#   are placed. Used for layout purposes in a custom graph plotting scenario.
#
# Params:
#   n           : Number of points to place.
#   center      : (x, y) center coordinate.
#   radius      : Initial radius of the first circle.
#   increase    : Radius increment for each subsequent circle.
#   maxdist     : Minimum distance between points (used in an area-based heuristic).
#
# Returns:
#   A matrix of (x, y) coordinates for each of the n points placed.
###############################################################################
arr_points <- function(n, center, radius, increase, maxdist){
  fin <- TRUE
  xm <- matrix(ncol = 2, nrow = n)
  startp <- 0
  i <- 1
  z <- 1
  print(n)
  
  # Place points in concentric circles until all are placed
  while(fin){
    nn <- calcular_numero_de_puntos_en_circunferencia(radius, maxdist)
    if(nn == 0) nn <- 1
    cat("Number of nodes in circunference", z, ":", nn, "\n")
    
    # If adding nn would exceed the total needed, place the remainder and stop
    if(n-nn < 1){
      cat("Last circunference...\n")
      fin <- FALSE
      nn <- n
    }else{
      n <- n-nn
    }
    
    # If the circle has only 1 node to place, put it at the center
    if(nn == 1 & fin){
      cat("Positionig as center(", center ,") as radius is", radius, "and maxdist is", maxdist, "in index:", i, "\n")
      xm[i,] <- center
      nn <- 1
    }else{
      f <- i+nn-1
      print(i:f)
      xm[i:f,] <- find_coords(center = center, radius = radius, n.nodes = nn, startp = startp)
    }
    
    # Increase radius for the next circle
    radius <- radius + increase   
    startp <- startp + round(runif(1, min = 0, max = 360))
    i <- i+nn  
    z <- z+1 
  }
  return(xm)
}

###############################################################################
# Function: calcular_numero_de_puntos_en_circunferencia
# -----------------------------------------------------
# Description:
#   Based on a circle's area and a minimal distance between points, estimates
#   how many points can fit inside a circle. Uses area-based heuristics (not
#   exact geometry).
#
# Params:
#   radio                      : Radius of the circle.
#   distancia_minima_entre_puntos : Minimal desired distance between points.
#
# Returns:
#   An integer estimate of how many points can be placed in the circle.
###############################################################################
calcular_numero_de_puntos_en_circunferencia <- function(radio, distancia_minima_entre_puntos) {
  diametro_punto <- distancia_minima_entre_puntos
  area_circulo_punto <- pi * (diametro_punto / 2)^2
  area_circunferencia_objetivo <- pi * radio^2
  num_puntos <- floor(area_circunferencia_objetivo / area_circulo_punto)
  return(num_puntos)
}

###############################################################################
# Function: find_coords
# ---------------------
# Description:
#   Distributes 'n.nodes' points evenly around a circle of given 'radius' and 
#   offsets their angles by 'startp'. Optionally can vary the radius if 
#   radius_var is TRUE.
#
# Params:
#   center       : (x, y) for the circle center.
#   radius_var   : Boolean for whether the radius differs for each point.
#   radius       : The circle's radius (or max radius, if radius_var is TRUE).
#   n.nodes      : Number of points to place.
#   startp       : Angle offset in degrees.
#
# Returns:
#   A matrix with columns = 2 (x, y) for each point placed.
###############################################################################
find_coords <- function(center = c(0,0), radius_var = FALSE, radius, n.nodes, startp) {
  
  cat("Number of nodes to place:", n.nodes, "\n")
  if(n.nodes > 10 & radius_var)
    radius <- runif(n = n.nodes, min = 0.1, max = radius)
  
  radius <- radius / sqrt(3)
  angle <- 360 / n.nodes
  angles <- seq((0+startp), (360+startp)-angle, by = angle)
  
  # Compute x,y positions from angles
  x_coordinates <- radius * cos(pi * angles / 180)
  y_coordinates <- radius * sin(pi * angles / 180)
  
  # Shift to the given center
  x_coordinates <- x_coordinates + center[1]
  y_coordinates <- y_coordinates + center[2]
  
  coordinates <- matrix(c(x_coordinates, y_coordinates), ncol = 2)
  return(coordinates)
}

###############################################################################
# Function: construct.archetype.signature.profile
# -----------------------------------------------
# Description:
#   Given a SingleCellExperiment object with reducedDims(sce)[["S_r"]] from
#   ACTIONet, reconstruct the full gene x archetype matrix. Returns a matrix
#   of size (genes x archetypes).
#
# Params:
#   sce             : A SingleCellExperiment containing rowData(sce)$rotation 
#                     and reducedDims(sce)[["S_r"]].
#   ACTIONet.out    : Output from run.ACTIONet, specifically reconstruct.out 
#                     which has C_stacked, etc.
#   reduction_slot  : Name of the reduced dimension slot, default "S_r".
#
# Returns:
#   archetype.signature.profile: A gene x archetype matrix.
###############################################################################
construct.archetype.signature.profile <- function(sce, ACTIONet.out, reduction_slot = "S_r") {
  
  # Multiply the SC reduced dimension (S_r) by the stacked coefficient matrix (C_stacked)
  reduced.archetype.profile = (t(reducedDims(sce)[[reduction_slot]]) %*% ACTIONet.out$reconstruct.out$C_stacked)
  
  V <- rowData(sce)$rotation
  rownames(V) <- rownames(rowData(sce))
  
  # Re-order columns by numeric PC index
  perm = order(sapply(colnames(V), function(str) as.numeric(stringr::str_replace(str, "PC", ""))))
  V = V[, perm]
  
  # Multiply V by the reduced archetype profile to get gene x archetype matrix
  archetype.signature.profile = V %*% reduced.archetype.profile
  rownames(archetype.signature.profile) = rownames(sce)
  
  return(archetype.signature.profile)
}

###############################################################################
# Function: build_sc_net
# ----------------------
# Description:
#   Similar to build_sc_network, constructs a single-cell (SC) network based on
#   a reference network and gene activity matrix, but returning only the network
#   and a marker gene plot (no raw data frame). Possibly an older or alternate
#   version.
###############################################################################
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

###############################################################################
# Function: gmm
# -------------
# Description:
#   Fits a 2-component Gaussian Mixture Model to a vector of edge weights, 
#   identifies outliers using boxplot.stats, and returns the mixture parameters
#   plus posterior probabilities for each edge weight.
#
# Params:
#   edge.weight : Named numeric vector of edge weights.
#   maxit       : Maximum iterations for the EM algorithm in normalmixEM.
#
# Returns:
#   A list containing the mixture model parameters (mu, lambda, sigma) and
#   a data frame "probabilities" with the raw edge.weight and posterior comps.
###############################################################################
gmm <- function(edge.weight, maxit = 1000){
  stats <- boxplot.stats(edge.weight)
  outliers <- which(edge.weight %in% stats$out)
  noutliers <- length(outliers)
  cat('Number of outliers:', noutliers, '\n')
  
  if( noutliers == 0) outliers <- length(edge.weight)+1
  
  cat('Performing GMM with', maxit, 'max iterations...\n')
  gaussian.dist <- normalmixEM(edge.weight[-outliers], k=2, maxit=maxit)
  cat('Done...\n')
  
  probs <- as.data.frame(
    gaussian.dist$posterior,
    row.names = names(edge.weight)[-outliers]
  )
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
  
  cat('Building dataframe...\n')
  decomposed.dist <- cbind(edge.weight, probs)
  
  return(list(
    mu=gaussian.dist$mu, 
    lambda=gaussian.dist$lambda,
    sigma=gaussian.dist$sigma, 
    probabilities=decomposed.dist))
}

###############################################################################
# Function: plot_group_proportion
# -------------------------------
# Description:
#   Creates a bar chart showing the proportion of each group label ("Connected" 
#   vs "NotConnected", or any assigned "group") in a decomposed distribution.
###############################################################################
plot_group_proportion <- function(dist.dat, observation.type, plot.title){
  plot <- ggplot(dist.dat, aes(x=factor(group))) +
    geom_bar() +
    xlab(element_blank()) +
    ylab(paste("Number of", observation.type)) +
    labs(title=plot.title) +
    theme_classic()
  
  return(plot)
}

###############################################################################
# Function: plot_dist
# -------------------
# Description:
#   Plots a histogram of edge weights with two fitted Gaussian mixture curves
#   overlaid (one for the "connected" group, one for "unconnected").
#
# Params:
#   dist.dat       : A list from gmm() containing mixture parameters and posterior.
#   plot.title     : Title for the plot.
#   connected.group: Which mixture component is considered "connected" (comp.1 or comp.2).
#   median         : A numeric value to show as a vertical dotted line on the plot.
###############################################################################
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
    theme_classic()
  
  return(plot)
}

###############################################################################
# Function: plot_mix_comps
# ------------------------
# Description:
#   Helper function for stat_function() calls in plot_dist(), returning the 
#   density of a normal distribution scaled by lam (mixture proportion).
###############################################################################
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

###############################################################################
# Function: normzero
# ------------------
# Description:
#   Normalizes a numeric vector to the [0,1] range by subtracting the min and
#   dividing by max-min (ignoring NA). Returns a vector of same length as input.
###############################################################################
normzero <- function(x){
  res <- rep(NA, length(x))
  inds <- which(!is.na(x))
  x <- as.numeric(na.omit(x))
  res[inds] <- (x-min(x))/(max(x)-min(x))
  return(res)
}

###############################################################################
# Main Script Execution
# ---------------------
# 1. Read a SingleCellExperiment object.
# 2. Run reduce.sce(), run.ACTIONet(), and derive a signature profile (W, H, A).
# 3. Clear memory of unneeded objects.
# 4. Load a reference network, read marker genes, etc.
# 5. Demonstrate custom layout logic for backbone or PLT nodes.
# 6. Construct single-cell networks for multiple cell types and save results.
# 7. Perform GMM-based decomposition of edge weights.
# (Commented code at the end shows further exploratory steps.)
###############################################################################
sce <- readRDS('others/sce.counts.RDS')

# Reduce SCE to 50 dims and store rotation (for constructing V in the archetype profile)
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

# Construct the archetype signature profile
ACTIONet.out$signature.profile <- construct.archetype.signature.profile(sce.reduced, ACTIONet.out)

saveRDS(ACTIONet.out, file = "others/actinet.out.RDS")

# Load results, build the gene activity matrix A = W %*% H
ACTIONet.out <- readRDS("others/actinet.out.RDS")
W <- ACTIONet.out$signature.profile
H <- ACTIONet.out$reconstruct.out$H_stacked
A <-  W %*% H

# Clean up environment
sce.reduced <- readRDS("others/sce.reduced.RDS")
genes <- rownames(sce.reduced)
cell.types <- as.character(sce.reduced$)  # The user might specify a column name
if("Unknown" %in% cell.types) cell.types <- cell.types[cell.types != "Unknown"]

rm(sce.reduced)
rm(ACTIONet.out)
rm(W)
rm(H)
gc()

marker.genes <- readLines("others/plt_list.txt", warn = FALSE)

# Load a reference network (3-column table: source, target, weight)
ref.net <- read.csv("output/Arabidopsis_thaliana/NETWORK/interaction.table.sc.csv")
ref.net <- ref.net[,1:3]
head(ref.net)

# Convert to igraph
ref.net.graph <- as.matrix(ref.net[,c(1,2)])
ref.net.graph <- graph_from_edgelist(
  ref.net.graph, 
  directed = FALSE
)
ref.net.graph
plot(ref.net.graph, layout = create_layout(ref.net.graph, layout = 'igraph', algorithm = 'kk'))

# (Several lines of code generating layout logic and computing closeness, 
#  building a "backbone" network, normalizing degrees, etc.)

###############################################################################
# Example of building and saving single-cell networks for each cell type
###############################################################################
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

###############################################################################
# GMM-based decomposition of SC edge weights across conditions or cell types
###############################################################################
ref.net <- read.csv("output/Arabidopsis_thaliana/NETWORK/interaction.table.sc.csv")

ref.net.weights <- ref.net[,-c(1:3)]
dim(ref.net.weights)
ref.net.weights[1:6,1:4]

# Normalize each column to 0..1
ref.net.weights.scaled <- apply(ref.net.weights, 2, normzero)
dim(ref.net.weights.scale)
ref.net.weights.scaled[1:6,1:4]

ref.net.scaled <- cbind(ref.net[,1:3], ref.net.weights.scaled)
write.csv(ref.net.scales, "interaction_table.centered.csv", quote=FALSE, row.names = FALSE, na = "")

###############################################################################
# Function: decompose_gmm
# -----------------------
# Description:
#   Wrapper to apply GMM decomposition (via gmm function) to a vector of 
#   (possibly scaled or smoothed) edge weights, identify connected group,
#   and produce a mixture plot.
###############################################################################
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

ggsave("gmm.pdf")

plot(ref.net.graph.)

###############################################################################
# Additional commented out code for exploration (constructing networks for 
# "meristem" cells, etc.)
###############################################################################
# (Remains unchanged from the original script)



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
