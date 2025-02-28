suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))

option_list = list(
    make_option(c("-p", "--plts"), type = "character", default = NULL,
                help = "PLTs table ids (two-column text file: plt.name, plt.id)", 
                metavar = "character"),
    make_option(c("-r", "--rpkm"), type = "character", default = NULL,
                help = "RPKM matrix expression data file path", 
                metavar = "character"),  
    make_option(c("-o", "--output"), type = "character", default = NULL,
                help = "Output SVG plot file path", 
                metavar = "character")           
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

###############################################################################
# Read PLT IDs (mapping PLT name to gene ID)
###############################################################################
plt.ids <- read.delim(opt$plts, header = FALSE , col.names = c('plt.name', 'plt.id'))

###############################################################################
# Read and transpose the RPKM expression matrix
###############################################################################
rpkm.matrix <- read.delim(opt$rpkm, check.names = FALSE)
rpkm.matrix <- as.data.frame(t(rpkm.matrix))

# Save a vector of gene names for reshaping
genes <- colnames(rpkm.matrix)

###############################################################################
# Add zone information (extracted from rownames) to the expression data
###############################################################################
rpkm.matrix$zone <- rownames(rpkm.matrix)

###############################################################################
# Reshape from wide to long format, to have columns: zone, Gene, RPKM
###############################################################################
rpkm.reshaped <- reshape(
  rpkm.matrix,
  direction = "long", 
  varying = list(genes),
  v.names = "RPKM",
  idvar = "zone",
  timevar = "Gene", 
  times = genes
)

###############################################################################
# Simplify the 'zone' field by splitting at '-' and retaining the first part
###############################################################################
rpkm.reshaped$zone <- sapply(rpkm.reshaped$zone, function(x) strsplit(x, '-')[[1]][1])

###############################################################################
# Filter the reshaped data to retain only the PLT IDs of interest
###############################################################################
plt.rpkm <- rpkm.reshaped %>%
  filter(Gene %in% plt.ids$plt.id)

###############################################################################
# Map the PLT gene IDs to their readable names
###############################################################################
plt.rpkm$plt <- sapply(plt.rpkm$Gene, function(x) plt.ids$plt.name[plt.ids$plt.id == x])

###############################################################################
# Order the data by PLT name for a consistent plot ordering
###############################################################################
plt.rpkm <- plt.rpkm[with(plt.rpkm, order(plt)),]

###############################################################################
# Generate a box plot of log2(RPKM+1) for each PLT gene, split by root zone
###############################################################################
plt.rpkm.plot <- ggplot(plt.rpkm, aes(x=plt, y=log2(RPKM+1), color=zone)) +
  geom_boxplot() +
  geom_point(alpha = 0.8) +
  theme(axis.ticks.x = element_blank(), legend.box = "horizontal") +
  labs(color = "Root zone", x="", y=expression("log"[2]*"(RPKM+1)")) +
  theme_classic() +
  scale_color_brewer(palette = "Paired")

###############################################################################
# Save the plot as an SVG to the specified output path
###############################################################################
ggsave(plot = plt.rpkm.plot, filename = opt$output)
