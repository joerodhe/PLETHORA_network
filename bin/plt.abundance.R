suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))

option_list = list(
    make_option(c("-p", "--plts"), type = "character", default = NULL,
              help = "PLTs table ids,", metavar = "character"),
    make_option(c("-r", "--rpkm"), type = "character", default = NULL,
              help = "RPKM matrix expression data.", metavar = "character"),  
    make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output svg plot file.", metavar = "character")           
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

plt.ids <- read.delim(opt$plts, header = FALSE , col.names = c('plt.name', 'plt.id'))

rpkm.matrix <-read.delim(opt$rpkm, check.names = FALSE)
rpkm.matrix <- as.data.frame(t(rpkm.matrix))

genes <- colnames(rpkm.matrix)
rpkm.matrix$zone <- rownames(rpkm.matrix)
rpkm.reshaped <- reshape(rpkm.matrix,
                    direction = "long", 
                    varying = list(genes),
                    v.names = "RPKM",
                    idvar = "zone",
                    timevar = "Gene", 
                    times = genes)
rpkm.reshaped$zone <- sapply(rpkm.reshaped$zone, function(x) strsplit(x, '-')[[1]][1])
plt.rpkm <- rpkm.reshaped %>%
  filter(Gene %in% plt.ids$plt.id)
plt.rpkm$plt <- sapply(plt.rpkm$Gene, function(x) plt.ids$plt.name[plt.ids$plt.id == x])
plt.rpkm <- plt.rpkm[with(plt.rpkm, order(plt)),]
  
plt.rpkm.plot <-  ggplot(plt.rpkm, aes(x=plt, y=log2(RPKM+1), color=zone)) +
  geom_boxplot() +
  geom_point(alpha = 0.8) +
  #coord_flip() +
  theme(axis.ticks.x = element_blank(), legend.box = "horizontal") +
  labs(color = "Root zone", x="", y=expression("log"[2]*"(RPKM+1)")) +
  theme_classic() +
  scale_color_brewer(palette = "Paired")

ggsave(plot = plt.rpkm.plot, filename = opt$output)
