suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))

option_list = list(
    make_option(c("-f", "--files"), type = "character", default = NULL,
              help = "RDS files with pca data. (comma-separated)", metavar = "character"),
    make_option(c("-o", "--organisms"), type = "character", default = NULL,
              help = "Organism list (comma-separated)", metavar = "character"),
    make_option(c("-p", "--plot"), type = "character", default = NULL,
              help = "Plot output file.", metavar = "character")         
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

files <- strsplit(opt$files, ",")[[1]]
organisms <- strsplit(opt$organisms, ",")[[1]]

organisms <- c(
    "Arabidopsis_thaliana",
    "Solanum_lycopersicum",
    "Cucumis_sativus",
    "Oryza_sativa_Japonica_Group",
    "Glycine_max",
    "Zea_mays"
)

files <- file.path("temp", organisms, "pca.RDS")

pca.dat <- as.data.frame(do.call(rbind, lapply(files, function(f) readRDS(f)$x)))

pca.dat$Sample <- as.factor(sapply(rownames(pca.dat), function(sample.id) strsplit(sample.id, "-")[[1]][1]))
pca.dat$Organism <- as.factor(rep(organisms, each=9))

for(organism in organisms){
    pca.dat$alpha <- organism == pca.dat$Organism
    ggplot(pca.dat, aes(x = PC1, y = PC2, color = Organism, shape = Sample, alpha=alpha)) +
        geom_point(size=3) +
        scale_alpha_discrete(range = c(0.1, 1)) +
        theme_classic() +
        labs(title = "Principal Component Analysis", 
            x = "PC1", 
            y = "PC2", 
            color = "Organism", 
            shape = "Root zone"
            )

    ggsave(file.path("temp", paste0("pca.", organism, ".svg")), width = 6, height = 4)
}
