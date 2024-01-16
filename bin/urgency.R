library(ggplot2)

Zs <- t(data.frame(
    'MZ_Zea'=4.9,
    'MZ_Solanum'=27.8,
    'MZ_Oryza'=11.2,
    'MZ_Glycine'=4.5,
    'MZ_Cucumis'=12.1,
    'EZ_Zea'=9.2,
    'EZ_Solanum'=8.3,
    'EZ_Oryza'=2.3,
    'EZ_Glycine'=9.1,
    'EZ_Cucumis'=9.9,
    'EZ_Arabidopsis'=29.2,
    'DZ_Zea'=6.2,
    'DZ_Solanum'=9.3,
    'DZ_Oryza'=10.1,
    'DZ_Glycine'=7.2,
    'DZ_Cucumis'=9.7,
    'DZ_Arabidopsis'=16.9
))

colnames(Zs) <- "Zsummary"

Zones <- sapply(rownames(Zs), function(x) strsplit(x, '_')[[1]][1])
Organism <- sapply(rownames(Zs), function(x) strsplit(x, '_')[[1]][2])

Zs <- as.data.frame(cbind(Zs, Zones, Organism))
Zs$Network <- rownames(Zs)
Zs$Zsummary <- as.numeric(Zs$Zsummary)

for(organism in Organism){
    Zs$alpha <- organism == Organism
ggplot(Zs, aes(x=Zones, y=Zsummary, fill=Organism, alpha = alpha)) +
    geom_bar(stat="identity", position = "dodge",  color = "grey") +
    geom_hline(yintercept = 10, linetype = "dashed", color="red") +
    geom_hline(yintercept = 2, linetype = "dashed", color="orange") +
    labs(title = "Zsummary. Reference: Arabidopsis MZ network") +
    scale_alpha_discrete(range = c(0.1, 1)) +
    theme_classic()

ggsave(file.path("temp", paste0(organism, ".zsum.svg")))

}


blast.files <- 