suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.At.tair.db))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(RColorBrewer))

option_list = list(
    make_option(c("-c", "--compendium"), type = "character", default = NULL,
              help = "Santurari's PLT-genes-regulated compendium", metavar = "character"),
    make_option(c("-d", "--degs"), type = "character", default = NULL,
              help = "MZ DEGs list", metavar = "character"),
    make_option(c("-t", "--tfbs"), type = "character", default = NULL,
              help = "Matrix-scan rsat results of PLT TFBS mz-degs-promoter-scanning", metavar = "character"),
    make_option(c("-f", "--tfs"), type = "character", default = NULL,
              help = "Arabidopsis tfs ids.", metavar = "path"),
    make_option(c("-p", "--plts"), type = "character", default = NULL,
              help = "Arabidopsis PLTs ids.", metavar = "path"),
    make_option(c("-e", "--edges"), type = "character", default = NULL,
              help = "Interaction table output file.", metavar = "path"),
    make_option(c("-n", "--nodes"), type = "character", default = NULL,
              help = "Nodes table output file.", metavar = "path"),
    make_option(c("-g", "--plots"), type = "character", default = NULL,
              help = "Plots path.", metavar = "path")               
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

opt$compendium <- "input/Arabidopsis_thaliana/PLTs/TPC2016-00656-RAR2_Supplemental_Data_Set_1.xlsx"
opt$degs <- "output/Arabidopsis_thaliana/DEA/upregulated.DEGs.txt"
opt$tfbs <- "output/Arabidopsis_thaliana/PLTs/degs.plt.compendium.ft"
opt$tfbs <- "output/Arabidopsis_thaliana/PLTs/mz.upregulated.degs.plt_regulated.ft"
opt$plts <- 'others/plt.ids.txt'
opt$tfs <- "input/Arabidopsis_thaliana/GENOME/TAIR10.TF_list.txt"
opt$edges <- "output/Arabidopsis_thaliana/NETWORK/network_table.csv"
opt$nodes <- "output/Arabidopsis_thaliana/NETWORK/prop_table.csv"
opt$plots <- "output/Arabidopsis_thaliana/NETWORK"

compendium <- read_excel(opt$compendium,
  sheet = "Compendium", 
  range = "A8:AB20534", 
  na = "NA", 
  col_names = c(
    "AGI", 
    "ProbeID", 
    "Gene Symbol", 
    "Description", 
    "PLT1_FC", 
    "PLT2_FC", 
    "PLT3_FC", 
    "PLT4_FC", 
    "PLT5_FC", 
    "PLT7_FC", 
    "QC_1h_FC", 
    "QC_4h_FC", 
    "PLT1_p", 
    "PLT2_p", 
    "PLT3_p", 
    "PLT4_p", 
    "PLT5_p", 
    "PLT7_p", 
    "QC_1h_p", 
    "QC_4h_p", 
    "PLT1", 
    "PLT2", 
    "PLT3", 
    "PLT4", 
    "PLT5", 
    "PLT7", 
    "QC_1h_plt", 
    "QC_4h_plt"))

compendium.regulated.gen <- apply(compendium[,21:26], MARGIN = 1, function(gen) 1%in%gen | -1%in%gen)
compendium.plt.targets <- compendium$AGI[compendium.regulated.gen]

mz.degs <- readLines(opt$degs)
plt.mz.degs.tfbs <- read.delim(opt$tfbs, comment.char = ';')
plt.mz.degs.tfbs.genes <- unique(plt.mz.degs.tfbs$X.seq_id)

pdf(file.path(opt$plots, "ntw.venn.pdf"))
grid.newpage()
draw.triple.venn(area1 = length(compendium.plt.targets), 
                 area2 = length(mz.degs), 
                 area3= length(plt.mz.degs.tfbs.genes), 
                 n12 = length(intersect(compendium.plt.targets, mz.degs)), 
                 n13 = length(intersect(compendium.plt.targets, plt.mz.degs.tfbs.genes)), 
                 n23 = length(intersect(mz.degs, plt.mz.degs.tfbs.genes)),
                 n123 = (length(
                   intersect(
                     intersect(compendium.plt.targets, mz.degs),
                     plt.mz.degs.tfbs.genes))),
                 category = c("Compendium of\nPLT regulated genes", 
                              "Meristem\nupregulated DEGs", 
                              "Meristem upregulated DEGs\nwith PLT binding motif"),
                 scaled = T, 
                 fill = c("#1B9E77", "#D95F02", "#E6AB02"), 
                 lty = "blank", cat.cex = 0.9, cat.pos = 180
                 )
dev.off()

c <- lapply(plt.mz.degs.tfbs.genes, function(query){
  temp <- plt.mz.degs.tfbs[plt.mz.degs.tfbs$X.seq_id == query,]
  return(temp[which.min(temp$Pval),])
})
filt_gen_ms_sant <- c[[1]]
for(row in c[2:length(c)]){
  filt_gen_ms_sant <- rbind(filt_gen_ms_sant, row)
}

santuari_evidence <- sapply(filt_gen_ms_sant$X.seq_id, function(gen){
  regulated <- compendium[compendium$AGI == gen,]
  plts <- paste(na.omit(sapply(colnames(regulated), function(plt){
    if(is.na(regulated[,plt])) return(NA)
    else if(regulated[,plt] == 1) return(plt)
    else if(regulated[,plt] == -1) return(paste(plt, "(repressive)"))
    else return(NA)
  })), collapse = ";")
  return(plts)
})

filt_gen_ms_sant <- cbind(filt_gen_ms_sant, santuari_evidence)


ego_f <- enrichGO(
    gene = filt_gen_ms_sant$X.seq_id, 
    keyType = 'TAIR', 
    OrgDb = org.At.tair.db, 
    ont = "BP", 
    pAdjustMethod = "BH", 
    pvalueCutoff = 0.05, 
    qvalueCutoff = 0.01)

ego_f <- mutate(
    ego_f, 
    log2foldEnrichment = log2((as.numeric(sub("/\\d+", "", GeneRatio)) / as.numeric(sub("\\d+/", "", GeneRatio))) / (as.numeric(sub("/\\d+", "", BgRatio)) / as.numeric(sub("\\d+/", "", BgRatio)))))
  
ggplot(ego_f, showCategory = 10, aes(log2foldEnrichment, fct_reorder(Description, log2foldEnrichment))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
  trans = "log10",
  guide=guide_colorbar(reverse=TRUE, order=1)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_dose(12) +
  xlab("Log2(Fold-Enrichment)") +
  ylab(NULL) +
  ggtitle("BP GO enriched in PLT-driven RAM GRN") +
      scale_y_discrete(labels = function(x) str_wrap(x, 45)) +
  theme(axis.text.y = element_text(size = 10))

ggsave(file.path(opt$plots, "ntw.ora.svg"))

compendium_filt <- compendium[compendium$AGI %in% filt_gen_ms_sant$X.seq_id, ]

plt.ids <- read.delim(opt$plts, header = FALSE , col.names = c('plt.name', 'plt.id'))

int_table <- do.call(rbind, lapply(compendium_filt$AGI, function(gen){
  temp <- as.list(compendium_filt[compendium_filt$AGI == gen, c(21:26)])
  names(temp) <- sapply(names(temp), function(plt){
    return(plt.ids$plt.id[plt.ids$plt.name == plt])
    })
  plts <- na.omit(sapply(temp, function(plt){
    if(!is.na(plt)){
    if(plt != 0){
      if(plt == 1) return("A")
      else return("R")
    }
    else return(NA)
    }else return(NA)}))
  return(data.frame(
    source = names(plts),
    
    target = rep(gen, length(plts)),
    
    interaction = plts
  ))
}))

x <- unique(c(unique(int_table$source), unique(int_table$target)))

setdiff(x, prop_table$gen_id)

tfs <- read.delim(opt$tfs, comment.char = "#")

prop_table <- do.call(rbind, 
    lapply(filt_gen_ms_sant$X.seq_id, function(gen){
    return(data.frame(
        gen_id = gen,
        TF = gen%in%tfs$Gene_ID
    ))
    }
))

write.csv(int_table, file = opt$edges, row.names = F, quote = F)
write.csv(prop_table, file = opt$nodes, row.names = F, quote = F)

