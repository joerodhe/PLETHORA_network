suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.At.tair.db))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(stringr))

tfbs.scan.results <- read.delim(
    'output/Arabidopsis_thaliana/PLTs/mz.upregulated.degs.plt_regulated.ft',
    comment.char = ';')

length(unique(tfbs.scan.results$X.seq_id))

plt.mz.compendium <- readLines('output/Arabidopsis_thaliana/PLTs/mz.upregulated.degs.plt_regulated.txt')
meristem <- read.delim('output/Arabidopsis_thaliana/DEA/mz.degs.tsv')
plt.compendium <- read.delim('output/Arabidopsis_thaliana/PLTs/plt.compendium.tsv')

grid.newpage()
draw.triple.venn(area1 = length(plt.compendium$AGI), 
                 area2 = length(meristem$gene_id), 
                 area3= length(unique(tfbs.scan.results$X.seq_id)), 
                 n12 = length(plt.mz.compendium), 
                 n13 = length(intersect(plt.compendium$AGI, unique(tfbs.scan.results$X.seq_id))), 
                 n23 = length(intersect(meristem$gene_id, unique(tfbs.scan.results$X.seq_id))),
                 n123 = length(
                   intersect(
                     intersect(plt.mz.compendium,unique(tfbs.scan.results$X.seq_id)),
                     meristem$gene_id)),
                 category = c("Compendium of\nPLT regulated genes", 
                              "Meristem\nupregulated DEGs", 
                              "Meristem upregulated DEGs\nwith PLT binding motif"),
                 scaled = T, 
                 fill = c("skyblue", "pink1", "mediumorchid"), 
                 lty = "blank", cat.cex = 0.9, cat.pos = 180)

c <- lapply(unique(tfbs.scan.results$X.seq_id), function(query){
  temp <- tfbs.scan.results[tfbs.scan.results$X.seq_id == query,]
  return(temp[which.min(temp$Pval),])
})
filt_gen_ms_sant <- c[[1]]
for(row in c[2:length(c)]){
  filt_gen_ms_sant <- rbind(filt_gen_ms_sant, row)
}

l2fc_MzEz <- sapply(
    filt_gen_ms_sant$X.seq_id, 
    function(id) 
        meristem$log2FoldChange_MeristemElong[meristem$gene_id == id]
        )

l2fc_MzDz <- sapply(
    filt_gen_ms_sant$X.seq_id, 
    function(id) 
        meristem$log2FoldChange_MeristemDif[meristem$gene_id == id]
        )

plt_gen_mz <- read.delim('output/Arabidopsis_thaliana/PLTs/mz.upregulated.degs.plt.compendium.tsv')

santuari_evidence <- sapply(filt_gen_ms_sant$X.seq_id, function(gen){
  regulated <- plt_gen_mz[plt_gen_mz$gene_id == gen,seq(29,34)]
  plts <- paste(na.omit(sapply(colnames(regulated), function(plt){
    if(is.na(regulated[,plt])) return(NA)
    else if(regulated[,plt] == 1) return(plt)
    else if(regulated[,plt] == -1) return(paste(plt, "(repressive)"))
    else return(NA)
  })), collapse = ";")
  return(plts)
})

filt_gen_ms_sant <- cbind(filt_gen_ms_sant, l2fc_MzDz, l2fc_MzEz, santuari_evidence)
length(rownames(filt_gen_ms_sant))

ego_f <- enrichGO(
    gene = filt_gen_ms_sant$X.seq_id, 
    keyType = 'TAIR', 
    OrgDb = org.At.tair.db, 
    ont = "BP", 
    pAdjustMethod = "BH", 
    pvalueCutoff = 0.05, 
    qvalueCutoff = 0.05)

ego_f <- mutate(
    ego_f, 
    log2foldEnrichment = log2((as.numeric(sub("/\\d+", "", GeneRatio)) / as.numeric(sub("\\d+/", "", GeneRatio))) / (as.numeric(sub("/\\d+", "", BgRatio)) / as.numeric(sub("\\d+/", "", BgRatio)))))
  
plot <- ggplot(ego_f, showCategory = 10, aes(log2foldEnrichment, fct_reorder(Description, log2foldEnrichment))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
  trans = "log10",
  guide=guide_colorbar(reverse=TRUE, order=1)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_dose(12) +
  xlab("Log2(Fold-Enrichment)") +
  ylab(NULL) +
  ggtitle("GO enriched for meristem up-regulated DEGs\nwith PLT binding site") +
      scale_y_discrete(labels = function(x) str_wrap(x, 45)) +
  theme(axis.text.y = element_text(size = 10))

plt.compendium_filt <- plt.compendium[plt.compendium$AGI %in% filt_gen_ms_sant$X.seq_id, ]

plt.ids <- read.delim('others/plt.ids.txt', header = FALSE , col.names = c('plt.name', 'plt.id'))

int_table <- do.call(rbind, lapply(plt.compendium_filt$AGI, function(gen){
  temp <- as.list(plt.compendium_filt[plt.compendium_filt$AGI == gen, c(21:26)])
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

write.csv(int_table, file = "output/Arabidopsis_thaliana/NETWORK/network_table.csv", row.names = F, quote = F)

tfs <- read.delim("input/Arabidopsis_thaliana/GENOME/TAIR10.TF_list.txt", comment.char = "#")

prop_table <- do.call(rbind, 
    lapply(filt_gen_ms_sant$X.seq_id, function(gen){
    return(data.frame(
        gen_id = gen,
        TF = gen%in%tfs$Gene_ID
    ))
    }
))

write.csv(prop_table, file = "output/Arabidopsis_thaliana/NETWORK/prop_table.csv", row.names = F, quote = F)
