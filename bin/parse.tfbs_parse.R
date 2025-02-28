###############################################################################
# Description:
#   This script focuses on identifying meristem up-regulated DEGs (differentially 
#   expressed genes) that contain a PLT (PLETHORA) binding site motif, integrating 
#   multiple data sources: 
#   1) A scanning result of TFBS (Transcription Factor Binding Site) matches, 
#   2) A list of meristem up-regulated DEGs, 
#   3) A compendium table of PLT-regulated genes, and 
#   4) A reference table of TFs.
#   It then:
#   - Creates a triple Venn diagram showing overlaps among these gene sets.
#   - Extracts the "best" (lowest p-value) TFBS match per gene.
#   - Merges log2 fold-change values from DEA results.
#   - Performs GO enrichment analysis on the final list of genes.
#   - Exports a small gene regulatory network table (interactions) as well as
#     a property table (whether each gene is a TF) for further use.
###############################################################################

###############################################################################
# Load libraries necessary for the analysis
###############################################################################
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.At.tair.db))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(DOSE))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(stringr))

###############################################################################
# Step 1: Read input data
###############################################################################

# 1.1: TFBS scanning results, which list genes (seq_id) that have PLT TF binding
tfbs.scan.results <- read.delim(
    'output/Arabidopsis_thaliana/PLTs/mz.upregulated.degs.plt_regulated.ft',
    comment.char = ';'
)
cat("Number of unique seq_ids in TFBS scan:", length(unique(tfbs.scan.results$X.seq_id)), "\n")

# 1.2: Genes identified as "meristem upregulated DEGs with PLT-regulation"
plt.mz.compendium <- readLines('output/Arabidopsis_thaliana/PLTs/mz.upregulated.degs.plt_regulated.txt')

# 1.3: The table of meristem up-regulated DEGs, including log2FC columns
meristem <- read.delim('output/Arabidopsis_thaliana/DEA/mz.degs.tsv')

# 1.4: A compendium table describing which genes are regulated by which PLTs
plt.compendium <- read.delim('output/Arabidopsis_thaliana/PLTs/plt.compendium.tsv')

###############################################################################
# Step 2: Create a triple Venn diagram showing overlaps among:
#         (1) PLT-regulated compendium, (2) meristem upregulated DEGs,
#         (3) Genes from TFBS scanning results (with PLT binding motifs).
###############################################################################
grid.newpage()
draw.triple.venn(
  area1 = length(plt.compendium$AGI), 
  area2 = length(meristem$gene_id), 
  area3 = length(unique(tfbs.scan.results$X.seq_id)), 
  n12 = length(plt.mz.compendium), 
  n13 = length(intersect(plt.compendium$AGI, unique(tfbs.scan.results$X.seq_id))), 
  n23 = length(intersect(meristem$gene_id, unique(tfbs.scan.results$X.seq_id))),
  n123 = length(
    intersect(
      intersect(plt.mz.compendium, unique(tfbs.scan.results$X.seq_id)),
      meristem$gene_id
    )
  ),
  category = c(
    "Compendium of\nPLT regulated genes", 
    "Meristem\nupregulated DEGs", 
    "Meristem upregulated DEGs\nwith PLT binding motif"
  ),
  scaled = TRUE, 
  fill = c("skyblue", "pink1", "mediumorchid"), 
  lty = "blank", cat.cex = 0.9, cat.pos = 180
)

###############################################################################
# Step 3: Extract the best TFBS match (lowest P-value) per gene (seq_id)
###############################################################################
c <- lapply(unique(tfbs.scan.results$X.seq_id), function(query){
  temp <- tfbs.scan.results[tfbs.scan.results$X.seq_id == query,]
  return(temp[which.min(temp$Pval),])
})

# Combine these best matches into one data frame
filt_gen_ms_sant <- c[[1]]
for(row in c[2:length(c)]){
  filt_gen_ms_sant <- rbind(filt_gen_ms_sant, row)
}

###############################################################################
# Step 4: Merge in log2 fold-change info for MZ-EZ and MZ-DZ comparisons
###############################################################################
l2fc_MzEz <- sapply(
    filt_gen_ms_sant$X.seq_id, 
    function(id) meristem$log2FoldChange_MeristemElong[meristem$gene_id == id]
)
l2fc_MzDz <- sapply(
    filt_gen_ms_sant$X.seq_id, 
    function(id) meristem$log2FoldChange_MeristemDif[meristem$gene_id == id]
)

###############################################################################
# Step 5: Parse the santuari_evidence, i.e., which PLTs regulate each gene
###############################################################################
plt_gen_mz <- read.delim('output/Arabidopsis_thaliana/PLTs/mz.upregulated.degs.plt.compendium.tsv')

santuari_evidence <- sapply(filt_gen_ms_sant$X.seq_id, function(gen){
  # Subset rows for this gene
  regulated <- plt_gen_mz[plt_gen_mz$gene_id == gen, seq(29,34)]
  
  # Gather which PLTs (columns) are either 1 or -1
  plts <- paste(
    na.omit(sapply(colnames(regulated), function(plt){
      if(is.na(regulated[, plt])) return(NA)
      else if(regulated[, plt] == 1)  return(plt)  # Activating
      else if(regulated[, plt] == -1) return(paste(plt, "(repressive)"))
      else return(NA)
    })), 
    collapse = ";"
  )
  return(plts)
})

# Append new columns to the best TFBS match data frame
filt_gen_ms_sant <- cbind(filt_gen_ms_sant, l2fc_MzDz, l2fc_MzEz, santuari_evidence)
cat("Number of genes in final TFBS/DEA-l2FC combined data:", nrow(filt_gen_ms_sant), "\n")

###############################################################################
# Step 6: Perform GO enrichment on the final gene set
###############################################################################
ego_f <- enrichGO(
    gene = filt_gen_ms_sant$X.seq_id, 
    keyType = 'TAIR', 
    OrgDb = org.At.tair.db, 
    ont = "BP", 
    pAdjustMethod = "BH", 
    pvalueCutoff = 0.05, 
    qvalueCutoff = 0.05
)

# Calculate log2(Fold-Enrichment) for each term
ego_f <- mutate(
    ego_f, 
    log2foldEnrichment = log2(
      (as.numeric(sub("/\\d+", "", GeneRatio)) / as.numeric(sub("\\d+/", "", GeneRatio))) /
      (as.numeric(sub("/\\d+", "", BgRatio)) / as.numeric(sub("\\d+/", "", BgRatio)))
    )
)

# Create a dot plot showing the top 10 enriched GO terms
plot <- ggplot(ego_f, showCategory = 10, aes(log2foldEnrichment, fct_reorder(Description, log2foldEnrichment))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_gradientn(
    colours=c("#f7ca64", "#46bac2", "#7e62a3"),
    trans = "log10",
    guide=guide_colorbar(reverse=TRUE, order=1)
  ) +
  scale_size_continuous(range=c(2, 10)) +
  theme_dose(12) +
  xlab("Log2(Fold-Enrichment)") +
  ylab(NULL) +
  ggtitle("GO enriched for meristem up-regulated DEGs\nwith PLT binding site") +
  scale_y_discrete(labels = function(x) str_wrap(x, 45)) +
  theme(axis.text.y = element_text(size = 10))

###############################################################################
# Step 7: Export the network of PLT -> Gene interactions
###############################################################################
plt.compendium_filt <- plt.compendium[plt.compendium$AGI %in% filt_gen_ms_sant$X.seq_id, ]

plt.ids <- read.delim('others/plt.ids.txt', header = FALSE, col.names = c('plt.name', 'plt.id'))

# For each gene, find which PLTs regulate it (activation or repression)
int_table <- do.call(rbind, lapply(plt.compendium_filt$AGI, function(gen){
  temp <- as.list(plt.compendium_filt[plt.compendium_filt$AGI == gen, c(21:26)])
  # Replace column names by the PLT gene IDs from plt.ids
  names(temp) <- sapply(names(temp), function(plt){
    return(plt.ids$plt.id[plt.ids$plt.name == plt])
  })
  
  # Keep only non-zero PLT relationships
  plts <- na.omit(sapply(temp, function(plt_val){
    if(!is.na(plt_val)){
      if(plt_val != 0){
        if(plt_val == 1)  return("A") # Activating
        else              return("R") # Repressing
      }
      else return(NA)
    }else return(NA)
  }))
  
  return(data.frame(
    source      = names(plts),
    target      = rep(gen, length(plts)),
    interaction = plts
  ))
}))

write.csv(
  int_table, 
  file = "output/Arabidopsis_thaliana/NETWORK/network_table.csv", 
  row.names = FALSE, 
  quote = FALSE
)

###############################################################################
# Step 8: Build a property table marking which of these genes are TFs
###############################################################################
tfs <- read.delim("input/Arabidopsis_thaliana/GENOME/TAIR10.TF_list.txt", comment.char="#")

prop_table <- do.call(rbind, lapply(filt_gen_ms_sant$X.seq_id, function(gen){
  data.frame(
    gen_id = gen,
    TF     = gen %in% tfs$Gene_ID
  )
}))

write.csv(
  prop_table, 
  file = "output/Arabidopsis_thaliana/NETWORK/prop_table.csv", 
  row.names = FALSE, 
  quote = FALSE
)
