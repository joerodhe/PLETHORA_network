###############################################################################
# Load a saved differential expression analysis (DEA) results object
###############################################################################
dea <- readRDS("others/dea.ranksum.RDS")

###############################################################################
# Identify genes that are significantly differentially expressed (adjusted p<0.05)
###############################################################################
mz.sc.genes <- rownames(dea[which(dea$p_val_adj < 0.05), ])

###############################################################################
# Print the number of significant genes
###############################################################################
length(mz.sc.genes)

###############################################################################
# Read a CSV file containing reference network nodes with metadata
###############################################################################
net.nodes <- read.csv("/media/external/Joel/PLTs/others/prop_table.csv")

###############################################################################
# Add a new logical column (SC.Meristem) indicating which genes are in the DE list
###############################################################################
net.nodes$SC.Meristem <- net.nodes$gen_id %in% mz.sc.genes

###############################################################################
# Write the updated table to a new CSV file without quotes or row names
###############################################################################
write.csv(
    net.nodes, 
    "output/Arabidopsis_thaliana/NETWORK/prop_table_sc.csv",
    quote = FALSE,
    row.names = FALSE, 
    na = ""
)
