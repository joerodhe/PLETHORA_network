dea <- readRDS("others/dea.ranksum.RDS")

mz.sc.genes <- rownames(dea[which(dea$p_val_adj < 0.05), ])

length(mz.sc.genes)

net.nodes <- read.csv("/media/external/Joel/PLTs/others/prop_table.csv")

net.nodes$SC.Meristem <- net.nodes$gen_id %in% mz.sc.genes

write.csv(
    net.nodes, 
    "output/Arabidopsis_thaliana/NETWORK/prop_table_sc.csv",
    quote = FALSE,
    row.names = FALSE, 
    na = ""
    )
