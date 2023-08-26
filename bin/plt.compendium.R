suppressPackageStartupMessages(library(readxl))

santuari <- read_excel(
  "input/Arabidopsis_thaliana/PLTs/TPC2016-00656-RAR2_Supplemental_Data_Set_1.xlsx",
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

regulated_gen <- apply(santuari[,21:26], MARGIN = 1, function(gen) 1%in%gen | -1%in%gen)
plt_gens <- santuari[regulated_gen,]

write.table(plt_gens,
  file = 'output/Arabidopsis_thaliana/PLTs/plt.compendium.tsv',
  sep = '\t', quote = FALSE)

meristem <- read.delim('output/Arabidopsis_thaliana/DEA/mz.degs.tsv')

plt_gen_mz <- meristem[meristem$gene_id%in%plt_gens$AGI,]
plt_gen_mz_comp <- sapply(meristem$gene_id, function(locus) which(locus == plt_gens$AGI))
plt_gen_mz_comp <- unlist(Filter(length,plt_gen_mz_comp))
plt_gen_mz_comp <-plt_gens[plt_gen_mz_comp,]
plt_gen_mz <- cbind(plt_gen_mz, plt_gen_mz_comp)


writeLines(plt_gen_mz$gene_id, 'output/Arabidopsis_thaliana/PLTs/mz.upregulated.degs.plt_regulated.txt')
