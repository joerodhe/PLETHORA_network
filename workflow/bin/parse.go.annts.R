# run on a R session located in root work directory

# SOLANUM LYCOPERSICUM
go.annts <- read.delim("input/Solanum_lycopersicum/GENOME/ITAG4.0_goterms.txt", sep = "\t", col.names = c("gen.id", "go.annts"), header = FALSE)
go.annts$gen.id <- sapply(go.annts$gen.id, function(x) strsplit(x, "\\.")[[1]][1])
go.annts <- do.call(rbind, lapply(unique(go.annts$gen.id), function(gen){
    annts <- go.annts$go.annts[which(go.annts$gen.id == gen)]
    annts.df <- data.frame()
    if(nchar(annts) > 0){
        annts <- strsplit(annts, ",")[[1]]
        annts.df <- data.frame(
            go.id = annts,
            gen.id = rep(gen, length(annts))
        )
    }
    return(annts.df)
}))

write.table(go.annts, "input/Solanum_lycopersicum/GENOME/TERM2GENE.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

#CUCUMIS SATIVUS
go.annts <- read.delim("input/Cucumis_sativus/GENOME/Gy14_GO_anno_v1.txt", header=FALSE, col.names=c("gen.id","go.id","category","term"))
go.annts <- go.annts[which(go.annts$category == "biological_process"),]
write.table(go.annts[,c("go.id", "gen.id")], "input/Cucumis_sativus/GENOME/TERM2GENE.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


#GLYCINE MAX
go.annts <- read.delim("input/Glycine_max/GENOME/soybase_genome_annotation_v4.0_04-07-2023.txt", skip = 13)
go.annts <- go.annts[,c(1,12)]
go.annts <- do.call(rbind, lapply(unique(go.annts$Gmax_508_Wm82.a4.v1.Locus..1.), function(gen){
    annts <- go.annts$GO.Biological.Process.IDs..5.[which(go.annts$Gmax_508_Wm82.a4.v1.Locus..1. == gen)]
    gen <- strsplit(gen, "\\.")[[1]][2]
    gen <- toupper(gen)
    gen <- paste0("Glyma.", gen)
    annts.df <- data.frame()
    if(!is.na(annts)){
        annts <- strsplit(annts, " ")[[1]]
        annts.df <- data.frame(
            go.id = annts,
            gen.id = rep(gen, length(annts))
        )
    }
    return(annts.df)
}))
write.table(go.annts[,c("go.id", "gen.id")], "input/Glycine_max/GENOME/TERM2GENE.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


#ORYZA SATIVA
go.annts <- read.delim("input/Oryza_sativa_Japonica_Group/GENOME/GO-Annotations.tsv", header = FALSE)
go.annts <- go.annts[,c(1,6)]
names(go.annts) <- c("gen.id","go.id")
write.table(go.annts[,c("go.id", "gen.id")], "input/Oryza_sativa_Japonica_Group/GENOME/TERM2GENE.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
