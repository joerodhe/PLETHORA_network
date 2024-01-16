suppressPackageStartupMessages(library(rbioapi))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(optparse))

map_panther <- function(gene.list, taxon.id) {
  mapping.df <- data.frame()

  try({
  mapping <- rba_panther_mapping(gene.list, organism = as.numeric(taxon.id))
  gen.mapped <- c(mapping$mapped_genes$gene)
  gen.unmapped <- mapping$unmapped_list$unmapped
  gene.list <-gene.list[!gene.list %in% gen.unmapped]
  #names(gen.mapped) <- gene.list

  mapping.df <- do.call(rbind, lapply(seq(gen.mapped), function(i){
    gen <- gen.mapped[i]
    gen.id <- gene.list
    annot <- gen$annotation_type_list$annotation_data_type
    annot <- do.call(rbind, lapply(annot, function(annot.data){
      if("content" %in% names(annot.data)){
      if(annot.data$content == "GO:0008150"){
        annot.data <- annot.data$annotation
        annot.data <- do.call(rbind, lapply(annot.data, function(a){
          if(!"name" %in% names(a)){
            x<-do.call(rbind, lapply(a, function(b) {
              x<-data.frame(
              gen.id = gen.id,
              name = b$name,
              go.id = b$id
              )
              }))
          }else{
            x<-data.frame(
              gen.id = gen.id,
              name = a$name,
              go.id = a$id
              )
          }
          return(x)
        }))
        return(annot.data)
      }}
    }))
    return(annot)
  }))
  })
  return(mapping.df)
}

option_list = list(
    make_option(c("-c", "--expressionMatrix"), type = "character", default = NULL,
              help = "Path to parsed featurecounts gene quantification matrix.", metavar = "character"),
    make_option(c("-t", "--taxonID"), type = "integer", default = NULL,
              help = "Taxonomic ID.", metavar = "integer"),
    make_option(c("-o", "--outputPath"), type = "character", default = NULL,
              help = "Path to save TERM2GENE and TERM2NAME files.", metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


count.matrix <- read.delim(opt$expressionMatrix)
gff.gen <- rownames(count.matrix)


start = 1
n.iter <- ceiling(length(gff.gen)/1) - start
#n.iter <- 5
go.annt <- data.frame()
for(i in seq(n.iter)){
  if((start) < length(gff.gen)) end <- start
  else end <- length(gff.gen)
  gen.map.list <- gff.gen[start:end]
  map.result <- map_panther(gen.map.list, opt$taxonID)
  go.annt <- rbind(go.annt, map.result)
  start = start+1
}

write.table(go.annt, file = file.path(opt$outputPath, "go.annotation.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

TERM2GENE <- go.annt[,c("go.id","gen.id")]
TERM2NAME <- go.annt[,c("go.id", "name")]
TERM2NAME <- TERM2NAME[!duplicated(TERM2NAME$id),]
write.table(TERM2GENE, file = file.path(opt$outputPath, "TERM2GENE.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(TERM2NAME, file = file.path(opt$outputPath, "TERM2NAME.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
