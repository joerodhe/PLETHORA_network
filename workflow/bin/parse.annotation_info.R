suppressPackageStartupMessages(library(optparse))

option_list = list(
    make_option(c("-a", "--annotation"), type = "character", default = NULL,
              help = "annotation_info from JGI", metavar = "path"),
    make_option(c("-s", "--go"), type = "character", default = NULL,
              help = "TERM2GENE output file", metavar = "path")    
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


#opt$annotation <- 'input/Zea_mays/GENOME/Zm-B73-REFERENCE-NAM-5.0.55.annotation_info.txt'
#opt$go <- 'input/Zea_mays/GENOME/TERM2GENE.tsv'

annotation <- read.delim(opt$annotation)

annotation.go <- do.call(
    rbind,
    lapply(
        unique(annotation[,2]), function(gen){
            go.ids <- annotation[which(annotation[,2] == gen),10]
            go.ids <- go.ids[go.ids != '']
            if(length(go.ids)){
            go.ids <- sapply(go.ids, function(go.list) strsplit(go.list, '[, ]+')[[1]])
            go.ids <- unlist(as.list(go.ids))
            x <- data.frame(
                go.id = unique(go.ids),
                gen.id = gen
            )
            return(x)
            }
        })
)

write.table(annotation.go, opt$go, sep = '\t', quote = FALSE, row.names = FALSE)
