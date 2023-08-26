suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(optparse))

rlog_norm <- function(counts, metadata, minCounts = 1, NminSamples = 3){
    # this function is for normalize the raw counts for PCA visualization
    x <- counts[ which( apply( cpm(DGEList(counts = counts)), 1,
        function(y) sum(y>=minCounts)) >= NminSamples), ]
    x <- DESeqDataSetFromMatrix(countData=x, colData = metadata, design=~dex)
    x <- estimateSizeFactors(x)
    x <- rlog(x, blind=TRUE)
    x <- assay(x)
    return(x)
}


option_list = list(
    make_option(c('-s', '--sampleIDs'), type = 'character', default = NULL,
              help = 'Path to sampleIDs.tsv file.', metavar = 'character'),
    make_option(c('-f', '--featurecounts'), type = 'character', default = NULL,
              help = 'Path to featurecounts result file.', metavar = 'character'),
    make_option(c('-g', '--gff'), type = 'character', default = NULL,
              help = 'Path to gff3 annotation file.', metavar = 'character'),
    make_option(c('-m', '--metadata'), type = 'character', default = NULL,
              help = 'Path to save metadata file.', metavar = 'character'),
    make_option(c('-c', '--counts'), type = 'character', default = NULL,
              help = 'Path to save count matrix file.', metavar = 'character'),
    make_option(c('-r', '--rpkm'), type = 'character', default = NULL,
              help = 'Path to save rpkm matrix file.', metavar = 'character'),
    make_option(c('-n', '--rlog'), type = 'character', default = NULL,
              help = 'Path to rlog normalized counts matrix.', metavar = 'character')
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

cat('Reading sample IDs table:', opt$sampleIDs, '\n')
cat('Reading GFF3 file:', opt$gff, '\n')
cat('Saving metadata file:', opt$metadata, '\n')
cat('Parsing featureCounts file:', opt$featurecounts, '\n')
cat('Saving parsed gene counts file:', opt$counts, '\n')
cat('Saving RPKM converted counts matrix file:', opt$rpkm, '\n')
cat('Saving rlog converted counts matrix file:', opt$rlog, '\n\n')

sampleIDs <- read.delim(opt$sampleIDs, 
    header = FALSE, col.names = c('Sample', 'SampleID'))

gff <- read.delim(opt$gff, comment.char = '#', header = FALSE)
gff <- gff[gff[,3] == 'gene', ]
gff.attributes <- do.call(rbind, apply(gff, 1, function(row){
    attrs <- strsplit(row[9], ';')[[1]]
    attrs <- gsub('ID=', '', attrs)
    attrs <- gsub('Name=', '', attrs)
    if(opt$gff == 'input/Solanum_lycopersicum/GENOME/ITAG4.gff3'){
        attrs[2] <- strsplit(attrs[2], '\\.')[[1]][1]
    }
    df <- data.frame(
        ID = attrs[1],
        Name = attrs[2]
    )
    return(df)
    }))

#write.table(gff.attributes, file = opt$gen.ids, quote = FALSE, sep = '\t')

metadata <- data.frame(
    dex = sapply(sampleIDs$SampleID, function(ID) strsplit(ID, '-')[[1]][1]),
    sample = rep.int(c(1,2,3), 3),
    stringsAsFactors = TRUE
    )

write.table(metadata, file = opt$metadata, quote = FALSE, sep = '\t')

gene.count <- read.delim(opt$featurecounts,
    header = 1, comment.char = '#')

gene.length <- gene.count$Length

gene.count <- gene.count[,c(1, 7:15)]

rownames(gene.count) <- sapply(gene.count$Geneid, function(id)
 gff.attributes$Name[which(gff.attributes$ID == id)])
gene.count <- gene.count[,-1]

colnames(gene.count) <- sapply(colnames(gene.count), function(sample){
    sample <- tail(strsplit(sample, '\\.')[[1]], 2)[1]
    sampleID <- sampleIDs$SampleID[which(sample == sampleIDs$Sample)]
    return(sampleID)
    })

cat('Dimensions of count matrix:', dim(gene.count), '\n')

write.table(gene.count, file = opt$counts, quote = FALSE, sep = '\t')

cat('Converting gene counts to RPKM...\n')

gene.rpkm <- rpkm(DGEList(counts = as.matrix(gene.count)), gene.length = gene.length)
cat('Dimensions of RPKM matrix:', dim(gene.rpkm), '\n')
write.table(round(gene.rpkm, 3), file = opt$rpkm, quote = FALSE, sep = '\t')

cat('Converting gene counts to rlog values...\n')

gene.count.normalized <- rlog_norm(gene.count, metadata)
cat('Dimensions of rlog matrix:', dim(gene.count.normalized), '\n')
write.table(round(gene.count.normalized,3), file = opt$rlog, quote = FALSE, sep = '\t')