suppressPackageStartupMessages(library(edgeR))  # For differential expression analysis
suppressPackageStartupMessages(library(DESeq2)) # For normalization and transformation
suppressPackageStartupMessages(library(optparse)) # For command-line argument parsing

# Function to normalize raw counts using rlog transformation for PCA visualization
rlog_norm <- function(counts, metadata, minCounts = 1, NminSamples = 3) {
    # Filter genes with at least 'minCounts' in 'NminSamples' samples
    x <- counts[which(apply(cpm(DGEList(counts = counts)), 1, 
                            function(y) sum(y >= minCounts)) >= NminSamples), ]
    
    # Create DESeq2 dataset with expression data and metadata
    x <- DESeqDataSetFromMatrix(countData = x, colData = metadata, design = ~dex)
    
    # Estimate size factors for normalization
    x <- estimateSizeFactors(x)
    
    # Apply rlog transformation (regularized log transformation)
    x <- rlog(x, blind = TRUE)
    
    # Extract transformed expression matrix
    x <- assay(x)
    
    return(x)
}

option_list = list(
    make_option(c('-s', '--sampleIDs'), type = 'character', default = NULL,
                help = 'Path to sampleIDs.tsv file.', metavar = 'character'),
    make_option(c('-f', '--featurecounts'), type = 'character', default = NULL,
                help = 'Path to featureCounts result file.', metavar = 'character'),
    make_option(c('-g', '--gff'), type = 'character', default = NULL,
                help = 'Path to GFF3 annotation file.', metavar = 'character'),
    make_option(c('-m', '--metadata'), type = 'character', default = NULL,
                help = 'Path to save metadata file.', metavar = 'character'),
    make_option(c('-c', '--counts'), type = 'character', default = NULL,
                help = 'Path to save count matrix file.', metavar = 'character'),
    make_option(c('-r', '--rpkm'), type = 'character', default = NULL,
                help = 'Path to save RPKM matrix file.', metavar = 'character'),
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

# Read sample IDs from file
sampleIDs <- read.delim(opt$sampleIDs, header = FALSE, col.names = c('Sample', 'SampleID'))

# Read GFF3 annotation file (skipping commented lines)
gff <- read.delim(opt$gff, comment.char = '#', header = FALSE)

# Extract only gene features from GFF3 file
gff <- gff[gff[, 3] == 'gene', ]

# Parse gene attributes from GFF file
gff.attributes <- do.call(rbind, apply(gff, 1, function(row) {
    attrs <- strsplit(row[9], ';')[[1]]  # Split attributes field
    attrs <- gsub('ID=', '', attrs)  # Remove 'ID=' prefix
    attrs <- gsub('Name=', '', attrs)  # Remove 'Name=' prefix
    
    # Special handling for ITAG4 GFF3 format (for Solanum lycopersicum)
    if (opt$gff == 'input/Solanum_lycopersicum/GENOME/ITAG4.gff3') {
        attrs[2] <- strsplit(attrs[2], '\\.')[[1]][1]  # Remove versioning suffix
    }
    
    df <- data.frame(ID = attrs[1], Name = attrs[2])
    return(df)
}))

# Create metadata dataframe for samples
metadata <- data.frame(
    dex = sapply(sampleIDs$SampleID, function(ID) strsplit(ID, '-')[[1]][1]), # Extract condition
    sample = rep.int(c(1, 2, 3), 3), # Assign sample numbers
    stringsAsFactors = TRUE
)

# Save metadata file
write.table(metadata, file = opt$metadata, quote = FALSE, sep = '\t')

# Read featureCounts gene count matrix
gene.count <- read.delim(opt$featurecounts, header = TRUE, comment.char = '#')

# Extract gene length for RPKM calculation
gene.length <- gene.count$Length

# Keep only relevant columns: Gene ID and sample columns
gene.count <- gene.count[, c(1, 7:15)]

# Rename row names to gene names from GFF3 attributes
rownames(gene.count) <- sapply(gene.count$Geneid, function(id)
    gff.attributes$Name[which(gff.attributes$ID == id)])
gene.count <- gene.count[, -1] # Remove Gene ID column

# Rename column names to sample IDs from metadata
colnames(gene.count) <- sapply(colnames(gene.count), function(sample) {
    sample <- tail(strsplit(sample, '\\.')[[1]], 2)[1]  # Extract last part of sample name
    sampleID <- sampleIDs$SampleID[which(sample == sampleIDs$Sample)]
    return(sampleID)
})

cat('Dimensions of count matrix:', dim(gene.count), '\n')

write.table(gene.count, file = opt$counts, quote = FALSE, sep = '\t')

# Convert raw counts to RPKM (Reads Per Kilobase per Million mapped reads)
cat('Converting gene counts to RPKM...\n')
gene.rpkm <- rpkm(DGEList(counts = as.matrix(gene.count)), gene.length = gene.length)
cat('Dimensions of RPKM matrix:', dim(gene.rpkm), '\n')
write.table(round(gene.rpkm, 3), file = opt$rpkm, quote = FALSE, sep = '\t')

# Convert raw counts to rlog-normalized values for visualization
cat('Converting gene counts to rlog values...\n')
gene.count.normalized <- rlog_norm(gene.count, metadata)
cat('Dimensions of rlog matrix:', dim(gene.count.normalized), '\n')
write.table(round(gene.count.normalized, 3), file = opt$rlog, quote = FALSE, sep = '\t')
