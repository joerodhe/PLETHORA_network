###############################################################################
# Description:
#   This script is designed for parsing and processing raw gene count data from
#   featureCounts, annotating the genes using a GFF3 file, generating RPKM and
#   rlog-normalized matrices for downstream analysis (e.g., PCA visualization).
#   It also outputs a metadata file for the samples, facilitating differential
#   expression workflows.
###############################################################################

###############################################################################
# Load necessary libraries (suppress startup messages for cleaner logs)
###############################################################################
suppressPackageStartupMessages(library(edgeR))   # For differential expression analysis and RPKM function
suppressPackageStartupMessages(library(DESeq2))  # For DESeq2-based normalization (rlog)
suppressPackageStartupMessages(library(optparse))# For command-line argument parsing

###############################################################################
# Function: rlog_norm
# -------------------
# Description:
#   Filters genes with at least `minCounts` counts in at least `NminSamples` 
#   samples, constructs a DESeq2 dataset, then applies rlog transformation.
#
# Params:
#   counts      : Raw read counts matrix (rows=genes, cols=samples).
#   metadata    : Data frame describing each sample (e.g., conditions).
#   minCounts   : Minimum count threshold.
#   NminSamples : Minimum number of samples meeting minCounts threshold.
#
# Returns:
#   A matrix of rlog-transformed counts with genes as rows and samples as columns.
###############################################################################
rlog_norm <- function(counts, metadata, minCounts = 1, NminSamples = 3) {
    
    # Step 1: Filter genes based on minimum count threshold in enough samples
    x <- counts[which(apply(cpm(DGEList(counts = counts)), 1, 
                            function(y) sum(y >= minCounts)) >= NminSamples), ]
    
    # Step 2: Create a DESeq2 data set from the filtered matrix and metadata
    x <- DESeqDataSetFromMatrix(countData = x, colData = metadata, design = ~dex)
    
    # Step 3: Estimate size factors for normalization
    x <- estimateSizeFactors(x)
    
    # Step 4: Apply the rlog transformation (regularized log)
    x <- rlog(x, blind = TRUE)
    
    # Step 5: Extract the transformed expression matrix
    x <- assay(x)
    
    return(x)
}

###############################################################################
# Command-line argument parsing
###############################################################################
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

###############################################################################
# Print out argument information for logging/troubleshooting
###############################################################################
cat('Reading sample IDs table:', opt$sampleIDs, '\n')
cat('Reading GFF3 file:', opt$gff, '\n')
cat('Saving metadata file:', opt$metadata, '\n')
cat('Parsing featureCounts file:', opt$featurecounts, '\n')
cat('Saving parsed gene counts file:', opt$counts, '\n')
cat('Saving RPKM converted counts matrix file:', opt$rpkm, '\n')
cat('Saving rlog converted counts matrix file:', opt$rlog, '\n\n')

###############################################################################
# Step 1: Read sample IDs (sample name to sample ID mapping)
###############################################################################
sampleIDs <- read.delim(opt$sampleIDs, header = FALSE, col.names = c('Sample', 'SampleID'))

###############################################################################
# Step 2: Read GFF3 annotation (skipping commented lines)
###############################################################################
gff <- read.delim(opt$gff, comment.char = '#', header = FALSE)

# Keep only gene features
gff <- gff[gff[, 3] == 'gene', ]

###############################################################################
# Step 3: Parse gene attributes from GFF
###############################################################################
gff.attributes <- do.call(rbind, apply(gff, 1, function(row) {
    attrs <- strsplit(row[9], ';')[[1]]  # Example: "ID=GENE;Name=GENE_NAME"
    attrs <- gsub('ID=', '', attrs)      # Remove 'ID=' prefix
    attrs <- gsub('Name=', '', attrs)    # Remove 'Name=' prefix
    
    # Special handling for ITAG4 GFF3 format (for Solanum lycopersicum)
    if (opt$gff == 'input/Solanum_lycopersicum/GENOME/ITAG4.gff3') {
        attrs[2] <- strsplit(attrs[2], '\\.')[[1]][1]  # Remove version suffix from second field
    }
    
    df <- data.frame(ID = attrs[1], Name = attrs[2])
    return(df)
}))

###############################################################################
# Step 4: Build metadata data frame for samples
# (This example assigns "dex" from the first part of SampleID and a sample
#  number repeating across conditions.)
###############################################################################
metadata <- data.frame(
    dex = sapply(sampleIDs$SampleID, function(ID) strsplit(ID, '-')[[1]][1]), # e.g. "Condition"
    sample = rep.int(c(1, 2, 3), 3),  # or some repeated pattern; customize as needed
    stringsAsFactors = TRUE
)

# Save metadata
write.table(metadata, file = opt$metadata, quote = FALSE, sep = '\t')

###############################################################################
# Step 5: Parse featureCounts output
###############################################################################
gene.count <- read.delim(opt$featurecounts, header = TRUE, comment.char = '#')

# Extract gene length vector (for RPKM)
gene.length <- gene.count$Length

# Keep only gene ID column + columns with counts (drop extra featureCounts columns)
gene.count <- gene.count[, c(1, 7:15)]

###############################################################################
# Step 6: Rename row names to GFF3 gene name, drop old Geneid column
###############################################################################
rownames(gene.count) <- sapply(gene.count$Geneid, function(id)
    gff.attributes$Name[which(gff.attributes$ID == id)]
)
gene.count <- gene.count[, -1] # Remove the 'Geneid' column

###############################################################################
# Step 7: Rename column headers to match the user-defined sample IDs
###############################################################################
colnames(gene.count) <- sapply(colnames(gene.count), function(sample) {
    # Extract final element from something like "SampleName.bam" -> "SampleName"
    sample <- tail(strsplit(sample, '\\.')[[1]], 2)[1]
    
    # Match the sample code with the sampleIDs table
    sampleID <- sampleIDs$SampleID[which(sample == sampleIDs$Sample)]
    return(sampleID)
})

cat('Dimensions of count matrix:', dim(gene.count), '\n')
# Save the raw count matrix
write.table(gene.count, file = opt$counts, quote = FALSE, sep = '\t')

###############################################################################
# Step 8: Convert raw counts to RPKM (Reads Per Kilobase of transcript per Million)
###############################################################################
cat('Converting gene counts to RPKM...\n')
gene.rpkm <- rpkm(DGEList(counts = as.matrix(gene.count)), gene.length = gene.length)
cat('Dimensions of RPKM matrix:', dim(gene.rpkm), '\n')
write.table(round(gene.rpkm, 3), file = opt$rpkm, quote = FALSE, sep = '\t')

###############################################################################
# Step 9: Apply rlog normalization for PCA/visualization
###############################################################################
cat('Converting gene counts to rlog values...\n')
gene.count.normalized <- rlog_norm(gene.count, metadata)
cat('Dimensions of rlog matrix:', dim(gene.count.normalized), '\n')
write.table(round(gene.count.normalized, 3), file = opt$rlog, quote = FALSE, sep = '\t')
