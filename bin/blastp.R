suppressPackageStartupMessages(library(orthologr))  # For performing BLAST and ortholog prediction.
suppressPackageStartupMessages(library(seqinr))     # For reading and processing FASTA files.
suppressPackageStartupMessages(library(optparse))   # For command-line option parsing.

# protein_to_locus:
# This function reads a protein FASTA file and extracts the locus information
# from the annotation line of each sequence.
# Input: fasta_file - path to a protein FASTA file.
# Output: A vector of locus identifiers.
protein_to_locus <- function(fasta_file){
    # Read the FASTA file with amino acid sequences.
    ref.aa <- read.fasta(fasta_file, seqtype = 'AA')
    # For each sequence, extract the 4th element from the annotation (after splitting by space)
    # and remove the "locus=" prefix.
    x <- sapply(ref.aa, function(id_line){
        attrs <- attr(id_line, 'Annot')
        locus <- strsplit(attrs, ' ')[[1]][4]
        locus <- gsub('locus=', '', locus)
        return(locus)
    })
    return(x)
}

# parse_orthologs:
# This function replaces the numeric identifiers in the BLAST output (orthologs)
# with the corresponding locus identifiers from the query and subject FASTA files.
# It also removes duplicated entries.
# Inputs:
#   query   - path to the query protein FASTA file.
#   subject - path to the subject protein FASTA file.
#   orthologs - data frame containing BLAST results with numeric query and subject IDs.
# Output: A filtered data frame with locus names as identifiers.
parse_orthologs <- function(query, subject, orthologs){
    # Get locus names for the query and subject proteins.
    query.locus <- protein_to_locus(query)
    subject.locus <- protein_to_locus(subject)
    # Replace numeric IDs with locus names.
    orthologs$query_id <- query.locus[orthologs$query_id]
    orthologs$subject_id <- subject.locus[orthologs$subject_id]
    # Remove duplicate entries based on query and subject locus IDs.
    orthologs <- orthologs[!duplicated(orthologs$query_id), ]
    orthologs <- orthologs[!duplicated(orthologs$subject_id), ]
    return(orthologs)
}

# Define a list of command-line options.
option_list = list(
    make_option(c("-q", "--query"), type = "character", default = NULL,
              help = "Query protein FASTA file.", metavar = "path"),
    make_option(c("-s", "--subject"), type = "character", default = NULL,
              help = "Subject protein FASTA file.", metavar = "path"),
    make_option(c("-x", "--orgSubject"), type = "character", default = NULL,
              help = "Subject organism.", metavar = "name"),
    make_option(c("-z", "--orgQuery"), type = "character", default = NULL,
              help = "Query organism.", metavar = "name"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Path to save output files.", metavar = "path"),
    make_option(c("-b", "--blast"), type = "character", default = "best",
              help = "Type of BLASTp: 'best' for unidirectional best hit or 'bbh' for reciprocal best hit.", 
              metavar = "option"),
    make_option(c("-c", "--cores"), type = "integer", default = 1,
              help = "Number of cores to use in BLASTp parallel processing.", metavar = "integer"),
    make_option(c("-f", "--force"), type = "character", default = "True",
              help = "Force re-run (True/False).", metavar = "character")     
)

# Create an option parser object and parse the provided command-line arguments.
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# Print the provided query and subject FASTA file paths.
cat("Query:", opt$query, "\n")
cat("Subject:", opt$subject, "\n")


# Construct file paths for saving BLAST results and ortholog predictions.
output.file <- file.path(opt$output, paste0("orthologs.", opt$blast, ".all.tsv"))
blast.file <- file.path(opt$output, paste0("blastp.", opt$blast, ".all.tsv"))

# Initialize data frames for BLAST results and orthologs.
blast.result <- data.frame()
orthologs <- data.frame()

### Perform BLASTp or Read Existing Results ###

# If the BLAST result file does not exist or if force re-run is enabled:
if(!file.exists(blast.file) | (opt$force == "True")){

    cat("Performing BLASTp...\n")
    
    # Only perform BLAST if the query and subject files are different.
    if(opt$orgQuery != opt$orgSubject){
    
        # Depending on the selected BLAST type, run the appropriate BLAST function.
        if(opt$blast == "best"){
            cat("Selected BLASTp is unidirectional best hit.\n")
            blast.result <- blast_best(
                query_file = opt$query, 
                subject_file = opt$subject, 
                seq_type = 'protein', 
                comp_cores = opt$cores, 
                clean_folders = TRUE)
        } else if (opt$blast == "bbh"){
            cat("Selected BLASTp is bidirectional best hit.\n")
            blast.result <- blast_rec(
                query_file = opt$query, 
                subject_file = opt$subject, 
                seq_type = 'protein', 
                comp_cores = opt$cores, 
                clean_folders = TRUE)
        } else {
            stop("Select a valid BLAST option ('best' or 'bbh').")
        }
        
        # Parse the BLAST results to convert numeric IDs to locus names.
        orthologs <- parse_orthologs(opt$query, opt$subject, blast.result)
    } else {
        cat("Query and subject are the same file.\n")
    }
} else {
    # If the BLAST file exists and force re-run is not requested, bypass BLAST.
    cat("Bypassing files...\n")
    if(file.size(blast.file) > 1L){
        blast.result <- read.delim(blast.file)
        orthologs <- parse_orthologs(opt$query, opt$subject, blast.result)
    }
}

# Write the BLAST results to file.
write.table(
    blast.result, 
    file = blast.file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# Write the parsed orthologs to file.
write.table(
    orthologs, 
    file = output.file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
