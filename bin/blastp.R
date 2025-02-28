###############################################################################
# Description:
#   This script performs a BLASTp search between a query and a subject protein 
#   FASTA file, detects either unidirectional best hits or reciprocal best hits, 
#   then converts numeric IDs in the BLAST results to locus names. The final 
#   orthologs are saved to TSV files for downstream usage.
###############################################################################

###############################################################################
# Load required libraries, suppressing startup messages
###############################################################################
suppressPackageStartupMessages(library(orthologr))  # For BLAST and ortholog predictions
suppressPackageStartupMessages(library(seqinr))     # For reading and processing FASTA files
suppressPackageStartupMessages(library(optparse))   # For command-line option parsing

###############################################################################
# Function: protein_to_locus
# --------------------------
# Description:
#   Reads a protein FASTA file and extracts the locus information from each
#   sequence's annotation line (the 'Annot' attribute). For each sequence,
#   it assumes the 4th space-delimited field after ">" is "locus=XXXXX".
#
# Params:
#   fasta_file: A file path to a protein FASTA.
#
# Returns:
#   A character vector of locus identifiers, in the same order as the sequences.
###############################################################################
protein_to_locus <- function(fasta_file){
    # Read the FASTA with amino acid sequences
    ref.aa <- read.fasta(fasta_file, seqtype = 'AA')
    
    # Extract the 4th field from each annotation string, stripping out 'locus='
    x <- sapply(ref.aa, function(id_line){
        attrs <- attr(id_line, 'Annot')
        locus <- strsplit(attrs, ' ')[[1]][4]
        locus <- gsub('locus=', '', locus)
        return(locus)
    })
    return(x)
}

###############################################################################
# Function: parse_orthologs
# -------------------------
# Description:
#   Takes a BLAST/ortholog output data frame and replaces numeric IDs with 
#   the corresponding locus names from the query and subject FASTA files. 
#   Removes any duplicated rows by query or subject locus ID. 
#
# Params:
#   query    : Query protein FASTA path.
#   subject  : Subject protein FASTA path.
#   orthologs: Data frame containing BLAST/ortholog results with columns 
#              'query_id' and 'subject_id' referencing numeric IDs.
#
# Returns:
#   A data frame (orthologs) with updated 'query_id' and 'subject_id' 
#   replaced by locus names, minus duplicates.
###############################################################################
parse_orthologs <- function(query, subject, orthologs){
    # Retrieve locus names from both query and subject
    query.locus <- protein_to_locus(query)
    subject.locus <- protein_to_locus(subject)
    
    # Map numeric IDs in orthologs to the actual locus strings
    orthologs$query_id <- query.locus[orthologs$query_id]
    orthologs$subject_id <- subject.locus[orthologs$subject_id]
    
    # Remove duplicates by query or subject locus
    orthologs <- orthologs[!duplicated(orthologs$query_id), ]
    orthologs <- orthologs[!duplicated(orthologs$subject_id), ]
    return(orthologs)
}

###############################################################################
# Define command-line options
###############################################################################
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
                help = "BLASTp mode: 'best' for unidirectional best hits, 'bbh' for reciprocal best hits.", 
                metavar = "option"),
    make_option(c("-c", "--cores"), type = "integer", default = 1,
                help = "Number of CPU cores for BLASTp parallelization.", metavar = "integer"),
    make_option(c("-f", "--force"), type = "character", default = "True",
                help = "Force re-run (True/False).", metavar = "character")     
)

###############################################################################
# Parse command-line arguments
###############################################################################
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

###############################################################################
# Print the provided query and subject FASTA paths for confirmation
###############################################################################
cat("Query:", opt$query, "\n")
cat("Subject:", opt$subject, "\n")

# Construct file paths for the BLAST results and ortholog output
output.file <- file.path(opt$output, paste0("orthologs.", opt$blast, ".all.tsv"))
blast.file <- file.path(opt$output, paste0("blastp.", opt$blast, ".all.tsv"))

# Initialize empty data frames
blast.result <- data.frame()
orthologs <- data.frame()

###############################################################################
# BLASTp Execution or Bypass
###############################################################################
# If the BLAST results file doesn't exist or if force re-run is "True"
if(!file.exists(blast.file) | (opt$force == "True")) {
    
    cat("Performing BLASTp...\n")
    
    # Only run BLAST if query and subject are not the same organism
    if(opt$orgQuery != opt$orgSubject){
    
        # Choose which BLAST function to call based on the user's blast option
        if(opt$blast == "best"){
            cat("Selected BLASTp mode: unidirectional best hit.\n")
            blast.result <- blast_best(
                query_file = opt$query, 
                subject_file = opt$subject, 
                seq_type = 'protein', 
                comp_cores = opt$cores, 
                clean_folders = TRUE
            )
        } else if (opt$blast == "bbh"){
            cat("Selected BLASTp mode: bidirectional best hit.\n")
            blast.result <- blast_rec(
                query_file = opt$query, 
                subject_file = opt$subject, 
                seq_type = 'protein', 
                comp_cores = opt$cores, 
                clean_folders = TRUE
            )
        } else {
            stop("Invalid BLAST option. Use 'best' or 'bbh'.")
        }
        
        # Convert numeric IDs in BLAST results to actual locus names
        orthologs <- parse_orthologs(opt$query, opt$subject, blast.result)
        
    } else {
        cat("Query and subject appear to be the same file. Skipping BLAST.\n")
    }

} else {
    # If the file is already present and no force-run, skip BLAST
    cat("Bypassing files...\n")
    # Make sure blast.file is non-empty
    if(file.size(blast.file) > 1L){
        # Read the stored BLAST results from disk
        blast.result <- read.delim(blast.file)
        # Convert numeric IDs to locus names
        orthologs <- parse_orthologs(opt$query, opt$subject, blast.result)
    }
}

###############################################################################
# Write BLAST results and ortholog pairs to output
###############################################################################
write.table(
    blast.result, 
    file = blast.file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    orthologs, 
    file = output.file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
