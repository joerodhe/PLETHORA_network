 suppressPackageStartupMessages(library(orthologr))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(optparse))


protein_to_locus <- function(fasta_file){
    ref.aa <- read.fasta(fasta_file, seqtype = 'AA')
    x <- sapply(ref.aa, function(id_line){
        attrs <- attr(id_line, 'Annot')
        locus <- strsplit(attrs, ' ')[[1]][4]
        locus <- gsub('locus=', '', locus)
        return(locus)
    })
    return(x)
}

parse_orthologs <- function(query, subject, orthologs){
    query.locus <- protein_to_locus(query)
    subject.locus <- protein_to_locus(subject)
    orthologs$query_id <- query.locus[orthologs$query_id]
    orthologs$subject_id <- subject.locus[orthologs$subject_id]
    orthologs <- orthologs[!duplicated(orthologs$query_id), ]
    orthologs <- orthologs[!duplicated(orthologs$subject_id), ]
    return(orthologs)
}

option_list = list(
    make_option(c("-q", "--query"), type = "character", default = NULL,
              help = "Query protein fasta file.", metavar = "path"),
    make_option(c("-s", "--subject"), type = "character", default = NULL,
              help = "Subject protein fasta file.", metavar = "path"),
    make_option(c("-x", "--orgSubject"), type = "character", default = NULL,
              help = "Subject organism.", metavar = "name"),
    make_option(c("-z", "--orgQuery"), type = "character", default = NULL,
              help = "Query organism.", metavar = "name"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Path to save output files.", metavar = "path"),
    make_option(c("-b", "--blast"), type = "character", default = "best",
              help = "Type of BLASTp best (unidirectional best hit) or rec (reciprocal BBH)", 
              metavar = "option"),
    make_option(c("-c", "--cores"), type = "integer", default = 1,
              help = "Cores to use in BLASTp parallel processing.", metavar = "integer"),
    make_option(c("-f", "--force"), type = "character", default = "True",
              help = "Force re-run", metavar = "character")     
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

cat("Query:", opt$query, "\n")
cat("Subject:", opt$subject, "\n")

output.file <- file.path(opt$output, paste0("orthologs.", opt$blast, ".all.tsv"))
blast.file <- file.path(opt$output, paste0("blastp.", opt$blast, ".all.tsv"))
blast.result <- data.frame()
orthologs <- data.frame()

if(!file.exists(blast.file) | (opt$force == "True")){

cat("Peforming BLASTp...\n")

if(opt$orgQuery != opt$orgSubject){



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

} else stop("Select a BLAST valid option.")

    orthologs <- parse_orthologs(opt$query, opt$subject, blast.result)

} else cat("Query and subject are the same file.")

} else{
    
     cat("Bypassing files...\n")
     if(file.size(blast.file) > 1L){
        blast.result <- read.delim(blast.file)
        orthologs <- parse_orthologs(opt$query, opt$subject, blast.result)
    }
}

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
