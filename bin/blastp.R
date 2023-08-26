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

option_list = list(
    make_option(c("-q", "--query"), type = "character", default = NULL,
              help = "Query protein fasta file.", metavar = "path"),
    make_option(c("-s", "--subject"), type = "character", default = NULL,
              help = "Subject protein fasta file.", metavar = "path"),
    make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Path to save output files.", metavar = "path"),
    make_option(c("-b", "--blast"), type = "character", default = "best",
              help = "Type of BLASTp best (unidirectional best hit) or rec (reciprocal BBH)", 
              metavar = "option"),
    make_option(c("-c", "--cores"), type = "integer", default = 1,
              help = "Cores to use in BLASTp parallel processing.", metavar = "integer")     
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(opt$blast == "best"){
    blast.result <- blast_best(
        query_file = opt$query, 
        subject_file = opt$subject, 
        seq_type = 'protein', 
        comp_cores = opt$cores, 
        clean_folders = TRUE)
}else if (opt$blast == "bbh") {
    blast.result <- blast_rec(
        query_file = opt$query, 
        subject_file = opt$subject, 
        seq_type = 'protein', 
        comp_cores = opt$cores, 
        clean_folders = TRUE)
}else{
    stop("Select a BLAST valid option.")
}

write.table(
    blast.result, 
    file = file.path(opt$output, paste0("blastp.", opt$blast, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE)

query.locus <- protein_to_locus(opt$query)
subject.locus <- protein_to_locus(opt$subject)

blast.result$query_id <- query.locus[blast.result$query_id]
blast.result$subject_id <- subject.locus[blast.result$subject_id]

blast.result <- blast.result[!duplicated(blast.result$query_id), ]
blast.result <- blast.result[!duplicated(blast.result$subject_id), ]

write.table(
    blast.result, 
    file = file.path(opt$output, paste0("orthologs.", opt$blast, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE)
