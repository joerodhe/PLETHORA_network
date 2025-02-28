###############################################################################
# Description:
#   This script reads a JSON file containing GO ontology information, extracts 
#   biological process (BP) terms, and then writes a two-column TSV file 
#   (go.id, term) for downstream usage, such as labeling GO terms in enrichment 
#   analyses.
###############################################################################

# Load the required package for JSON handling
library(rjson)

# Read the GO ontology JSON file
go.json <- fromJSON(file = "workflow/go.json")

# Initialize an empty data frame to store GO IDs and terms
go.df <- data.frame()

# Loop through each node in the JSON's first graph
for(node in go.json$graphs[[1]]$nodes){
  try({
    # Extract GO id, label, and type (if it's present in meta)
    g <- data.frame(
      go.id = tail(strsplit(node$id, "/")[[1]],1),
      term  = node$lbl,
      type  = node$meta$basicPropertyValues[[1]]$val
    )
    go.df <- rbind(go.df, g)
  }, silent = TRUE)
}

# Filter to only keep biological_process terms
go.df <- go.df[which(go.df$type == "biological_process"),]

# Reorder columns to keep just (go.id, term)
go.df <- go.df[, c("go.id", "term")]

# Convert underscores to colons in GO IDs (e.g., GO_0008150 -> GO:0008150)
go.df$go.id <- sapply(go.df$go.id, function(id){
  paste(strsplit(id, "_")[[1]], collapse=":")
})

# Write the resulting data frame to a tab-separated file
write.table(go.df, "workflow/TERM2NAME.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
