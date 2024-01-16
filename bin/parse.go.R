library(rjson)

go.json <- fromJSON(file = "workflow/go.json")


go.df <- data.frame()

    for(node in go.json$graphs[[1]]$nodes){
        try({
        g <- data.frame(
            go.id = tail(strsplit(node$id, "/")[[1]],1),
            term = node$lbl,
            type = node$meta$basicPropertyValues[[1]]$val
        )
        go.df <- rbind(go.df, g)
        }, silent = TRUE)

    }

go.df <- go.df[which(go.df$type == "biological_process"),]
go.df <- go.df[,c("go.id", "term")]
go.df$go.id <- sapply(go.df$go.id, function(id){
    paste(strsplit(id, "_")[[1]], collapse=":")
})
write.table(go.df, "workflow/TERM2NAME.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
