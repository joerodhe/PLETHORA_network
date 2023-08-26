suppressPackageStartupMessages(library(WGCNA))
options(stringsAsFactors = FALSE)
enableWGCNAThreads(60)

labeledBarplot2 <- function(x, names, g = NULL, horiz = FALSE,
                           srt = if (horiz) 0 else 45,
                           adj = if (horiz) c(1, 0.5) else c(1, 0.5),
                           namesColor = 1, cex.names = 1,
                           addGuide = TRUE, guideColor = "grey30", guideLty = 2, ...){
  
  if (!is.null(g)) {
    mp = verboseBarplot(x, g, names.arg = NA, horiz = horiz, ...);
    heights = attr(mp, "height");
  } else {
    mp = barplot(x, names.arg = NA, horiz= horiz, ...);
    heights = x;
    attr(mp, "height") = x;
    attr(mp, "stdErr") = rep(0, length(x));
  }
  
  box = par("usr");
  
  if (horiz)
  {
    yText = mp;
    xText = rep(box[1] - 0.01 * (box[2] - box[1]), length(mp));
  } else {
    xText = mp;
    yText = rep(box[3] - 0.02 * (box[4] - box[3]), length(mp));
  }
  
  text(xText, yText, names, srt = srt, adj = adj, col = namesColor, cex = cex.names, xpd = TRUE);
  if (addGuide) for (i in 1:length(mp)) # nolint
    if (horiz)
    {
      lines(c(min(0, heights[i]), box[1] - 0.02 * (box[2] - box[1])),  rep(mp[i], 2), 
            col = guideColor, lty = guideLty);
    } else {
      lines(rep(mp[i], 2), c(min(0, heights[i]), box[3] - 0.02 * (box[4] - box[3])), 
            col = guideColor, lty = guideLty);
    }
  
  invisible(mp);
  
}

# Sets that we work with
expression.matrices <- c(
    "output/Arabidopsis_thaliana/FEATURECOUNTS/TAIR10.gene.rlog.tsv",
    "output/Cucumis_sativus/FEATURECOUNTS/Gy14.gene.rlog.tsv",
    "output/Glycine_max/FEATURECOUNTS/Wm82.gene.rlog.tsv",
    "output/Oryza_sativa_Japonica_Group/FEATURECOUNTS/MSUv7.gene.rlog.tsv",
    "output/Solanum_lycopersicum/FEATURECOUNTS/ITAG4.gene.rlog.tsv",
    "output/Zea_mays/FEATURECOUNTS/Zm-B73-REFERENCE-NAM-5.0.55.gene.rlog.tsv"
)

organisms <- c(
    "Arabidopsis_thaliana",
    "Cucumis_sativus",
    "Glycine_max",
    "Oryza_sativa_Japonica_Group",
    "Solanum_lycopersicum",
    "Zea_mays"
)

metadata <- c(
    "temp/Arabidopsis_thaliana/metadata.tsv",
    "temp/Cucumis_sativus/metadata.tsv",
    "temp/Glycine_max/metadata.tsv",
    "temp/Oryza_sativa_Japonica_Group/metadata.tsv",
    "temp/Solanum_lycopersicum/metadata.tsv",
    "temp/Zea_mays/metadata.tsv"
)

samples <- c("MZ","DZ","EZ")

nOrgs <- length(organisms)
nSamples <- length(samples)
nSets <- length(organisms)*nSamples

organism.ref <- "Arabidopsis_thaliana"
orthologs <- read.delim("output/interspecie/orthologs.tsv")

parse_nodes <- function(organism, orthologs, TFBS = FALSE){

    wgcna.path <- file.path("output", organism, "WGCNA")

    module.colors <- list.files(wgcna.path, pattern = ".ft")
    module.colors <- gsub(x = module.colors, pattern = ".nodes.ft", replacement = "")

    genes <- do.call(
        rbind,
        lapply(
            module.colors,
            function(color){
                read.delim(file.path(wgcna.path, paste0(color,".nodes.tsv")))
            }
        ))


    orthologs <- orthologs[orthologs[,organism] %in% genes$nodeName, ]

    gen.plt.regulated <- sapply(
        module.colors,
        function(color){
            ft.file <- file.path(wgcna.path, paste0(color,".nodes.ft"))
            if(file.size(ft.file) > 1L){
                plt.tfbs <- read.delim(ft.file, comment.char = ";")
                return(unique(plt.tfbs$X.seq_id))
            }
        })
    gen.plt.regulated <- unlist(gen.plt.regulated)

    genes$PLT.TFBS <- genes$nodeName %in% gen.plt.regulated

    if(TFBS)
        mz.genes.plt.regulated <- genes$nodeName[which(genes$MZ.expressed & genes$PLT.TFBS)]
    else
        mz.genes.plt.regulated <- genes$nodeName[which(genes$MZ.expressed)]

    mz.nodes.plt.regulated <- orthologs[orthologs[,organism] %in% mz.genes.plt.regulated, 1]

    return(mz.nodes.plt.regulated)
}

mz.nodes.plt.regulated <- list()
idx <- 1
for(organism in organisms){
    mz.nodes.plt.regulated[[idx]] <- parse_nodes(organism, orthologs, FALSE)
    idx <- idx + 1
}

length(mz.nodes.plt.regulated)
mz.nodes.plt.regulated <- Reduce(intersect, mz.nodes.plt.regulated)
length(mz.nodes.plt.regulated)

mz.nodes.plt.regulated.ara <- parse_nodes(organism.ref, orthologs, TFBS = TRUE)
length(mz.nodes.plt.regulated.ara)

mz.nodes.plt.regulated <- intersect(mz.nodes.plt.regulated, mz.nodes.plt.regulated.ara)
length(mz.nodes.plt.regulated)

mz.nodes.plt.regulated.module <- setNames(
    object = sample(c(0,1), nrow(orthologs), replace = TRUE),
    orthologs[,1]
)
mz.nodes.plt.regulated.module[mz.nodes.plt.regulated] <- 2
table(mz.nodes.plt.regulated.module)

ref.net <- read.csv("output/Arabidopsis_thaliana/NETWORK/prop_table.csv")
ref.net.genes <- ref.net$gen_id[which(ref.net$gen_id %in% orthologs[,organism.ref])]
ref.net.nodes <- orthologs[orthologs[,organism.ref] %in% ref.net.genes,1]
length(ref.net.nodes)
mz.nodes.plt.regulated.module[ref.net.nodes] <- 2
mz.nodes.plt.regulated <- intersect(mz.nodes.plt.regulated, ref.net.nodes)
length(mz.nodes.plt.regulated)
mz.nodes.plt.regulated.module[mz.nodes.plt.regulated] <- 2
table(mz.nodes.plt.regulated.module)
length(mz.nodes.plt.regulated.module)

# Object that will contain the expression data
multiExpr <- list()
moduleList <- list()
sampleLabels = c()
idx <- 1
for(i in 1:nOrgs){
    expression.matrix <- read.delim(expression.matrices[i], check.names=FALSE)
    #expression.matrix <- log2(expression.matrix + 0.1)
    #if(!all(orthologs[,organisms[i]] %in% rownames(expression.matrix))){
    #    print(head(orthologs[!orthologs[,organisms[i]] %in% rownames(expression.matrix), ]))
    #    }
    expression.matrix <- expression.matrix[orthologs[,organisms[i]], ]
    expression.matrix[is.na(expression.matrix)] <- 0
    expression.matrix <- t(expression.matrix)
    colnames(expression.matrix) <- orthologs[,1]
    sample.metadata <- read.delim(metadata[i])
    for(z in 1:nSamples){
        sample <- samples[z]
        z.samples <- row.names(sample.metadata[sample.metadata$dex == sample,])
        sample.expression <- expression.matrix[z.samples, ]
        multiExpr[[idx]] <- list(data = sample.expression)
        organism <- strsplit(organisms[i], "_")[[1]][1]
        sampleLabel <- paste(sample, organism, sep = "_")
        sampleLabels <- c(sampleLabels, sampleLabel)
        moduleList[[idx]] <- mz.nodes.plt.regulated.module
        idx = idx+1
    }
}

length(multiExpr)

#sampleLabels <- sapply(sampleLabels, function(x) strsplit(x,"_")[[1]][1])

names(multiExpr) <- sampleLabels
names(moduleList) <- sampleLabels

lapply(multiExpr, lapply, dim)
lapply(moduleList, length)
lapply(moduleList, table)

ggs = goodSamplesGenesMS(multiExpr, minNSamples = 1, minNGenes = 1)
for (set in 1:nSets){
    multiExpr[[set]]$data = multiExpr[[set]]$data[, ggs$goodGenes]
    moduleList[[set]] = moduleList[[set]][ggs$goodGenes]
}

mp <- modulePreservation(
    multiExpr, 
    moduleList,
    referenceNetworks = 1,
    loadPermutedStatistics = FALSE,
    corFnc = "cor",
    nPermutations = 200,
    verbose = 4,
    maxGoldModuleSize = 1000, 
    maxModuleSize = 1000,
    parallelCalculation = TRUE, 
    checkData = FALSE)

multiMEs = multiSetMEs(multiExpr, universalColors = moduleList[[1]])

mz.nodes.plt.regulated <- mz.nodes.plt.regulated[mz.nodes.plt.regulated %in% names(moduleList[[set]])]
ref.net.nodes <- ref.net.nodes[ref.net.nodes %in% names(moduleList[[set]])]
nGenes <- length(ref.net.nodes)

KMEpathway = matrix(0, nGenes, nSets)
for (set in 1:nSets){
    KMEpathway[, set] <- cor(
        multiExpr[[set]]$data[,ref.net.nodes], 
        multiMEs[[set]]$data[, 3], 
        use = "p")
}

row.names(KMEpathway) <- ref.net.nodes
colnames(KMEpathway) <- names(multiExpr)

ref <- 1
preservation.results <- mp$preservation$Z[[1]][-ref]
names(preservation.results) <- sub(
    names(preservation.results),
    pattern = "inColumnsAlsoPresentIn.",
    replacement = "",
    fixed = TRUE
)

Zs <- do.call(
    rbind,
    lapply(
        preservation.results, function(result)
            data.frame(
                Zsummary = result['0',]$Zsummary.pres,
                Zdensity = result['0',]$Zdensity.pres,
                Zconn = result['0',]$Zconnectivity.pres
            )
    ))

Zs <- Zs[order(row.names(Zs), decreasing = TRUE),]

order = c(1:18);

order2 = c(1,3,2,6,4,7,10,);
order2 <- c(1,3,2,5,4,7,6,9,10,8,11,12,15,14,13,17,16,18)

sizeGrWindow(10, 8)
pdf(file = "indivPathway-allgenes-GOCholesterolBiosynthesis-summaryForPaper.pdf", w = 10, h=18)
layout(matrix(c(1:21), 7, 3, byrow = TRUE));
#par(mar = c(5.3, 3.3, 4, 0.5));
#par(mgp = c(1.9, 0.6, 0))
# Plot of summary Z statistics
for (i in 1:3){
    pres <- Zs[,i]
    labeledBarplot2(pres, names = row.names(Zs),
    main = spaste(LETTERS[i], ". ", colnames(Zs)[i], "\nReference: MZ Arabidopsis"),
    cex.names = 1.2,
    cex.axis = 1.2, ylab = colnames(Zs)[i], cex.main = 1.4, cex.lab = 1.2)
    abline(h = 10, col = "red", lty = 2)
}
# KME scatterplots
par(mar = c(3.3, 3.3, 4, 0.5));
par(mgp = c(1.9, 0.6, 0))
ref = 1
ind = 5
for (set in 1:nSets){
if (set!=ref){
verboseScatterplot(KMEpathway[, ref], KMEpathway[, set],
xlab = spaste("KME ", colnames(KMEpathway)[set]),
ylab = "KME MZ_Arabidopsis",  abline = TRUE)
#main = spaste(LETTERS[ind], ". KME in ", organisms[order2[set]], "\nvs. ",
#organisms[ref], "\n"),
#cex.main = 1.4, cex.lab = 1.2, cex.axis = 1.2
ind = ind + 1;
} else {
plot(c(0,1), type ="n", axes = FALSE, xlab = "", ylab ="")
}
}
# If plotting into a file, close it
dev.off()


## log2(rlog + 0.1) All mz expressed and all PLT TFBS (42 genes)
#                Zsummary   Zdensity     Zconn
#DZ_Arabidopsis  8.623043  5.1459132 12.100173
#EZ_Arabidopsis 15.712200 12.5364511 18.887949
#MZ_Cucumis     11.249632 10.0757283 12.423536
#DZ_Cucumis     10.435592  9.3429399 11.528245
#EZ_Cucumis     12.167795 10.4930727 13.842517
#MZ_Glycine      8.203014  6.7627370  9.643292
#DZ_Glycine      3.109191  1.8908642  4.327519
#EZ_Glycine      7.642925  6.5778551  8.707995
#MZ_Oryza        7.444263  6.5707264  8.317800
#DZ_Oryza        4.764301  2.5090091  7.019594
#EZ_Oryza        3.799626  2.4228707  5.176381
#MZ_Solanum     19.433505 15.7646656 23.102344
#DZ_Solanum     10.891124  9.4468852 12.335363
#EZ_Solanum      4.135676  2.4922416  5.779110
#MZ_Zea          2.192255  0.8991507  3.485359
#DZ_Zea          4.190139  3.2928087  5.087470
#EZ_Zea          5.017914  3.7953345  6.240493

# log(rlog + 0.1) All mz expressed TFBS in A thaliana (248 genes)
#                Zsummary   Zdensity     Zconn
#DZ_Arabidopsis  5.976255  3.3454466  8.607064
#EZ_Arabidopsis 13.919538 11.3380683 16.501008
#MZ_Cucumis      9.917341  8.3661675 11.468514
#DZ_Cucumis      5.060283  2.6801194  7.440447
#EZ_Cucumis      9.008231  7.7993891 10.217073
#MZ_Glycine      5.908547  3.5108668  8.306227
#DZ_Glycine      4.243834  3.2196507  5.268016
#EZ_Glycine      4.667122  3.5036650  5.830579
#MZ_Oryza        7.844304  6.8117208  8.876888
#DZ_Oryza        4.313809  2.9151574  5.712460
#EZ_Oryza        4.091242  2.6437304  5.538754
#MZ_Solanum     16.712589 14.7390001 18.686178
#DZ_Solanum      9.081776  7.8830167 10.280536
#EZ_Solanum      6.939055  5.9864880  7.891622
#MZ_Zea          2.195394  0.9364931  3.454294
#DZ_Zea          4.235797  2.7686975  5.702897
#EZ_Zea          3.749348  2.2238026  5.274894

# log(rlog + 0.1) All mz expressed Manual net as ref net (90 genes)
#                Zsummary   Zdensity     Zconn
#DZ_Arabidopsis  9.410903  6.7242692 12.097537
#EZ_Arabidopsis 18.610396 14.9210954 22.299697
#MZ_Cucumis      9.662665  6.9866600 12.338670
#DZ_Cucumis     10.927949 10.0126763 11.843223
#EZ_Cucumis     14.361193 12.8186509 15.903736
#MZ_Glycine      8.238428  6.8499378  9.626918
#DZ_Glycine      5.200005  4.0844553  6.315555
#EZ_Glycine      4.926837  3.5728254  6.280848
#MZ_Oryza        9.434156  7.5385278 11.329784
#DZ_Oryza        7.348191  6.1404646  8.555918
#EZ_Oryza        1.448590  0.4728998  2.424280
#MZ_Solanum     19.815989 16.7264863 22.905492
#DZ_Solanum     10.176244  9.4420001 10.910488
#EZ_Solanum      5.880304  4.3020399  7.458568
#MZ_Zea          3.144651  1.0663481  5.222954
#DZ_Zea          2.572438  1.5463779  3.598498
#EZ_Zea          6.511062  4.6307016  8.391423

# rlog Ara ref net (277 genes)
#                Zsummary  Zdensity     Zconn
#DZ_Arabidopsis 17.742067 12.183208 23.300926
#EZ_Arabidopsis 43.255858 35.884971 50.626746
#MZ_Cucumis     12.927381  9.670542 16.184220
#DZ_Cucumis      9.645396  7.226562 12.064229
#EZ_Cucumis     10.098909  8.253608 11.944210
#MZ_Glycine      6.440306  4.331749  8.548863
#DZ_Glycine     11.772110  8.379467 15.164753
#EZ_Glycine      6.499274  4.156377  8.842171
#MZ_Oryza       11.855119  8.991717 14.718521
#DZ_Oryza        4.818942  2.098337  7.539547
#EZ_Oryza        3.357615  1.146744  5.568487
#MZ_Solanum     21.458302 17.595319 25.321286
#DZ_Solanum      9.935214  7.431170 12.439257
#EZ_Solanum      6.895197  4.446847  9.343546
#MZ_Zea          5.058125  2.101936  8.014314
#DZ_Zea          4.309930  3.435078  5.184782
#EZ_Zea          6.175596  3.913725  8.437467

## rlog  All mz expressed and all PLT TFBS (42 genes)
#                Zsummary  Zdensity     Zconn
#DZ_Arabidopsis 19.977804 13.533004 26.422604
#EZ_Arabidopsis 37.889987 31.927226 43.852748
#MZ_Cucumis      6.645665  4.288800  9.002531
#DZ_Cucumis     14.102583 11.588691 16.616475
#EZ_Cucumis      9.382936  7.800779 10.965093
#MZ_Glycine      3.343264  1.361540  5.324988
#DZ_Glycine      8.218757  5.835408 10.602105
#EZ_Glycine      7.557719  4.874728 10.240709
#MZ_Oryza       12.197060  9.126606 15.267513
#DZ_Oryza        6.978280  4.549592  9.406967
#EZ_Oryza        3.861792  2.178007  5.545578
#MZ_Solanum     32.444020 26.694587 38.193453
#DZ_Solanum      8.793021  5.676448 11.909594
#EZ_Solanum     10.316235  7.125908 13.506562
#MZ_Zea          4.535489  1.775825  7.295153
#DZ_Zea          3.727648  1.372713  6.082583
#EZ_Zea          6.127573  3.616684  8.638461