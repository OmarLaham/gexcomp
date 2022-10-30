#load libraries

library('DESeq2');

#install ones for homomorphic
#install.packages(c("Rcpp", "RcppParallel", "gmp"), repos = "http://cran.us.r-project.org")
#install.packages("http://www.louisaslett.com/HomomorphicEncryption/dl/HomomorphicEncryption_0.3.5.tar.gz", repos = NULL, type="source")

#library("HomomorphicEncryption")

library(homomorpheR)
keyPair <- PaillierKeyPair$new(modulusBits = 1024)
print("KeyPair:")
print(keyPair)
print("=============================")

#You can enable this for faster computing on the cloud for e.g.
#library("BiocParallel")
#register(MulticoreParam(4))


#commandArgs picks up the variables you pass from the command line
#args <- commandArgs(trailingOnly = TRUE);
runID <- 1;
winStart <- 776#args[[2]];
winEnd <- 917#args[[3]];
#dsPath <- file.path("../media/runs", runID, "data.csv");

enc <- function(x) keyPair$pubkey$encrypt(x)

dec <- function(x) keyPair$getPrivateKey()$decrypt(x)

encryptAndDecrypt <- function(x) keyPair$getPrivateKey()$decrypt(keyPair$pubkey$encrypt(x))

# DEG analysis with DESeq2

# read norm counts
normCountsPath <- file.path("../media", "runs", runID, paste0("data_norm_counts_win_", winStart, "_to_", winEnd, ".csv"))
cts <- as.matrix(read.csv(normCountsPath, header = T, sep = ",", row.names='gene'))
#encrypt cts
enc_cts <- enc(cts)#enc(k$pk, cts)

# raead experiment annotation
annoPath <- file.path("../media", "runs", runID, paste0("data_deg_annotation.csv"))
coldata <- as.matrix(read.csv(annoPath, header=T, sep=',', row.names=1))
#encrypt coldata
enc_coldata <- enc(coldata)#enc(k$pk, coldata)

#assert
rownames(enc_coldata) <- sub("fb", "", rownames(enc_coldata))
all(rownames(enc_coldata) %in% colnames(enc_cts))

enc_cts <- enc_cts[, rownames(enc_coldata)]
all(rownames(enc_coldata) == colnames(enc_cts))

## Create DESeq object
dds <- DESeqDataSetFromMatrix(countData = enc_cts, colData = enc_coldata, design = ~ group)

## Run analysis
dds <- DESeq(dds)

# DEG results
res <- results(dds)

# Note: you should use shrinkage if you want to visualize..

##remove rows with NA padj or log2FoldChange
res <- res[!is.na(res$padj),]
res <- res[!is.na(res$log2FoldChange),]
#res <- res[complete.cases(res[,c("padj", "log2FoldChange")]),]

# order by padj
res <- res[order(res$padj),]

#filter results to significant genes

res <- res[res$padj <= 0.05,]
res <- res[abs(res$log2FoldChange) >= 1.5,]
#
res["regulation"] <- "up"
#res[res$log2FoldChange < 0, "regulation"] <- "down"

#decrypt res
dec_res <- dec(res)#dec(k$sk, res)

# export results to continue with Python
savingPath <- file.path("../media", "runs", runID, paste0("homomorphic_deg_group_original_vs_prediction_win_", winStart, "_to_", winEnd, ".csv"))
write.csv(as.data.frame(dec_res), file= savingPath)













#chooseCRANmirror(ind=33); #Germany
#

#install.packages("devtools") # unnecessary if you have it already
#devtools::install_github("Bioconductor/BiocManager", ref="ghost-binary-repo")

#
# if (!require("BiocManager", quietly = TRUE)) {
# 	chooseCRANmirror(ind=77);
#     install.packages("BiocManager");
# }
# BiocManager::install("DESeq2", version = "3.15");

#
# ## ----load-libraries, message = FALSE------------------------------------------
# library("enrichR")
# library("ggplot2")
#
# #'
# #' # Load data
# #'
# #' If you have downloaded the `DESeq2_DEG.txt` file with `wget`:
# #'
# ## ----load-data----------------------------------------------------------------
# data <- data.table::fread("DESeq2_DEG.txt")
# data$GeneID <- substr(data$GeneID, 1, 18)
#
# #'
# #' If you like to donwload the file in `R` now:
# #'
# ## ----eval = FALSE-------------------------------------------------------------
# ## data <- data.table::fread("https://raw.githubusercontent.com/ycl6/GO-Enrichment-Analysis-Demo/master/DESeq2_DEG.txt")
# ## data$GeneID <- substr(data$GeneID, 1, 18)
#
# #'
# ## -----------------------------------------------------------------------------
# data
#
# #'
# #' # Define significance threshold
# #'
# ## ----significance-threshold---------------------------------------------------
# up.idx <- which(data$padj < 0.05 & data$log2fc > 0)	# FDR < 0.05 and logFC > 0
# dn.idx <- which(data$padj < 0.05 & data$log2fc < 0)	# FDR < 0.05 and logFC < 0
#
# #'
# ## -----------------------------------------------------------------------------
# dim(data)
# length(up.idx)
# length(dn.idx)
#
# #'
# #' # Define significant genes
# #'
# ## ----significant-genes--------------------------------------------------------
# up.genes <- data[up.idx,]$GeneSymbol
# dn.genes <- data[dn.idx,]$GeneSymbol
#
# #'
# ## -----------------------------------------------------------------------------
# head(up.genes, 10)
# head(dn.genes, 10)
#
# #'
# #' Alternatively, if you only have Ensembl gene ID
# #'
# ## -----------------------------------------------------------------------------
# up.genes <- data[up.idx,]$GeneID
# dn.genes <- data[dn.idx,]$GeneID
#
# #'
# ## -----------------------------------------------------------------------------
# head(up.genes, 10)
# head(dn.genes, 10)
#
# #'
# #' # Prepare input data
# #'
# #' We would need to convert any other identifier format to `SYMBOL` which is the required input identifier format. This can be done by using the `select` function from `AnnotationDbi` that we saw in [Part 1](1_Organism_DB.html) of this demo, or by using the "Biological Id TRanslator" `bitr` function from `clusterProfiler` which is a wrapper function of `AnnotationDbi::select`.
# #'
# #' Here, we will use `bitr` here to see how this can be done.
# #'
# ## ----use-bitr-----------------------------------------------------------------
# # Use fromType = "ENSEMBL" if your input identifier is Ensembl gene ID
# up.genes.df = clusterProfiler::bitr(up.genes, fromType = "ENSEMBL", toType = "SYMBOL",
# 				    OrgDb = "org.Mm.eg.db")
# head(up.genes.df, 10)
#
# dn.genes.df = clusterProfiler::bitr(dn.genes, fromType = "ENSEMBL", toType = "SYMBOL",
# 				    OrgDb = "org.Mm.eg.db")
# head(dn.genes.df, 10)
#
# #'
# #' # Connecting to Enrichr web service
# #'
# #' List available databases from Enrichr
# #'
# ## ----list-dbs-----------------------------------------------------------------
# dbs <- listEnrichrDbs()
# dbs <- dbs[order(dbs$libraryName),]
#
# class(dbs)
# dim(dbs)
# head(dbs)
#
# #'
# #' Show all database names.
# #'
# ## ----show-dbs-----------------------------------------------------------------
# dbs$libraryName
#
# #'
# #' Search for mouse databases with keyword `"Mouse"`
# #'
# ## ----show-mouse-dbs-----------------------------------------------------------
# dbs[grep("Mouse",dbs$libraryName),]$libraryName
#
# #'
# #' # Select databases
# #'
# ## ----select-dbs---------------------------------------------------------------
# dbs_go <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018")
# dbs_pw <- c("KEGG_2019_Mouse", "WikiPathways_2019_Mouse", "BioPlanet_2019")
# dbs_dd <- c("PheWeb_2019", "ClinVar_2019")
#
# #'
# #' # Perform enrichment analysis
# #'
# #' ## GO analysis
# #'
# ## ----run-enrichr-go-----------------------------------------------------------
# upEnriched_go <- enrichr(genes = up.genes.df$SYMBOL, databases = dbs_go)
# dnEnriched_go <- enrichr(genes = dn.genes.df$SYMBOL, databases = dbs_go)
#
# #'
# ## -----------------------------------------------------------------------------
# class(upEnriched_go)
# names(upEnriched_go)
#
# # View top 5 terms in the first element of the list
# head(upEnriched_go[[1]], 5)
#
# #'
# #' ## Pathway analysis
# #'
# ## ----run-enrichr-pw-----------------------------------------------------------
# upEnriched_pw <- enrichr(genes = up.genes.df$SYMBOL, databases = dbs_pw)
# dnEnriched_pw <- enrichr(genes = dn.genes.df$SYMBOL, databases = dbs_pw)
#
# #'
# ## -----------------------------------------------------------------------------
# class(upEnriched_pw)
# names(upEnriched_pw)
#
# # View top 5 terms in the first element of the list
# head(upEnriched_pw[[1]], 5)
#
# #'
# #' ## Diseases/Drugs analysis
# #'
# ## ----run-enrichr-dd-----------------------------------------------------------
# upEnriched_dd <- enrichr(genes = up.genes.df$SYMBOL, databases = dbs_dd)
# dnEnriched_dd <- enrichr(genes = dn.genes.df$SYMBOL, databases = dbs_dd)
#
# #'
# ## -----------------------------------------------------------------------------
# class(upEnriched_dd)
# names(upEnriched_dd)
#
# # View top 5 terms in the first element of the list
# head(upEnriched_dd[[1]], 5)
#
# #'
# #' # Plot enrichment
# #'
# #' Demonstrate using different paramters to plot enrichment using the `plotEnrich` function.
# #'
# ## ----plot-results, fig.width = 8, fig.height = 6, fig.align = "center", dpi = 100----
# plotEnrich(upEnriched_go[[3]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value")
# plotEnrich(upEnriched_pw[[1]], showTerms = 15, numChar = 40, y = "Ratio", orderBy = "P.value")
# plotEnrich(upEnriched_dd[[2]], showTerms = 10, numChar = 30, y = "Count", orderBy = "Combined.Score")
#
# #'
# #' # Output results to files
# #'
# #' Use the `printEnrich` function to output Enrichr results to tab-delimited text files.
# #'
# ## ----output-results-----------------------------------------------------------
# printEnrich(upEnriched_go, prefix = "enrichr-GO-up", showTerms = 20)
# printEnrich(dnEnriched_go, prefix = "enrichr-GO-dn", showTerms = 20)
# printEnrich(upEnriched_pw, prefix = "enrichr-PW-up", showTerms = 20)
# printEnrich(dnEnriched_pw, prefix = "enrichr-PW-dn", showTerms = 20)
# printEnrich(upEnriched_dd, prefix = "enrichr-DD-up", showTerms = 20)
# printEnrich(dnEnriched_dd, prefix = "enrichr-DD-dn", showTerms = 20)
#
#
# print("Done!")
