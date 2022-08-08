#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("DESeq2")

#commandArgs picks up the variables you pass from the command line
args <- commandArgs(trailingOnly = TRUE)
print(args)
