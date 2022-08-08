# R Script to export counts with metadata from .rds

# For questions: Please contact Omar Laham

#import DESeq2 to normalize counts
library(DESeq2)


dds <- readRDS("002_dds.rds")

#pre-normalized
countsData <- dds@assays@data$counts
meta <- dds@colData


# remove: no need to do again, just use counts() from DESeq2
### Check that sample names match in both files
#all(colnames(countsData) %in% rownames(meta))
#all(colnames(countsData) == rownames(meta))
#normalize counts using DESeq2
#dds_matrix <- DESeqDataSetFromMatrix(countData = countsData, colData = meta, design = ~ Condition)
#dds_matrix <- estimateSizeFactors(dds_matrix)
# /remove

normCounts <- as.data.frame(counts(dds, normalized=TRUE))

#transpose
genes <- rownames(dds)
normCounts <- as.data.frame(t(normCounts))
colnames(normCounts) <- genes

#add metadata to normCounts
#cell_line
normCounts$cell_line <- dds@colData@listData$Cell_line
#batch
normCounts$batch <- dds@colData@listData$Batch
#condition
normCounts$condition <- dds@colData@listData$Condition

#export tsv
write.table(normCounts, file='ips_norm_expression.tsv', sep='\t', row.names = T, col.names = T)

