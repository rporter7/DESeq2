#DESeq2 for differential gene expression
#This code takes a count matrix and will allow for differential gene expression analysis and create useful visualizations
#Help and troubleshooting: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

#Author: Robert Porter
#Date : 7/1/20

#load DESeq2 and read count matrix
library(DESeq2)
readcounts <- read.csv("CountMatrix_input.csv", header=TRUE)

#set row names as needed per sample amount
row.names(readcounts) <- readcounts$Gene_id
readcounts <- readcounts[,2:7]
dim(readcounts)
head(readcounts, n=3) 

#Simplify Column names and set conditions
sample_info <- data.frame (row.names=names(readcounts), condition=gsub("[0-9]+", "", names(readcounts)))
sample_info

# generate the DESeqDataSet (with required fields: countData, ColData, design)
DESeq.ds <- DESeqDataSetFromMatrix (countData = readcounts, colData = sample_info, design = ~condition)

# DESeq2's default method to normalize read counts to account for differences in sequencing depths is implemented in estimateSizeFactors()
# calculate the size factor and add it to the data set
DESeq.ds <- estimateSizeFactors (DESeq.ds)

# if you check colData() again, you can see it now contains the sizeFactors
colData (DESeq.ds)

# DESeq2's rlog() function shrinks the variance of low read counts and returns values that are both normalized for sequencing depth and transformed to the log2 scale
DESeq.rlog <- rlog (DESeq.ds, blind = TRUE) #blind = FALSE, if large differences in a large proportion of genes.
rlog.counts <- assay(DESeq.rlog) #counts() will cause error because after rlog it is not counts(integer) any more
head(rlog.counts)

# set WT as the first-level reference (to be used as the denominator for the fold change calculation) 
colData(DESeq.ds)$condition <- relevel (colData(DESeq.ds)$condition, "CTRL")
colData(DESeq.ds)$condition  #note that only the level "order" is changed, while colData(DESeq.ds) is not changed

# running the DGE analysis using the DESeq() function
# DESeq() function is basically a wrapper around the following three individual functions:
# (1) DESeq.ds <- estimateSizeFactors (DESeq.ds) # sequencing depth normalization between the samples
# (2) DESeq.ds <- estimateDispersions (DESeq.ds) # gene-wise dispersion estimates 
# (3) DESeq.ds <- nbinomWaldTest (DESeq.ds) # this fits a negative binomial GLM (Logistic Regression) and applies Wald statistics to each gene

DESeq.ds <- DESeq (DESeq.ds) # output = DESeqDataSet object
#Note that the input for the DEseq () function are the raw read counts (non-normalized, untransformed), as the function will perform normalizations and transformations under the hood; supplying anything but raw read counts will result in nonsensical results.


#The results() function extracts the base means across samples, moderated log2 fold changes, standard errors, test statistics etc. for every gene.
DGE.results <- results (DESeq.ds, independentFiltering = TRUE , alpha = 0.05) #alpha is the adjusted p-value cutoff (FDR)
summary (DGE.results)


#Differential gene expression lists
upregulated <- (subset(DGE.results, padj < 0.05 & log2FoldChange > 0))
upregulated <- as.data.frame(upregulated)
write.csv(upregulated, file="Upregulated_genes.csv", quote=FALSE, row.names = TRUE)

downregulated <- (subset(DGE.results, padj < 0.05 & log2FoldChange < 0))
downregulated <- as.data.frame(downregulated)
write.csv(downregulated, file = "Downregulated_genes.csv", quote=FALSE, row.names = TRUE)


# MA plot
# y-axis: the expression change between conditions (log ratios, M); x-axis: the average expression strength of genes (average mean, 'A')
plotMA (DGE.results, alpha = 0.05 , main = "Exp vs. CTRL", ylim = c(-5,5)) # genes that pass the significance threshold (adjusted p-value < 0.05) are colored in red

#PCA Plot of the samples
library (ggplot2)
plotPCA (DESeq.rlog)

#volcano plot of the gene expression
library(EnhancedVolcano)
EnhancedVolcano(DGE.results,
                lab = rownames(DGE.results),
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-8, 8),
                title = 'Exp vs CTRL',
                pCutoff = .05,
                FCcutoff = 0.0,
                pointSize = 1.0,
                labSize = 2.0,
                col=c('black', 'black', 'black', 'royalblue'),
                colAlpha = 0.75)


#HTML widget that allows you to have an interactable plot for a target gene and see the normalized counts between samples.
library(htmlwidgets)
df <- plotCounts(DESeq.ds, gene="TC009469", intgroup="condition", returnData= TRUE)

library(ggplot2)
library(plotly)
p <- ggplot(df, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))
ggplotly(p)

#Noninteractable version of the same gene of interest normalized plot.
d <- plotCounts(DESeq.ds, gene="TC009469", intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))