#DESeq2 for differential gene expression
#This code takes a count matrix and will allow for differential gene expression analysis and create volcano plots.
#Help and troubleshooting: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

#Author: Robert Porter
#Date : 7/1/20

#load DESeq2 and read count matrix
library(DESeq2)
readcounts <- read.csv("CountMatrix_Tribolium.csv", header=TRUE)

#set row names
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

#Alternate gene labels from a csv file
labels <- read.csv("gene_names.csv", header=TRUE)
res1 <- as.data.frame(DGE.results)
volcano_dataframe <- cbind(res1,labels)
volcano_dataframe_ordered <- volcano_dataframe[order(volcano_dataframe$padj), ]
  
#basic volcano plot
library(EnhancedVolcano)
EnhancedVolcano(volcano_dataframe_ordered,
                lab = volcano_dataframe_ordered$description,
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-8, 8),
                title = 'Eve vs CTRL',
                pCutoff = .05,
                FCcutoff = 0.0,
                pointSize = 1.0,
                labSize = 2.0,
                col=c('black', 'black', 'black', 'royalblue'),
                colAlpha = 0.75)


#more advanced volcano plot

#Load libraries
library(ggplot2)
library(ggrepel)
library(scales)
library(gganimate)
#Volcano Plots of Differentially expressed genes; signifance threshold can be changed as needed.
ggplot(volcano_dataframe_ordered) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = padj < 0.05)) +
  geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label = ifelse(padj<0.1E-50, description,""))) +
  ggtitle("Eve RNAi overexpression") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

#more simplified background and axis
ggplot(volcano_dataframe_ordered) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = padj < 0.05)) +
  geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label = ifelse(padj<1E-50, description,""))) +
  ggtitle("Eve RNAi overexpression") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey92"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey92"),
        legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

#transition_reveal() (slider up and down)
#transitions need omitted cells gone
omitted_volcanodf <-na.omit(volcano_dataframe_ordered)
ggplot(data = omitted_volcanodf,
       aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = padj < 0.05)) +
  geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label = ifelse(padj<1E-50, description,""))) +
  ggtitle("Eve RNAi overexpression") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  transition_time(-log10(padj)) +
  shadow_mark()
anim_save(filename = "volcano_reveal.gif")
