#PCA Plot of rlog transformed data from the DESeq2
#This code takes a count matrix and will allow for qunatification and visualization of PCA contributions
#Help and troubleshooting:
#https://tavareshugo.github.io/data-carpentry-rnaseq/03_rnaseq_pca.html
#https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html

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

#plotPCA -> original graph
library (ggplot2)
plotPCA (DESeq.rlog)
###################################################

# load packages
library(tidyverse)
library(broom)

# Perform the PCA
#prcomp() creates a prcomp element with two dataframes: "x" and "rotation"
#The github tutorial uses transform dataframe, t(), and does not standardize the data.
p <- prcomp(t(rlog.counts))
summary(p)
head(p$rotation)
head(p$x)

# PC variances of standardized data (eigen values)
tidy(p, matrix = "pcs")

tidy(p, matrix = "pcs") %>% 
  ggplot(aes(x = factor(PC))) +
  geom_col(aes(y = percent)) +
  geom_line(aes(y = cumulative, group = 1)) + 
  geom_point(aes(y = cumulative)) +
  labs(x = "Principal component", y = "Fraction variance explained")

# extract values
p$rotation %>% 
  # convert to a tibble
  as_tibble(rownames = "gene")

# variable loadings (eigen vectors) for genes that were used as variables in PCA
tidy(p, matrix = "variables")
head(p)

#extract most relevant genes

top_genes <- p %>% 
  # extract variable (gene) loadings
  tidy(matrix = "variables") %>%  
  # retain only PC1 and PC2
  filter(PC %in% c(1, 2)) %>%
  # for each PC
  group_by(PC) %>%
  # sort descending value
  arrange(desc(abs(value))) %>%
  # take top 5 rows of each PC group
  slice(1:5) %>%
  # extract the column (gene name) from the table
  pull(column) %>%
  # retain unique gene names only
  unique()

top_genes


#Plot loading data

gene_loadings <- p$rotation %>% 
  as_tibble(rownames = "gene") %>% 
  filter(gene %in% top_genes)

loadings_plot <- ggplot(gene_loadings) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               colour = "brown") +
  geom_text(data = gene_loadings, 
            aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))

loadings_plot
