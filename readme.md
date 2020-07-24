The goal of this program is to take the output of a RNAseq pipeline, which is a count matrix, and enable the user to analyze differentially expressed genes.

Files:
+ "CountMatrix_Tribolium.csv" is the csv file of gene counts output from featureCounts.
+ "DESeq2_DifferentialGeneExpression.R" is an R file that contains code to analyze gene expression across samples.
      This file also contains tools to make a basic volcano plot, PCA, and MAPlot.
+ "dsEve RNAi_PCAplot_PCAgenecontributions.R" is an R file that contains code to visualize gene contributions to the PCA.
+ "Triboliumsample_PCA.png" is the PCA output of "DESeq2_DifferentialGeneExpression.R".
+ "volcano_reveal.gif" is a GIF file that is the output of the "volcanoplot.R" file.
+ "volcano_total.png" is one of the outputs of "volcanoplot.R".
+ "volcanoplot.R" is an R file that contains the code for ggplot2 volcano plots and transitions.
