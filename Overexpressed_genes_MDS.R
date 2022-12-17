setwd("D:/Projects/Genes_level/Rhizocephala/MolSign_renewed/MDS/Clustering/")
### Libraries ###
library(vegan)
library(dplyr)
library(ggplot2)
library(ggpubr)

### Input data ###
pres_abs_matrix <- read.csv2("../Input_matrix/Preticulata_overexpr_pres_abs_matrix.tsv", header = T, row.names = 1, sep="\t")
cluster_num <- 2
species_tag <- "P.reticulata"

options(ggrepel.max.overlaps = Inf)
options(ggrepel.max.iter = 100000)
options(ggrepel.force = 100000)

set.seed(1234)

### Processing ###
inv_rows <- apply(pres_abs_matrix, 1, function(x) all(x == 0) | all(x == 1))
pres_abs_matrix_subset <- pres_abs_matrix[!inv_rows, ]
pres_abs_matrix_subset_tr <- t(pres_abs_matrix_subset) # cols = genes, rows = samples

### Analysis ###
mds <- metaMDS(comm = pres_abs_matrix_subset_tr, distance="manhattan", k=2, try=100, trymax=100000, autotransform = FALSE, binary=TRUE)
### Stress ###
# Stress. Not an assumption of the process, but the most important factor to consider after generating an MDS plot is the ‘stress’. 
# The stress provides a measure of the degree to which the distance between samples in reduced dimensional space (usually 2-dimensions) 
# corresponds with the actual multivariate distance between the samples. Lower stress values indicate greater conformity and therefore are desirable. 
# High stress values indicate that there was no 2-dimensional arrangement of your points that reflect their similarities. 
# A rule of thumb is that stress values should ideally be less than 0.2 or even 0.1.
mds$stress

### Visual ###
mds_xy <- data.frame(mds$points)
rownames(mds_xy)
mds_xy$Sample <- as.factor(rownames(mds_xy))
str(mds_xy)

clust <- as.factor(kmeans(mds_xy[,1:2], cluster_num)$cluster)
mds_xy$Clusters <- clust

# ggscatter #
pdf(file=sprintf("%s_sample_overexpr_MDS_ggscatter.%d_clusters.7x7.pdf", species_tag, cluster_num), width = 7, height = 7)
ggscatter(mds_xy, x = "MDS1", y = "MDS2", 
          label = "Sample",
          color = "Clusters",
          palette = "jco",
          size = 8, 
          font.label = c(16, "plain"),
          ellipse = TRUE,
          ellipse.type = "convex",
          repel = TRUE, label.rectangle = TRUE,
          title=sprintf("MDS plot for the over-expressed genes in %s female body parts", species_tag))  
dev.off()