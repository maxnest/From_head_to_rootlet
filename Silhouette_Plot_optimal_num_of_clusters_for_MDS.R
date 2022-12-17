setwd("D:/Projects/Genes_level/Rhizocephala/MolSign_renewed/MDS/Num_of_clusters/")
### Library ###
library(factoextra)
library(dplyr)

### Input ###
pres_abs_matrix <- read.csv2("../Input_matrix/Preticulata_expr_pres_abs_matrix.threshold_2TPM.tsv", header = T, row.names = 1, sep="\t")
species_tag <- "P.reticulata"

### Processing ###
inv_rows <- apply(pres_abs_matrix, 1, function(x) all(x == 0) | all(x == 1))
pres_abs_matrix_subset <- pres_abs_matrix[!inv_rows, ]
pres_abs_matrix_subset_tr <- t(pres_abs_matrix_subset) # cols = genes, rows = samples

### Analysis ###
set.seed(1234)
pdf(file=sprintf("%s_samples_silhouette_method.pdf", species_tag), width = 10, height = 10)
fviz_nbclust(pres_abs_matrix_subset_tr, kmeans, method = "silhouette", k.max = nrow(pres_abs_matrix_subset_tr) - 1) + 
  theme_minimal() + ggtitle(sprintf("The Silhouette Plot:: %s samples", species_tag))
dev.off()
