setwd("D:/Projects/Genes_level/Rhizocephala/MolSign_renewed/MDS/Input_matrix/")
### Library ###
library(dplyr)

### Input ###
diffexpr <- read.csv2("../../Expr_analysis/DiffExpr/Pret_RNentropy_significant_results.tsv", header = T, sep="\t", row.names = 1)
species_tag <- "Preticulata"

### Processing ###
diffexpr_selected <- diffexpr[, 3:ncol(diffexpr)]
# From df to matrix:
diffexpr_selected_matrix <- as.matrix(diffexpr_selected)
# Convert expression to abs/pres (0/1):
diffexpr_selected_matrix[diffexpr_selected_matrix < 1] = 0
# Filtering #
diffexpr_selected_matrix_filtered <- na.omit(diffexpr_selected_matrix)

### Output writing ###
write.table(diffexpr_selected_matrix_filtered, file=sprintf("%s_overexpr_pres_abs_matrix.tsv", species_tag), sep="\t", col.names = T, row.names = T)
