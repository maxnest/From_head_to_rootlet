setwd("D:/Projects/Genes_level/Rhizocephala/MolSign_renewed/MDS/Input_matrix/")
### Library ###
library(dplyr)

### Input ###
unaver_expr <- read.csv2("../../Set_of_genes/Preticulata_ref.genes_level.after_filters.unaveraged_TPMs.tsv", header = T, sep="\t", row.names = 1)
aver_expr <- read.csv2("../../Set_of_genes/Preticulata_ref.genes_level.after_filters.averaged_TPMs.tsv", header = T, sep="\t", row.names = 1)
species_tag <- "Preticulata"

### Processing ###
# Chr -> Num:
unaver_expr_mut <- mutate_all(unaver_expr[,1:ncol(unaver_expr)], function(x) as.numeric(as.character(x)))
aver_expr_mut <- mutate_all(aver_expr[,1:ncol(aver_expr)], function(x) as.numeric(as.character(x)))
# Select and rename :
colnames(unaver_expr_mut)
unaver_expr_mut_select <- select(unaver_expr_mut, "Externa_first", "Externa_second", "Growing_first", "Growing_second", 
                                 "Middle_first", "Middle_second", "Terminal_first", "Terminal_second")
colnames(aver_expr_mut)
aver_expr_mut_select <- select(aver_expr_mut, "Externa", "Growing", "Middle", "Terminal")

colnames(unaver_expr_mut_select) <- c("Externa_rep1", "Externa_rep2", "Growing_part_rep1", "Growing_part_rep2", 
                                      "Main_part_rep1", "Main_part_rep2", "Thoracic_part_rep1", "Thoracic_part2")
# Only genes with noticeable expression (>= 2TPM):
aver_expr_mut_select_subset <- subset(aver_expr_mut_select, Externa >= 2 | Growing >= 2 | Middle >= 2 | Terminal >= 2)

unaver_expr_mut_subset <- subset(unaver_expr_mut_select, rownames(unaver_expr_mut_select) %in% rownames(aver_expr_mut_select_subset))

# From df to matrix:
unaver_expr_mut_subset_matrix <- as.matrix(unaver_expr_mut_subset)
# Convert expression to abs/pres (0/1):
# Threshold = 2TPM
unaver_expr_mut_subset_matrix[unaver_expr_mut_subset_matrix < 2] = 0
unaver_expr_mut_subset_matrix[unaver_expr_mut_subset_matrix >= 2] = 1
### Output writing ###
write.table(unaver_expr_mut_subset_matrix, file=sprintf("%s_expr_pres_abs_matrix.threshold_2TPM.tsv", species_tag), 
            sep="\t", col.names = T, row.names = T)