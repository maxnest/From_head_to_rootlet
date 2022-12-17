setwd("D:/Projects/Genes_level/Rhizocephala/MolSign_renewed/Expr_analysis/DiffExpr/")
### Libraries ###
library(RNentropy)
library(dplyr)

### Input data ###
expr_unaveraged <- read.csv2("../../Set_of_genes/Preticulata_ref.genes_level.after_filters.unaveraged_TPMs.tsv", header = T, sep="\t", row.names = 1)
expr_averaged <- read.csv2("../../Set_of_genes/Preticulata_ref.genes_level.after_filters.averaged_TPMs.tsv", header = T, sep="\t", row.names = 1)
design <- read.csv2("Pret_design.externa_and_interna_only.txt", header = T, sep="\t", row.names = 1)
species_tag <- "Pret"

### Processing ###
# high expression - 2TPM at least in one sample #
expr_averaged_mutate <- mutate_all(expr_averaged, function(x) as.numeric(as.character(x)))
str(expr_averaged_mutate)
expr_averaged_mutate_selected <- select(expr_averaged_mutate, Externa, Growing, Middle, Terminal)
high_expr <- filter(expr_averaged_mutate_selected, Externa >= 2 | Growing >= 2 | Middle >= 2 | Terminal >= 2)
# filters applying #
expr_unaveraged_filtered <- subset(expr_unaveraged, rownames(expr_unaveraged) %in% rownames(high_expr))
expr_unaveraged_filtered_selected <- select(expr_unaveraged_filtered, rownames(design))
# Mutate #
expr_mut <- mutate_all(expr_unaveraged_filtered_selected, function(x) as.numeric(as.character(x)))
expr_mut$GeneIDs <- as.factor(rownames(expr_mut))
expr_mut_reorder <- select(expr_mut, GeneIDs, rownames(design))
design_matrix <- data.matrix(design, rownames.force = T)

### Analysis ###
# compute statistics and p-values
RNresults <- RN_calc(expr_mut_reorder, design = design_matrix)
# select only genes with significant changes of expression
# select sequences with global p-value lower than an user defined threshold 
# and provide a summary of over- and under-expression accoding to local p-values
RNresults_selected <- RN_select(RNresults, gpv_t = 0.01, lpv_t = 0.01, method = 'BH')
######
# NB! Был добавлен новый df - selected:
# Transcripts/genes with a corrected global p-value lower than gpv_t. 
# For each condition it will contain a column where values can be -1,0,1 or NA. 
# 1 means that all the replicates of this condition have expression value higher than the
# average and local p-value <= lpv_t (thus the corresponding gene will be over-expressed in this condition). 
# -1 means that all the replicates of this condition have expression value lower than the average and local p-value <= lpv_t (thus
# the corresponding gene will be under-expressed in this condition). 
# 0 means that at least one of the replicates has a local p-value > lpv_t. 
# NA means that the local p-values of the replicates are not consistent for this condition, that is, at least
# one replicate results to be over-expressed and at least one results to be under-expressed.
#######
RNresults_significant <- RNresults_selected$selected
write.table(RNresults_significant, file=sprintf("%s_RNentropy_significant_results.tsv", species_tag),
            sep="\t", col.names = T, row.names = F)

### Subset over-expressed genes only ###
colnames(RNresults_significant)

for (sample in c("Externa", "Growing_part", "Main_part", "Thoracic_part")){
  sample_overexpr <- subset(RNresults_significant, RNresults_significant[[sample]] == "1")
  write.table(sample_overexpr, file=sprintf("%s_%s_RNentropy_overexp.tsv", species_tag, sample),  sep="\t", col.names = T, row.names = F)
}

### Compute point mutual information matrix for the experimental conditions ###
RNresults_pmi <- RN_pmi(RNresults)
RNresults_npmi_matrix <- RNresults_pmi$npmi
write.table(RNresults_npmi_matrix, file=sprintf("%s_RNentropy_npmi.tsv", species_tag), sep="\t", col.names = T, row.names = T)
