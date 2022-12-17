setwd("D:/Projects/Genes_level/Rhizocephala/MolSign_renewed/ExSec/Expr_analysis/OverExpr/")
### Libraries ###
library(dplyr)
library(phylotools)

### Input files ###
class_es <- get.fasta.name("../../Pret_class_ExSec_with_noticeable_expr.fasta")
nonclass_es <- get.fasta.name("../../Pret_nonclass_ExSec_with_noticeable_expr.fasta")
rnentropy_tab <- read.csv2("../../../Expr_analysis/DiffExpr/Pret_RNentropy_significant_results.tsv", header = T, sep="\t")

### Processing ###
class_es <-  gsub("_i.*", "", class_es)
nonclass_es <-  gsub("_i.*", "", nonclass_es)

### Analysis ###
colnames(rnentropy_tab)
samples <- c("Externa", "Growing_part", "Main_part", "Thoracic_part")

# Classical ES: 
for (sample in samples){
  class_over <- subset(rnentropy_tab, rnentropy_tab[sample] == 1 & GeneIDs %in% class_es)
  write.table(class_over$GeneIDs, file=sprintf("Pret_classical_ExSec_in_%s_overexpr.2TPM.txt", sample), col.names = F, row.names = F, sep="\t")
}

# Nonclassical ES:
for (sample in samples){
  nonclass_over <- subset(rnentropy_tab, rnentropy_tab[sample] == 1 & GeneIDs %in% nonclass_es)
  write.table(nonclass_over$GeneIDs, file=sprintf("Pret_nonclassical_ExSec_in_%s_overexpr.2TPM.txt", sample), col.names = F, row.names = F, sep="\t")
}