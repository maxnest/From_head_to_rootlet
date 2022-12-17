setwd("D:/Projects/Genes_level/Rhizocephala/MolSign_renewed/ExSec/Expr_analysis/ES_in_molsign/")
### Libraries ###
library(dplyr)
library(phylotools)

### Input files ###
class_es <- get.fasta.name("../../Pret_class_ExSec_with_noticeable_expr.fasta")
nonclass_es <- get.fasta.name("../../Pret_nonclass_ExSec_with_noticeable_expr.fasta")
molsign_tab <- read.csv2("../../../Expr_analysis/MolSign_2TPM/Pret_ref.molsign_summary.threshold_2TPM.tsv", header = T, sep="\t")

### Processing ###
class_es <-  gsub("_i.*", "", class_es)
nonclass_es <-  gsub("_i.*", "", nonclass_es)

### Analysis ###
colnames(molsign_tab)
samples <- c("Externa", "Growing", "Middle", "Terminal")

# Classical ES: 
for (sample in samples){
  class_sample <- subset(molsign_tab, molsign_tab[sample] == 1 & GeneIDs %in% class_es)
  write.table(class_sample$GeneIDs, file=sprintf("Pret_classical_ExSec_in_%s_molsign.2TPM.txt", sample), col.names = F, row.names = F, sep="\t")
}

# Nonclassical ES:
for (sample in samples){
  nonclass_sample <- subset(molsign_tab, molsign_tab[sample] == 1 & GeneIDs %in% nonclass_es)
  write.table(nonclass_sample$GeneIDs, file=sprintf("Pret_nonclassical_ExSec_in_%s_molsign.2TPM.txt", sample), col.names = F, row.names = F, sep="\t")
}