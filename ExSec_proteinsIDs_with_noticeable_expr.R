setwd("D:/Projects/Genes_level/Rhizocephala/MolSign_renewed/ExSec/")
### Libraries ###
library(dplyr)
library(phylotools)

### Input data ###
class_exsec <- get.fasta.name(infile = "TMHMM/Classic/Pret_ref.genes.potential_classical_exseq.without_mTP_and_tmhelices.fasta")
nonclass_exsec <- get.fasta.name(infile = "TMHMM/Nonclassic/Pret_ref.genes.potential_nonclassical_exseq.without_mTP_and_tmhelices.fasta")
expr <- read.csv2("../Set_of_genes/Preticulata_ref.genes_level.after_filters.averaged_TPMs.tsv", header = T, sep="\t")

### Processing ###
expr_mut <- mutate_all(expr, function(x) as.numeric(as.character(x)))
expr_mut$GeneIDs <- expr$GeneIDs
colnames(expr_mut)
expr_mut <- select(expr_mut, "GeneIDs", "Externa", "Growing", "Middle", "Terminal")
expr_mut_subset <- subset(expr_mut, Externa >= 2 | Growing >= 2 | Middle >= 2 | Terminal >= 2)

class_df <- data.frame("ProtIDs" = class_exsec, "GeneIDs" = gsub("_i.*", "", class_exsec))
nonclass_df <- data.frame("ProtIDs" = nonclass_exsec, "GeneIDs" = gsub("_i.*", "", nonclass_exsec))

class_exsec_expr <- select(subset(class_df, GeneIDs %in% expr_mut_subset$GeneIDs), ProtIDs)
nonclass_exsec_expr <- select(subset(nonclass_df, GeneIDs %in% expr_mut_subset$GeneIDs), ProtIDs)

### Output ###
write.table(class_exsec_expr, file="Pret_ref.potential_class_ExSec_with_noticeable_expr.protIDs.txt", col.names = F, row.names = F, sep="\t")
write.table(nonclass_exsec_expr, file="Pret_ref.potential_nonclass_ExSec_with_noticeable_expr.protIDs.txt", col.names = F, row.names = F, sep="\t")