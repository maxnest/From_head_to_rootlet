setwd("D:/Projects/Genes_level/Rhizocephala/MolSign_renewed/GSEA/")
### Library ###
library(dplyr)

### Input files ###
emapper <- read.csv2("../vs_DBs/eggNOG-mapper_v2/Pret_100aa_emapper_v2.annotations.modified.tsv", header = T, sep="\t") # First you need to remove all lines starting with '#'
averaged_expr <- read.csv2("../Set_of_genes/Preticulata_ref.genes_level.after_filters.averaged_TPMs.tsv", header = T, row.names = 1,  sep="\t")
species_tag <- "Pret"

### Processing ###
averaged_expr_mut <- mutate_all(averaged_expr, function(x) as.numeric(as.character(x)))
colnames(averaged_expr_mut)
averaged_expr_mut_selected <- select(averaged_expr_mut,  "Externa", "Growing", "Middle", "Terminal")
averaged_expr_mut_selected_subset <- subset(averaged_expr_mut_selected, Externa >= 2 | Growing >= 2 | Middle >= 2 | Terminal >= 2)
emapper$X.query <-  gsub("_i.*", "", emapper$X.query)

### Filtering ###
emapper_subset <- subset(emapper, emapper$X.query %in% rownames(averaged_expr_mut_selected_subset) & emapper$GOs != "-")
emapper_subset_reordered <- select(emapper_subset, X.query, GOs)

### Output ###
write.table(emapper_subset_reordered, file=sprintf("%s_GeneOntologyUniverse_filtered.tsv", species_tag), sep="\t", col.names = F, row.names = F)