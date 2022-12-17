setwd("D:/Projects/Genes_level/Rhizocephala/MolSign_renewed/GSEA_Dmelanogaster/MolSign_wordclouds/")
### Libraries ###
library(wordcloud)
library(dplyr)
library(viridis)

### Input files ###
externa_reduced <- read.csv2("../MolSign/P.reticulata_Externa_all_genes_reduced_GSEA_results.Dmelanogaster.tsv", header = T, sep="\t")
growing_reduced <- read.csv2("../MolSign/P.reticulata_Growing_all_genes_reduced_GSEA_results.Dmelanogaster.tsv", header = T, sep="\t")
main_reduced <- read.csv2("../MolSign/P.reticulata_Middle_all_genes_reduced_GSEA_results.Dmelanogaster.tsv", header = T, sep="\t")
thoracic_reduced <- read.csv2("../MolSign/P.reticulata_Terminal_all_genes_reduced_GSEA_results.Dmelanogaster.tsv", header = T, sep="\t")
species_tag <- "P.reticulata"
set_tag <- "molecular_signatures"

set.seed(1234)

### Processing ###
input_list <- list("externa" = externa_reduced,
                   "growing_part" = growing_reduced, 
                   "main_part" = main_reduced, 
                   "thoracic_part" = thoracic_reduced)

for (set in 1:length(input_list)){
  set_reduced_terms <- count(input_list[[set]], input_list[[set]]["parentTerm"], sort=T)
  colnames(set_reduced_terms) <- c("word", "freq")
  viridis_colors <- viridis(n = length(set_reduced_terms$word))
  # visual #
  pdf(file=sprintf("%s_%s_%s.wordcloud.pdf", species_tag, names(input_list[set]), set_tag), width = 9, height = 9)
  set_reduced_term_cloud <- wordcloud(words=set_reduced_terms$word, freq=set_reduced_terms$freq, min.freq = 1, 
                                      max.words=length(set_reduced_terms$word), scale = c(4.25, 0.25),
                                      random.order = F, rot.per = 0, colors = as.character(viridis_colors), ordered.colors = T)
  dev.off()
}