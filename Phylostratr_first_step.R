setwd("D:/Projects/Genes_level/Rhizocephala/MolSign_renewed/Phylostratr/")
### Libraries ###
library(phylostratr)
library(reshape2)
library(taxizedb)
library(dplyr)
library(readr)
library(magrittr)

### Peltogaster reticulata Shiino, 1943 NCBI Taxonomy ID: 2664287 ###
focal_taxid <- '2664287'
strata <- 
  # Get stratified relativies represented in UniProt (from=1 - cellular organisms)
  uniprot_strata(focal_taxid, from=1) %>% 
  # Customize clade tree: Get diverse sample of 5 species from each stratum
  strata_apply(f=diverse_subtree, n=5, weights=uniprot_weight_by_ref()) %>% 
  # Use prebuilt set of prokaryotes; add human and yeast
  use_recommended_prokaryotes %>%
  add_taxa(c('4932', '9606')) %>%
  # Download genomes, storing the filenames 
  uniprot_fill_strata

pdf(file="Preticulata_phylostratr_tree_plot.before_customize.pdf")
strata %>% strata_convert(target='all', to='name') %>% sort_strata %>% plot(cex=0.2, no.margin=TRUE, label.offset=1)
dev.off()

# Customize proteome: add P.reticulata data
strata@data$faa[[focal_taxid]] <- "/Projects/Genes_level/Rhizocephala/MolSign_renewed/Public_100aa/Preticulata_ref.genes_level.after_filters.prot.fasta"
# Customize proteome: add Armadillidium vulgare proteome from UniProt (100aa)
strata <- add_taxa(strata, "13347")
strata@data$faa[["13347"]] <- "/Projects/Genes_level/Rhizocephala/MolSign_renewed/Public_100aa/Armadillidium_vulgare_UP000288706_13347.100aa.fasta"
# Customize proteome: add Armadillidium nasatum proteome form UniProt (100aa):
strata <- add_taxa(strata, "96803")
strata@data$faa[["96803"]] <- "/Projects/Genes_level/Rhizocephala/MolSign_renewed/Public_100aa/Armadillidium_nasatum_UP000326759_96803.100aa.fasta"
# Customize proteome: add Daphnia magna proteome from UniProt (100aa)
strata <- add_taxa(strata, "35525")
strata@data$faa[["35525"]] <- "/Projects/Genes_level/Rhizocephala/MolSign_renewed/Public_100aa/Daphnia_magna_UP000076858_35525.100aa.fasta"
# Customize proteome: add Daphnia pulex proteome from UniProt (100aa):
strata <- add_taxa(strata, "6669")
strata@data$faa[["6669"]] <- "/Projects/Genes_level/Rhizocephala/MolSign_renewed/Public_100aa/Daphnia_pulex_UP000000305_6669.100aa.fasta"
# Customize proteome: add Penaeus vannamei proteome form UniProt (100aa):
strata <- add_taxa(strata, "6689")
strata@data$faa[["6689"]] <- "/Projects/Genes_level/Rhizocephala/MolSign_renewed/Public_100aa/Penaeus_vannamei_UP000283509_6689.100aa.fasta"
# Customize proteome: add Portunus trituberculatus proteome form UniProt (100aa):
strata <- add_taxa(strata, "210409")
strata@data$faa[["210409"]] <- "/Projects/Genes_level/Rhizocephala/MolSign_renewed/Public_100aa/Portunus_trituberculatus_UP000324222_210409.100aa.fasta"
# Customize proteome: add Tigriopus californicus proteome form UniProt (100aa):
strata <- add_taxa(strata, "6832")
strata@data$faa[["6832"]] <- "/Projects/Genes_level/Rhizocephala/MolSign_renewed/Public_100aa/Tigriopus_californicus_UP000318571_6832.100aa.fasta"
# Customize proteome: add Amphibalanus amphitrite proteome form UniProt (100aa):
strata <- add_taxa(strata, "1232801")
strata@data$faa[["1232801"]] <- "/Projects/Genes_level/Rhizocephala/MolSign_renewed/Public_100aa/Amphibalanus_amphitrite_UP000440578_1232801.100aa.fasta"

pdf(file="Preticulata_phylostratr_tree_plot.after_customize.pdf")
strata %>% strata_convert(target='all', to='name') %>% sort_strata %>% plot(cex=0.2, no.margin=TRUE, label.offset=1)
dev.off()

### Save image ###
save.image(file="Preticulata_phylostratr_new.RData")
