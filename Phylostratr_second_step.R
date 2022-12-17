setwd("/media/sf_Space/Crustacea_renewed/Phylostratr/")
load("Preticulata_phylostratr_new.RData")
### Library ###
library(dplyr)

### Customize datasets ###
# Customize proteome: add P.reticulata data
strata@data$faa[[focal_taxid]] <- "/media/sf_Space/Crustacea_renewed/Public_100aa/Preticulata_ref.genes_level.after_filters.prot.fasta"
# Customize proteome: add Armadillidium vulgare proteome from UniProt (100aa)
strata <- add_taxa(strata, "13347")
strata@data$faa[["13347"]] <- "/media/sf_Space/Crustacea_renewed/Public_100aa/Armadillidium_vulgare_UP000288706_13347.100aa.fasta"
# Customize proteome: add Daphnia magna proteome from UniProt (100aa)
strata <- add_taxa(strata, "35525")
strata@data$faa[["35525"]] <- "/media/sf_Space/Crustacea_renewed/Public_100aa/Daphnia_magna_UP000076858_35525.100aa.fasta"
# Customize proteome: add Daphnia pulex proteome from UniProt (100aa):
strata <- add_taxa(strata, "6669")
strata@data$faa[["6669"]] <- "/media/sf_Space/Crustacea_renewed/Public_100aa/Daphnia_pulex_UP000000305_6669.100aa.fasta"
# Customize proteome: add Penaeus vannamei proteome form UniProt (100aa):
strata <- add_taxa(strata, "6689")
strata@data$faa[["6689"]] <- "/media/sf_Space/Crustacea_renewed/Public_100aa/Penaeus_vannamei_UP000283509_6689.100aa.fasta"
# Customize proteome: add Portunus trituberculatus proteome form UniProt (100aa):
strata <- add_taxa(strata, "210409")
strata@data$faa[["210409"]] <- "/media/sf_Space/Crustacea_renewed/Public_100aa/Portunus_trituberculatus_UP000324222_210409.100aa.fasta"
# Customize proteome: add Tigriopus californicus proteome form UniProt (100aa):
strata <- add_taxa(strata, "6832")
strata@data$faa[["6832"]] <- "/media/sf_Space/Crustacea_renewed/Public_100aa/Tigriopus_californicus_UP000318571_6832.100aa.fasta"
# Customize proteome: add Armadillidium nasatum proteome form UniProt (100aa):
strata <- add_taxa(strata, "96803")
strata@data$faa[["96803"]] <- "/media/sf_Space/Crustacea_renewed/Public_100aa/Armadillidium_nasatum_UP000326759_96803.100aa.fasta"
# Customize proteome: add Amphibalanus amphitrite proteome form UniProt (100aa):
strata <- add_taxa(strata, "1232801")
strata@data$faa[["1232801"]] <- "/media/sf_Space/Crustacea_renewed/Public_100aa/Amphibalanus_amphitrite_UP000440578_1232801.100aa.fasta"

### Run BLASTp ###
strata <- strata_blast(strata, blast_args = list(nthreads=3)) %>% strata_besthits

# Merge results into a single hittable
results <- merge_besthits(strata)
write.table(results, file="Preticulata_100aa_phylostratr_besthits.tsv", sep="\t", row.names=F)

# Infer phylostrata
ph <- stratify(results)
write.table(ph, file="Preticulata_phylostratr_results.tsv", sep="\t", row.names = F)

# And get the number of genes in each phylostrata
ph$locus <- sub('\\.[0-9]+', '', ph$qseqid, perl=TRUE)
ph %>%
  dplyr::select(-qseqid) %>%
  dplyr::distinct() %>%
  dplyr::group_by(mrca_name, ps) %>%
  dplyr::summarize(n = length(ps))
write.table(table(ph$mrca_name), file="Preticulata_phylostratr_table.tsv", sep="\t", row.names = F, col.names = F)
