setwd("D:/Projects/Genes_level/Rhizocephala/Analysis/new_Quant")
### Libraries ###
library(tximport)
library(dplyr)

### Input data ###
# Salmon quant
Pret_quant_files <- c("Pret_ref_salmon_results/Pret_externa_first_salmon_quant/quant.sf",
                      "Pret_ref_salmon_results/Pret_externa_second_salmon_quant/quant.sf",
                      "Pret_ref_salmon_results/Pret_growing_first_salmon_quant/quant.sf",
                      "Pret_ref_salmon_results/Pret_growing_second_salmon_quant/quant.sf",
                      "Pret_ref_salmon_results/Pret_middle_first_salmon_quant/quant.sf",
                      "Pret_ref_salmon_results/Pret_middle_second_salmon_quant/quant.sf", 
                      "Pret_ref_salmon_results/Pret_terminal_first_salmon_quant/quant.sf",
                      "Pret_ref_salmon_results/Pret_terminal_second_salmon_quant/quant.sf",
                      "Pret_ref_salmon_results/Pret_whole_body_first_salmon_quant/quant.sf",
                      "Pret_ref_salmon_results/Pret_whole_body_second_salmon_quant/quant.sf")

# Contig-to-Gene Correspondence Table based on Trinity assembling results
Pret_transcripts_2_genes <- read.csv2("Pret_ref_decont_trans_2_genes/Pret_ref_decont.trans_gene_map",
                                      header = T, sep="\t")
### Analysis ###
# From contigs/transcripts expr levels to genes expr levels  #
Pret_genes <- tximport(Pret_quant_files, type="salmon", txIn=T, txOut = F, countsFromAbundance = "no", txIdCol = "Name",
                       abundanceCol = "TPM", countsCol = "NumReads", lengthCol = "Length",
                       tx2gene = Pret_transcripts_2_genes)
Pret_genes_TPMs <- as.data.frame(Pret_genes$abundance)
colnames(Pret_genes_TPMs) <- c("Externa_first", "Externa_second",
                               "Growing_first", "Growing_second", 
                               "Middle_first", "Middle_second",
                               "Terminal_first", "Terminal_second",
                               "Whole_body_first", "Whole_body_second")
Pret_genes_TPMs$GeneIDs <- rownames(Pret_genes_TPMs)
new_Pret_genes_TPMs <- select(Pret_genes_TPMs, GeneIDs, 
                              Externa_first, Externa_second,
                              Growing_first, Growing_second, 
                              Middle_first, Middle_second,
                              Terminal_first, Terminal_second, 
                              Whole_body_first, Whole_body_second)

### Output files creating ###
# All TPMs per genes #
write.table(new_Pret_genes_TPMs, file="good.Pret_ref.clustered_decont.genes_level.all_TPMs.tsv",
            sep="\t", quote = F, col.names = T, row.names = F)

# TPM values > 1 at least in one of the sample #
new_Pret_selected_genes_TPMs <- subset(new_Pret_genes_TPMs, Externa_first > 1 | Externa_second > 1 | Growing_first > 1 | Growing_second > 1 |
                                         Middle_first > 1 | Middle_second > 1 | Terminal_first > 1 | Terminal_second > 1 | 
                                         Whole_body_first > 1 | Whole_body_second > 1)
write.table(new_Pret_selected_genes_TPMs, file="good.Pret_ref.clustered_decond.genes_level.selected_TPMs.tsv", 
            sep="\t", quote = F, col.names = T, row.names = F)

# All averaged values #
new_Pret_genes_TPMs$Externa <- (new_Pret_genes_TPMs$Externa_first + new_Pret_genes_TPMs$Externa_second)/2
new_Pret_genes_TPMs$Growing <- (new_Pret_genes_TPMs$Growing_first + new_Pret_genes_TPMs$Growing_second)/2
new_Pret_genes_TPMs$Middle <- (new_Pret_genes_TPMs$Middle_first + new_Pret_genes_TPMs$Middle_second)/2
new_Pret_genes_TPMs$Terminal <- (new_Pret_genes_TPMs$Terminal_first + new_Pret_genes_TPMs$Terminal_second)/2
new_Pret_genes_TPMs$Whole_body <- (new_Pret_genes_TPMs$Whole_body_first + new_Pret_genes_TPMs$Whole_body_second)/2
new_Pret_averaged_TPMs <- select(new_Pret_genes_TPMs, GeneIDs, Externa, Growing, Middle, Terminal, Whole_body)
write.table(new_Pret_averaged_TPMs, file="good.Pret_ref.clustered_decont.genes_level.all_averaged_TPMs.tsv", 
            sep="\t", quote = F, col.names = T, row.names = F)

# Only averaged TPM values > 1 at least in one of the samples #
new_Pret_averaged_and_selected_TPMs <- subset(new_Pret_averaged_TPMs, Externa > 1 | Growing > 1 | Middle > 1 | Terminal > 1 | Whole_body > 1) 
write.table(new_Pret_averaged_and_selected_TPMs, file="good.Pret_ref.clustered_decont.genes_level.selected_averaged_TPMs.tsv",
            sep="\t", quote = F, col.names = T, row.names = F)