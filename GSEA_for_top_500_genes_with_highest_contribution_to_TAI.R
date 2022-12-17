setwd("D:/Projects/Genes_level/Rhizocephala/MolSign_renewed/Gene_2_TAI_contibution/top500_GSEA_Dmelanogaster/")
### Libraries ###
library(topGO)
library(dplyr)
library(rrvgo)
library(ggplot2)

options(ggrepel.max.overlaps = Inf)

### Input data ###
species_geneID2GO <- readMappings(file="../../GSEA_Dmelanogaster/Pret_GeneOntologyUniverse_filtered.tsv")
species_geneNames <- names(species_geneID2GO)
species_tag <- "P.reticulata"

gene_contribution <- read.csv2("../../myTAI/P.reticulata_filtered_logTPMs.Genes_to_TAI_contribution.tsv", header = T, sep="\t")
phylostatr_tab <- read.csv2("../../Phylostratr/Preticulata_phylostratr_results.geneIDs.tsv", header = T, sep="\t")
gene_exp_filtered <- read.csv2("../../myTAI/Preticulata_ref.genes_level.after_filters.averaged_TPMs.tsv", header = T, sep="\t")
top_num <- 500

### Processing ###
gene_exp_filtered_mut <- mutate_all(gene_exp_filtered, function(x) as.numeric(as.character(x)))
gene_exp_filtered_mut$GeneIDs <- gene_exp_filtered$GeneIDs
colnames(gene_exp_filtered_mut)
gene_exp_filtered_mut <- select(gene_exp_filtered_mut, "GeneIDs", "Externa", "Growing", "Middle", "Terminal")
gene_exp_filtered_mut_subset <- subset(gene_exp_filtered_mut, Externa >= 2 | Growing >= 2 | Middle >= 2 | Terminal >= 2)
gene_exp_filtered_mut_subset <- subset(gene_exp_filtered_mut_subset, GeneIDs %in% phylostatr_tab$qseqid)

gene_IDs_df <- data.frame("GeneID_wt" = gene_exp_filtered_mut_subset$GeneIDs, "GeneID_tolower" = tolower(gene_exp_filtered_mut_subset$GeneIDs))
gene_contribution_matrix <- mutate_all(gene_contribution, function(x) as.numeric(as.character(x)))

species_geneNames_tolowest <- tolower(species_geneNames)
# Only with annotation #
gene_contribution_matrix_subset <- subset(gene_contribution_matrix, rownames(gene_contribution_matrix) %in% species_geneNames_tolowest)

### Analysis ###
for (sample in colnames(gene_contribution_matrix_subset)){
  # Top genes with highest contribution #
  sample_ordered <- gene_contribution_matrix_subset[order(-gene_contribution_matrix_subset[[sample]]), ]
  sample_top <- top_n(sample_ordered, top_num, sample_ordered[sample])
  top_genes <- select(subset(gene_IDs_df, GeneID_tolower %in% rownames(sample_top)), GeneID_wt)
  write.table(top_genes, file=sprintf("%s_%s_top500_annotated_2_TAI.tsv", species_tag, sample), col.names = F, row.names = F, sep="\t")
  # GSEA #
  sample_all_genes <- factor(as.integer(species_geneNames %in% top_genes$GeneID_wt))
  names(sample_all_genes) <- species_geneNames
  ## Biological Processes (BP) ##
  GOdata_sample_BP <- new("topGOdata", ontology="BP", allGenes=sample_all_genes,
                          annot=annFUN.gene2GO, gene2GO=species_geneID2GO)
  resultsFisher_sample_BP <- runTest(GOdata_sample_BP, algorithm="classic", statistic="fisher")
  resultsFisher_sample_BP_df <- as.data.frame(score(resultsFisher_sample_BP))
  colnames(resultsFisher_sample_BP_df) <- "P-values"
  resultsFisher_sample_BP_df_subset <- subset(resultsFisher_sample_BP_df, resultsFisher_sample_BP_df$`P-values` < 0.01)
  results_sample_BP <- GenTable(GOdata_sample_BP, classicFisher=resultsFisher_sample_BP,
                                ranksOf="classicFisher", topNodes=length(resultsFisher_sample_BP_df_subset$`P-values`))
  results_sample_BP_short <- subset(results_sample_BP, Significant >= 10)
  output_file_name <- sprintf("%s_%s_TAI_top_%d_genes_GOenrichment_BP_Fisher.min_10_genes.tsv", species_tag, sample, top_num)
  write.table(results_sample_BP_short, file=output_file_name, sep="\t", quote=F, col.names = T, row.names = F)
  # Redundancy reducing and visual #
  set_simMatrix <- calculateSimMatrix(results_sample_BP_short$GO.ID, 
                                      orgdb = "org.Dm.eg.db", 
                                      ont="BP", method="Rel")
  set_classicFisher <- gsub("< ", "", results_sample_BP_short$classicFisher)
  set_scores <- setNames(-log10(as.numeric(as.character(set_classicFisher))), results_sample_BP_short$GO.ID)
  set_reducedTerms <- reduceSimMatrix(set_simMatrix, set_scores, threshold = 0.7, orgdb="org.Dm.eg.db")
  pdf(file=sprintf("%s_%s_TAI_top_%d_genes_reduced_GSEA_results.Dmelanogaster.scatterPlot.pdf", species_tag, sample, top_num), width = 10, height = 10)
  scatter <- scatterPlot(set_simMatrix, set_reducedTerms, size="score", addLabel = TRUE, labelSize = 5)
  print(scatter)
  dev.off()
  write.table(set_reducedTerms, file=sprintf("%s_%s_TAI_top_%d_genes_reduced_GSEA_results.Dmelanogaster.tsv", species_tag, sample, top_num), col.names = T, row.names = F, sep="\t")
}