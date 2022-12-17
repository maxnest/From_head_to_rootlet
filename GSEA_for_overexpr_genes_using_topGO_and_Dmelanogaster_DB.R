setwd("D:/Projects/Genes_level/Rhizocephala/MolSign_renewed/GSEA_Dmelanogaster/Overexpr/")
### Libraries ###
library(topGO)
library(dplyr)
library(rrvgo)
library(ggplot2)

options(ggrepel.max.overlaps = Inf)

### Input data ###
species_geneID2GO <- readMappings(file="../Pret_GeneOntologyUniverse_filtered.tsv")
species_geneNames <- names(species_geneID2GO)
species_tag <- "P.reticulata"
molsign_tab <- read.csv2("../../Expr_analysis/DiffExpr/Pret_RNentropy_significant_results.tsv", header = T, sep="\t")

stages <- colnames(molsign_tab[, 4:ncol(molsign_tab)])

### Analysis ###
for (stage in stages){
  stage_molsign <- select(subset(molsign_tab, molsign_tab[stage] == 1), GeneIDs)
  stage_all_genes <- factor(as.integer(species_geneNames %in% stage_molsign$GeneIDs))
  names(stage_all_genes) <- species_geneNames
  ## Biological Processes (BP) ##
  GOdata_stage_BP <- new("topGOdata", ontology="BP", allGenes=stage_all_genes,
                         annot=annFUN.gene2GO, gene2GO=species_geneID2GO)
  resultsFisher_stage_BP <- runTest(GOdata_stage_BP, algorithm="classic", statistic="fisher")
  resultsFisher_stage_BP_df <- as.data.frame(score(resultsFisher_stage_BP))
  colnames(resultsFisher_stage_BP_df) <- "P-values"
  resultsFisher_stage_BP_df_subset <- subset(resultsFisher_stage_BP_df, resultsFisher_stage_BP_df$`P-values` < 0.01)
  results_stage_BP <- GenTable(GOdata_stage_BP, classicFisher=resultsFisher_stage_BP,
                               ranksOf="classicFisher", topNodes=length(resultsFisher_stage_BP_df_subset$`P-values`))
  results_stage_BP_short <- subset(results_stage_BP, Significant >= 10)
  output_file_name <- sprintf("%s_%s_overexpr_genes_GOenrichment_BP_Fisher.min_10_genes.tsv", species_tag, stage)
  write.table(results_stage_BP_short, file=output_file_name, sep="\t", quote=F, col.names = T, row.names = F)
  ## Reduced Terms ##
  stage_simMatrix <- calculateSimMatrix(results_stage_BP_short[["GO.ID"]], 
                                        orgdb = "org.Dm.eg.db", 
                                        ont="BP", method="Rel")
  stage_classicFisher <- gsub("< ", "", results_stage_BP_short[["classicFisher"]])
  stage_scores <- setNames(-log10(as.numeric(as.character(stage_classicFisher))), results_stage_BP_short[["GO.ID"]])
  stage_reducedTerms <- reduceSimMatrix(stage_simMatrix, stage_scores, threshold = 0.7, orgdb="org.Dm.eg.db")
  stage_reducedTerms_selected <- select(stage_reducedTerms, parent, score, parentTerm)
  # duplication removing: #
  stage_reducedTerms_uniq <- stage_reducedTerms_selected %>% 
    group_by(parentTerm) %>%
    filter(row_number() == 1)
  # Preparation for visualization #
  stage_reducedTerms_uniq_df <- as.data.frame(stage_reducedTerms_uniq)
  stage_reducedTerms_uniq_df$GOid_and_desc <- paste(stage_reducedTerms_uniq_df$parent, "|", stage_reducedTerms_uniq_df$parentTerm)
  stage_reducedTerms_uniq_df$score_signif <- signif(stage_reducedTerms_uniq_df$score, digits = 3)
  # Visualization #
  # Plot #
  pdf(file=sprintf("%s_%s_overexpr_genes_reduced_GSEA_results.Dmelanogaster_parent_GOterms.pdf", species_tag, stage), width = 10, height = 12)
  result_plot <- stage_reducedTerms_uniq_df %>%
    arrange(score) %>%
    mutate(GOid_and_desc=factor(GOid_and_desc, levels=GOid_and_desc)) %>% 
    ggplot(aes(x=score, y=GOid_and_desc, fill=as.factor(score_signif))) + 
    geom_point(alpha=0.9, shape=21, size=5, color="black") + theme_bw() + 
    # legend.position = "right", axis.text=element_text(size=12, face="bold")
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5, size = 15, face="bold"), 
          axis.title.x = element_text(size=13, face="bold"), 
          axis.title.y = element_text(size=13, face="bold")) + 
    xlab("-log10(Fisher's Test p-values)") + ylab("Parental Gene Ontology terms (Fly)") + 
    ggtitle(sprintf("Reduced GSEA results: \n %s %s molecular signature \n (over-expressed genes)", species_tag, stage)) + labs(fill="Log-transformed scores")
  print(result_plot)
  dev.off()
  # ScatterPlot #
  pdf(file=sprintf("%s_%s_overexpr_genes_reduced_GSEA_results.Dmelanogaster.scatterPlot.pdf", species_tag, stage), width = 15, height = 15)
  scatter <- scatterPlot(stage_simMatrix, stage_reducedTerms, size="score", addLabel = TRUE, labelSize = 5)
  print(scatter)
  dev.off()
  # Output writing #
  write.table(stage_reducedTerms, file=sprintf("%s_%s_overexpr_genes_reduced_GSEA_results.Dmelanogaster.tsv", species_tag, stage),
              sep="\t", quote = F, col.names = T, row.names = F)
}