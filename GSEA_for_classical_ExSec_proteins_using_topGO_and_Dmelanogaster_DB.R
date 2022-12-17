setwd("D:/Projects/Genes_level/Rhizocephala/MolSign_renewed/GSEA_Dmelanogaster/classES_MolSign/")
### Libraries ###
library(topGO)
library(dplyr)
library(rrvgo)
library(ggplot2)

### Input data ###
species_geneID2GO <- readMappings(file="../Pret_GeneOntologyUniverse_filtered.tsv")
species_geneNames <- names(species_geneID2GO)
species_tag <- "P.reticulata"

class_externa_molsign <- read.csv2("../../ExSec/Expr_analysis/ES_in_molsign/Pret_classical_ExSec_in_Externa_molsign.2TPM.txt", header = F, sep="\t")
class_growing_molsign <- read.csv2("../../ExSec/Expr_analysis/ES_in_molsign/Pret_classical_ExSec_in_Growing_molsign.2TPM.txt", header = F, sep="\t")
class_main_molsign <- read.csv2("../../ExSec/Expr_analysis/ES_in_molsign/Pret_classical_ExSec_in_Middle_molsign.2TPM.txt", header = F, sep="\t")
class_thoracic_molsign <- read.csv2("../../ExSec/Expr_analysis/ES_in_molsign/Pret_classical_ExSec_in_Terminal_molsign.2TPM.txt", header = F, sep="\t")

es_tag <- "classical"

### Analysis ###
input_list <- list("externa" = class_externa_molsign$V1,
                   "growing_part" = class_growing_molsign$V1, 
                   "main_part" = class_main_molsign$V1, 
                   "thoracic_part" = class_thoracic_molsign$V1)

for (set in 1:length(input_list)){
  set_all_genes <- factor(as.integer(species_geneNames %in% input_list[[set]]))
  names(set_all_genes) <- species_geneNames
  ## Biological Processes (BP) ##
  GOdata_set_BP <- new("topGOdata", ontology="BP", allGenes=set_all_genes,
                       annot=annFUN.gene2GO, gene2GO=species_geneID2GO)
  resultsFisher_set_BP <- runTest(GOdata_set_BP, algorithm="classic", statistic="fisher")
  resultsFisher_set_BP_df <- as.data.frame(score(resultsFisher_set_BP))
  colnames(resultsFisher_set_BP_df) <- "P-values"
  resultsFisher_set_BP_df_subset <- subset(resultsFisher_set_BP_df, resultsFisher_set_BP_df$`P-values` < 0.01)
  results_set_BP <- GenTable(GOdata_set_BP, classicFisher=resultsFisher_set_BP,
                             ranksOf="classicFisher", topNodes=length(resultsFisher_set_BP_df_subset$`P-values`))
  results_set_BP_short <- subset(results_set_BP, Significant >= 10)
  output_file_name <- sprintf("%s_%s_ES_in_%s_molsign_GOenrichment_BP_Fisher.min_10_genes.tsv", species_tag, es_tag, names(input_list[set]))
  write.table(results_set_BP_short, file=output_file_name, sep="\t", quote=F, col.names = T, row.names = F)
  if (nrow(results_set_BP_short) > 0){
    ## Reduced Terms ##
    set_simMatrix <- calculateSimMatrix(results_set_BP_short[["GO.ID"]], 
                                        orgdb = "org.Dm.eg.db", 
                                        ont="BP", method="Rel")
    set_classicFisher <- gsub("< ", "", results_set_BP_short[["classicFisher"]])
    set_scores <- setNames(-log10(as.numeric(as.character(set_classicFisher))), results_set_BP_short[["GO.ID"]])
    set_reducedTerms <- reduceSimMatrix(set_simMatrix, set_scores, threshold = 0.7, orgdb="org.Dm.eg.db")
    set_reducedTerms_selected <- select(set_reducedTerms, parent, score, parentTerm)
    # duplication removing: #
    set_reducedTerms_uniq <- set_reducedTerms_selected %>% 
      group_by(parentTerm) %>%
      filter(row_number() == 1)
    # Preparation for visualization #
    set_reducedTerms_uniq_df <- as.data.frame(set_reducedTerms_uniq)
    set_reducedTerms_uniq_df$GOid_and_desc <- paste(set_reducedTerms_uniq_df$parent, "|", set_reducedTerms_uniq_df$parentTerm)
    set_reducedTerms_uniq_df$score_signif <- signif(set_reducedTerms_uniq_df$score, digits = 3)
    # Visualization #
    pdf(file=sprintf("%s_%s_ES_in_%s_molsign_reduced_GSEA_results.Dmelanogaster_parent_GOterms.pdf", species_tag, es_tag, names(input_list[set])), width = 10, height = 12)
    result_plot <- set_reducedTerms_uniq_df %>%
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
      ggtitle(sprintf("Reduced GSEA results: \n %s potential %s ES \n in the molecular signature of %s", species_tag, es_tag, names(input_list[set]))) + labs(fill="Log-transformed scores")
    print(result_plot)
    dev.off()
    # Output writing #
    write.table(set_reducedTerms, file=sprintf("%s_%s_ES_in_%s_molsign_reduced_GSEA_results.Dmelanogaster.tsv", species_tag, es_tag, names(input_list[set])),
                sep="\t", quote = F, col.names = T, row.names = F)
  } else {
    print(sprintf("There is not enough enriched GO-terms for %s", names(input_list[set])))
  }
}