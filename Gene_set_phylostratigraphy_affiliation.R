setwd("D:/Projects/Genes_level/Rhizocephala/MolSign_renewed/Phylostratigraphy_affiliation/Visual/")
### Library ###
library(dplyr)
library(ggplot2)
library(viridis)
library(reshape)
theme_set(theme_minimal())

### Input ###
## Phylostratr ##
phylostr_ref_set <- read.csv2("../../Phylostratr/Preticulata_phylostratr_table.with_percents.tsv", header = T, sep="\t")
## Expression ##
common_expr <- read.csv2("../Common_expr/Pret_ref_common_expr_2TPM.phylostrata_summary_table.tsv", header = T, sep="\t")
externa_over <- read.csv2("../OverExpr/Pret_Externa_RNentropy_overexp.phylostrata_summary_table.tsv", header = T, sep="\t")
growing_over <- read.csv2("../OverExpr/Pret_Growing_part_RNentropy_overexp.phylostrata_summary_table.tsv", header = T, sep="\t")
main_over <- read.csv2("../OverExpr/Pret_Main_part_RNentropy_overexp.phylostrata_summary_table.tsv", header = T, sep="\t")
thoracic_over <- read.csv2("../OverExpr/Pret_Thoracic_part_RNentropy_overexp.phylostrata_summary_table.tsv", header = T, sep="\t")
## ES ##
class_es <- read.csv2("../ExSec/Pret_class_ExSec_with_noticeable_expr.phylostrata_summary_table.tsv", header = T, sep="\t")
nonclass_es <- read.csv2("../ExSec/Pret_nonclass_ExSec_with_noticeable_expr.phylostrata_summary_table.tsv", header = T, sep="\t")

### Processing ###
phylostr_ref_set$Phylostrata <- common_expr$Phylostratum

set_merged <- data.frame("Phylostrata"=common_expr$Phylostratum, "Reference"=phylostr_ref_set$Percent, 
                         "Common_expr"=common_expr$Percent, "Externa_overexpr"=externa_over$Percent, 
                         "Growing_part_overexpr"=growing_over$Percent, "Main_part_overexpr"=main_over$Percent, "Thoracic_part_overexpr"=thoracic_over$Percent, 
                         "Classical_ESP"=class_es$Percent, "Nonclassical_ESP" = nonclass_es$Percent)
set_merged$Phylostrata
ordered_phylostrata <- c("1:cellular organisms", "2:Eukaryota", "3:Opisthokonta", "4:Metazoa", "5:Eumetazoa", 
                         "6:Bilateria", "7:Protostomia", "8:Ecdysozoa", "9:Panarthropoda", "10:Arthropoda", "11:Mandibulata", 
                         "12:Pancrustacea", "13:Crustacea", "14:Multicrustacea", "15:Hexanauplia", "16:Cirripedia", "17:P.reticulata")
colnames(set_merged)
ordered_sets <- c("Reference", "Common_expr", "Externa_overexpr", "Growing_part_overexpr", "Main_part_overexpr", "Thoracic_part_overexpr", "Classical_ESP", "Nonclassical_ESP")

set_merged_matrix <- mutate_all(set_merged, function(x) as.numeric(as.character(x)))
set_merged_matrix$Phylostrata <- set_merged$Phylostrata
set_merged_matrix$Phylostrata <- c("1:cellular organisms", "2:Eukaryota", "3:Opisthokonta", "4:Metazoa", "5:Eumetazoa", 
                                   "6:Bilateria", "7:Protostomia", "8:Ecdysozoa", "9:Panarthropoda", "10:Arthropoda", "11:Mandibulata", 
                                   "12:Pancrustacea", "13:Crustacea", "14:Multicrustacea", "15:Hexanauplia", "16:Cirripedia", "17:P.reticulata")

set_merged_matrix_melt <- melt.data.frame(set_merged_matrix, id.vars = "Phylostrata")
set_merged_matrix_melt$Phylostrata <- factor(set_merged_matrix_melt$Phylostrata, levels=rev(ordered_phylostrata))
set_merged_matrix_melt$variable <- factor(set_merged_matrix_melt$variable, levels=ordered_sets)

### Visual ###
pdf("Preticulata_ref_gene_set_plus_expression_and_ESP.phylostratigraphic_aff.percent_barplot.pdf", width = 7)
set_merged_barplot <- ggplot(set_merged_matrix_melt, aes(x=variable, y=value, fill=Phylostrata)) + 
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_viridis(discrete = TRUE, direction=1) + 
  scale_color_manual(values=c("black", "white")) +
  xlab("Gene sets") + ylab("Percent") + 
  theme(legend.key.size = unit(0.5, 'cm'), 
        legend.title = element_text(size=10, face="bold"), 
        legend.text = element_text(size=10),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 0.8, hjust=0.5, face="bold"),
        axis.text=element_text(size=10),
        axis.title = element_text(size=15, face="bold"))
print(set_merged_barplot)
dev.off()