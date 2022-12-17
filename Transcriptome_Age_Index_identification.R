setwd("/media/sf_Space/Crustacea_renewed/myTAI/")
### Libraries ###
library(edgeR)
library(myTAI)
library(dplyr)

### Input data ###
phylostratr <- read.csv2("Preticulata_phylostratr_results.geneIDs.tsv", sep="\t", header = T)
expr <- read.csv2("Preticulata_ref.genes_level.after_filters.averaged_TPMs.tsv", sep="\t", header = T)
species_tag <- "P.reticulata"

### Processing ###
expr_mut <- mutate_all(expr[, 1:length(colnames(expr))], function(x) as.numeric(as.character(x)))
expr_mut$GeneIDs <- expr$GeneIDs
colnames(expr_mut)
expr_mut <- select(expr_mut, "GeneIDs", "Externa", "Growing", "Middle", "Terminal")
colnames(expr_mut) <- c("GeneIDs", "Externa", "Growing_part", "Main_part", "Thoracic_part")

expr_mut <- subset(expr_mut, Externa >= 2 | Growing_part >= 2 | Main_part >= 2 | Thoracic_part >= 2)

colnames(phylostratr) <- c("GeneIDs", "MRCA", "Phylostratum", "MRCA_name")

phylostratr_subset <- subset(phylostratr, GeneIDs %in% expr_mut$GeneIDs)

### Transformation ###
## PhyloMap ##
phylomap <- select(phylostratr_subset, "Phylostratum", "GeneIDs")
## Mapping ##
phyloexp <- MatchMap(phylomap, expr_mut)
## Log-transformation ##
phyloexp_tf <- tf(phyloexp, function(x) log2(x+1))
### Analysis ###
## Расчет TAI (Transcriptome Age Index) ##
# The lower the TAI value the older the mean transcriptome age and the higher the TAI value the younger the mean transcriptome age
pdf(sprintf("%s_filtered_logTPMs_TAI.pdf", species_tag), width = 10)
PlotSignature(ExpressionSet = phyloexp_tf, measure = "TAI",
              TestStatistic = "FlatLineTest", 
              xlab=expression(paste(italic("P. reticulata"), " ", "female body parts")),
              ylab="Transcriptome Age Index")
dev.off()

# Due to the nature of the arithmetic mean, this value does not represent the true origin of individual genes, 
# and thus the TAI measure is only helpful to screen for stages that express (on average) older or younger genes. 
# Subsequent analyses such as mean expression of age categories, relative expression levels, 
# and gene expression level distributions for each age category will then reveal which exact genes or age categories generate the overall TAI value.

# Category-centered visualization of phylostrates
# specific expression level distributions (log-scale)
pdf(sprintf("%s_filtered_logTPMs_CategoryExpr.pdf", species_tag), width = 10)
PlotCategoryExpr(ExpressionSet = phyloexp_tf, legendName = "PS",
                 test.stat = TRUE, type="category-centered", 
                 distr.type = "boxplot",
                 log.expr=FALSE)
dev.off()

### Contibution ###
# How the final (global) TAI profile emerges from the cumulative TAI distribution of each Phylostratum?
pdf(sprintf("%s_filtered_logTPMs_Phylostrata_contribution.pdf", species_tag), width = 10)
PlotContribution(phyloexp_tf, legendName = "PS",
                 xlab=expression(paste(italic("P. reticulata"), " ", "female body parts")), 
                 ylab = "Transcriptome Age Index")
dev.off()

### Quantife mean expression of individual gene age categories ###
# PS1-13: 1) Cellular organism - 13) Crustacea
# PS13-*** (species-specific)

pdf(sprintf("%s_filtered_logTPMs_mean_expression.pdf", species_tag), width = 10)
PlotMeans(phyloexp_tf, Groups=list(c(1:13), c(14:max(phyloexp_tf$Phylostratum))),
          legendName = "PS", xlab=expression(paste(italic("P. reticulata"), " ", "female body parts")),
          adjust.range = TRUE)
dev.off()

### Relative expression (RE) of individual gene age categories ###
pdf(sprintf("%s_filtered_logTPMs_RE.pdf", species_tag), width = 10)
PlotRE(phyloexp_tf, Groups=list(c(1:13), c(14:max(phyloexp_tf$Phylostratum))),
       legendName = "PS", xlab=expression(paste(italic("P. reticulata"), " ", "female body parts")),
       adjust.range = TRUE)
dev.off()

### RE comparison ###
pdf(sprintf("%s_filtered_logTPMs_RE_comparison.pdf", species_tag), width = 10)
PlotBarRE(phyloexp_tf, Groups=list(c(1:13), c(14:max(phyloexp_tf$Phylostratum))),
          xlab=expression(paste(italic("P. reticulata"), " ", "female body parts")), 
          ylab="Mean Relative Expression",
          cex=1.5)
dev.off()

### Relative gene frequency distribution ###
pdf(sprintf("%s_filtered_logTPMs_relative_gene_frequency_distribution.pdf", species_tag), width = 10)
PlotDistribution(phyloexp_tf, as.ratio = TRUE, xlab = "Phylostrata")
dev.off()

### Outputs ###
# TAI #
write.table(as.data.frame(TAI(phyloexp_tf)), file=sprintf("%s_filtered_logTPMs.TAI.tsv", species_tag), sep="\t", col.names = T, row.names = T)
# Phylostrata contribution #
write.table(as.data.frame(pTAI(phyloexp_tf)), file=sprintf("%s_filtered_logTPMs.Phylostratum_to_TAI_contribution.tsv", species_tag), sep="\t", col.names = T, row.names = T)
# Partial contribution #
write.table(as.data.frame(pStrata(phyloexp_tf)), file=sprintf("%s_filtered_logTPMs.Phylostratum_to_TAI_partial_contribution.tsv", species_tag), sep = "\t", col.names = T, row.names = T)
# Mean expression #
write.table(as.data.frame(age.apply(ExpressionSet = phyloexp_tf, FUN = colMeans)), file=sprintf("%s_filtered_logTPMs.Phylostratum_mean_expression.tsv", species_tag), sep="\t", col.names = T, row.names = T)
# Relative expression #
write.table(as.data.frame(REMatrix(phyloexp_tf)), file=sprintf("%s_filtered_logTPMs.Phylostratum_relative_expression.tsv", species_tag), sep="\t", col.names = T, row.names = T)
# Genes contribution #
write.table(as.data.frame(pMatrix(phyloexp_tf)), file=sprintf("%s_filtered_logTPMs.Genes_to_TAI_contribution.tsv", species_tag), sep="\t", col.names = T, row.names = T)