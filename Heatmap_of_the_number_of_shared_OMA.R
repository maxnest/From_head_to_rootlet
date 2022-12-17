setwd("D:/Projects/Genes_level/Rhizocephala/MolSign_renewed/OMA_renewed/")
### Libraries ###
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

### Input data ###
SpeciesOverlap <- read.csv2("Crustacea.num_shared_OMAs.tsv", header=T, sep="\t", row.names = 1)
my_palette <- colorRampPalette(c("white", "light blue", "blue", "orange", "red"))(n=100)

## Processing #
names(SpeciesOverlap)
colnames(SpeciesOverlap) <- c("A.amphitrite", "A.nasatum", "A.vulgare", "D.magna" , "D.pulex", 
                              "P.reticulata", "P.vannamei", "P.trituberculatus", "T.californicus")
rownames(SpeciesOverlap) <- c("A.amphitrite", "A.nasatum", "A.vulgare", "D.magna" , "D.pulex", 
                              "P.reticulata", "P.vannamei", "P.trituberculatus", "T.californicus")

### Output files creating ###
pdf("Crustacea_100aa.OMA_species_overlaps.19x19.pdf")
pheatmap(SpeciesOverlap, cluster_cols = T, cluster_rows = T,  clustering_distance_cols = "correlation",
         clustering_distance_rows = "correlation", clustering_method = "complete", col=my_palette,
         scale="none", key=T, symkey = F, density.info = "none", trace="none", cexCol=1, cexRow=1,
         cellwidth = 19, cellheight = 19)
dev.off()
