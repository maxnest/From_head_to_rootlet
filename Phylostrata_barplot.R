setwd("D:/Projects/Genes_level/Rhizocephala/MolSign_renewed/Phylostratr/")
### Libraries ###
library(dplyr)
library(ggplot2)
library(viridis)
library(stringr)
library(reshape)
theme_set(theme_minimal())

### Input ###
phylostratr_tab <- read.csv2("Preticulata_phylostratr_table.tsv", header = F, sep="\t")

### Processing ###
colnames(phylostratr_tab) <- c("Phylostrata", "GeneCount")
phylostratr_tab$Phylostrata
phylostratr_tab$Phylostrata <- c("1:cellular organisms", "2:Eukaryota", "3:Opisthokonta", "4:Metazoa",                
                                 "5:Eumetazoa", "6:Bilateria", "7:Protostomia", "8:Ecdysozoa",              
                                 "9:Panarthropoda", "10:Arthropoda", "11:Mandibulata", "12:Pancrustacea",           
                                 "13:Crustacea", "14:Multicrustacea", "15:Hexanauplia", "16:Cirripedia", 
                                 "17:P.reticulata")
phylostratr_tab <- mutate(phylostratr_tab, Percent = (GeneCount / 12620) * 100)

ordered_phylostrata <- c("1:cellular organisms", "2:Eukaryota", "3:Opisthokonta", "4:Metazoa",                
                         "5:Eumetazoa", "6:Bilateria", "7:Protostomia", "8:Ecdysozoa",              
                         "9:Panarthropoda", "10:Arthropoda", "11:Mandibulata", "12:Pancrustacea",           
                         "13:Crustacea", "14:Multicrustacea", "15:Hexanauplia", "16:Cirripedia", 
                         "17:P.reticulata")

phylostratr_tab_percent <- select(phylostratr_tab, Phylostrata, Percent)
str(phylostratr_tab_percent)

phylostratr_tab_percent_melt <- melt.data.frame(phylostratr_tab_percent, id.vars = "Phylostrata")
phylostratr_tab_percent_melt$Phylostrata <- factor(phylostratr_tab_percent_melt$Phylostrata, levels=ordered_phylostrata)

str(phylostratr_tab_percent_melt)

### Visual ###
pdf("Preticulata_100aa_phylosummary.percent.barplot.pdf", width = 12)
species_barplot <- ggplot(phylostratr_tab_percent_melt, aes(x=Phylostrata, y=value, fill=Phylostrata)) + 
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_viridis(discrete = TRUE, direction=-1) + 
  scale_color_manual(values=c("black", "white")) +
  ggtitle(expression(paste(italic("P. reticulata"), " ", 'phylostratigraphy results summary'))) + 
  xlab("Phylostrata") + ylab("Percent") + 
  theme(legend.position="none",
        plot.title = element_text(size=15, face="bold", hjust = 0.5), 
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.5),
        axis.text=element_text(size=12),
        axis.title = element_text(size=15, face="bold")
  )
print(species_barplot)
dev.off()

### Output ###
write.table(phylostratr_tab, file="Preticulata_phylostratr_table.with_percents.tsv", col.names = T, row.names = F, sep="\t")
