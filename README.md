# From head to rootlet: comparative transcriptomic analysis of a rhizocephalan barnacle *Peltogaster reticulata* (Crustacea: Rhizocephala)

## Maksim A. Nesterenko1,2* and Aleksei A. Miroliubov2
### 1 Department of Invertebrate Zoology, St Petersburg State University, St Petersburg 199034, Russia.
### 2 Laboratory of parasitic worms and protists, Zoological Institute, Russian Academy of Sciences, St Petersburg 199034, Russia.
#### *Correspondence: maxnest.research@gmail.com

## Abstract
### Background: 
Rhizocephalan barnacles stand out in the diverse world of metazoan parasites. The body of a rhizocephalan female is modified beyond revealing any recognizable morphological features, consisting of the interna, a system of rootlets, and the externa, a sac-like reproductive body. Moreover, rhizocephalans have an outstanding ability to control their hosts, literally turning them into “zombies”. Despite all these amazing traits, there are no genomic or transcriptomic data about any Rhizocephala. 
### Methods: 
We collected transcriptomes from four body parts of an adult female rhizocephalan *Peltogaster reticulata* : the externa, and the main, growing, and thoracic parts of the interna. We used all prepared data for the de novo assembly of the reference transcriptome. Next, a set of encoded proteins was determined, the expression levels of protein-coding genes in different parts of the parasite’s body were calculated and lists of enriched bioprocesses were identified. We also in silico identified and analyzed sets of potential excretory / secretory proteins. Finally, we applied phylostratigraphy and evolutionary transcriptomics approaches to our data. 
### Results: 
The assembled reference transcriptome included transcripts of 12,620 protein-coding genes and was the first for any rhizocephalan. Based on the results obtained, the spatial heterogeneity of protein-coding gene expression in different regions of the  adult female body of *P. reticulata* was established. The results of both transcriptomic analysis and histological studies indicated the presence of germ-like cells in the lumen of the interna. The potential molecular basis of the interaction between the nervous system of the host and the parasite's interna was also determined. Given the prolonged expression of development-associated genes, we suggest that rhizocephalans “got stuck in their metamorphosis”, even at the reproductive stage.
### Conclusions: 
The results of the first comparative transcriptomic analysis for Rhizocephala not only clarified but also expanded the existing ideas about the biology of these extraordinary parasites.
### Keywords:
Rhizocephala, parasitic barnacles, evolutionary transcriptomics, host manipulation, coloniality

The full text of the article is published on the platform [F1000Research](https://f1000research.com/articles/11-583/v1).

## Below are brief descriptions of the scripts used in the study:
### Preparation of a set of protein-coding genes with a noticeable level of expression
```
1. From_transcripts_expr_levels_to_genes_expr_levels.R
2. Set_of_protein-coding_genes_with_high_expression.py
```
### Selection of orthologs groups detected using the [OMA standalone](https://omabrowser.org/standalone/) for subsequent visualization in the form of a heatmap
```
1. OMAstandalone_num_of_shared_OMA.py
2. Heatmap_of_the_number_of_shared_OMA.R
```
### Analysis of differential gene expression using the [RNentropy](https://cran.r-project.org/web/packages/RNentropy/index.html) R package
```
expressed_genes_identification_using_RNentropy.R
```
### Identification of potential excretory/secretory proteins and genes encoding these products
```
1. SignalP_parsing.py
2. Renaming_for_secretomeP.py
3. SecretomeP_results_parsing.py
4. TMHMM_parser.py
5. ExSec_proteinsIDs_with_noticeable_expr.R
6. GenesIDs_encoding_ExSec_proteins_included_in_molsignatures.R
7. GenesIDs_encoding_ExSec_proteins_with_overexpr.R
```
### Gene Set Enrichment Analysis (GSEA) and visualization of significant Gene Ontology (GO) terms as wordclouds  
```
1. From_eggNOG_to_topGO_input.R
2. GSEA_for_molsignatures_using_topGO_and_Dmelanogaster_DB.R
3. GSEA_for_overexpr_genes_using_topGO_and_Dmelanogaster_DB.R
4. GSEA_for_classical_ExSec_proteins_using_topGO_and_Dmelanogaster_DB.R
5. GOterms_wordclouds_for_molsignatures.R
6. GOterms_wordclouds_for_overexpr_genes.R
7. GOterms_wordclouds_for_classExSec_proteins.R
```
### Multidimensional scaling (MDS) of the considered samples
```
1. Molsignatures_pres_abs_matrix_for_MDS.R
2. Overexpr_genes_pres_abs_matrix_for_MDS.R
3. Silhouette_Plot_optimal_num_of_clusters_for_MDS.R
4. Molsignatures_MDS.R
5. Overexpressed_genes_MDS.R
```
### Phylostratigraphy using the [Phylostratr](https://github.com/arendsee/phylostratr) R framework
```
1. Phylostratr_first_step.R
2. run_BLASTP_for_phylostratr.sh
3. Phylostratr_second_step.R
4. Phylostrata_barplot.R
5. Gene_set_phylostratigraphy_affiliation.R
```
### Calculation of Transcriptome Age Indices (TAI) using the [myTAI](https://github.com/drostlab/myTAI) R package
```
1. Transcriptome_Age_Index_identification.R
2. GSEA_for_top_500_genes_with_highest_contribution_to_TAI.R
```
