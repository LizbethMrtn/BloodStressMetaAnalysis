#Gene Ontology using Brain.GMT and fgsea
#https://www.biorxiv.org/content/10.1101/2024.04.05.588301v1

#################################################################################
setwd("/Users/flandree/Documents/2023 2024 Sabbatical/01 BrainAlchemy")
file_path = "/Users/flandree/Documents/2023 2024 Sabbatical/01 BrainAlchemy"
figure_path = "/Users/flandree/Documents/2023 2024 Sabbatical/01 BrainAlchemy"

#################################################################################
library(readr)
library(stats)
library(data.table)
library(ggplot2)
#install.packages("BiocParallel")
library(BiocParallel)
library(plyr)
library(dplyr)
library(stats)
library(tidyverse)
library(writexl)
#install.packages("conflicted")
library(conflicted)
#################################################################################

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("fgsea")
library(fgsea)
#################################################################################


DEResults<-read.csv("DEResults.csv", header=TRUE, stringsAsFactors = FALSE)
sum(is.na(DEResults))
#[1] 946

#Remove rows of DE results that are missing gene symbol annotation
#Replace $gene_symbol in the code with the column name containing gene symbols in your DE output
DEResults_noNA<-DEResults[is.na(DEResults$gene_symbol)==FALSE & is.na(DEResults$Log2FC_estimate)==FALSE,]
sum(is.na(DEResults_noNA))

str(DEResults_noNA)

#The analysis only works if there is one effect size (e.g., log2 fold change or Log2FC) per gene symbol.
#One way to deal with multiple effect sizes mapping to the same gene (e.g., multiple transcripts or probes) is to average them:
#Replace $Log2FC in the code with the column name containing effect sizes in your DE output
#Replace $gene_symbol in the code with the column name containing gene symbols in your DE output
DEResults_Log2FC_forGSEA<-tapply(X=DEResults_noNA$Log2FC_estimate, INDEX=DEResults_noNA$gene_symbol, FUN=mean)
names(DEResults_Log2FC_forGSEA)<-names(table(DEResults_noNA$gene_symbol))

#The effect sizes should be ordered from smallest to largest:
#Replace $Log2FC in the code with the column name containing effect sizes in your DE output
DEResults_Log2FC_forGSEA_Ranked<-DEResults_Log2FC_forGSEA[order(DEResults_Log2FC_forGSEA)]
# sum(is.na(DEResults_Log2FC_forGSEA_Ranked))
# is.numeric(DEResults_Log2FC_forGSEA_Ranked)
# str(DEResults_Log2FC_forGSEA_Ranked)
# head(DEResults_Log2FC_forGSEA_Ranked)

#DEResults_Log2FC_forGSEA_Ranked

#Read in Brain.GMT for your species of interest ###############--------------------------------........----------this takes a while<<<<<<
library(fgsea)
BrainGMT<-gmtPathways("BrainGMTv2_wGO_MouseOrthologs.gmt.txt")


#Run fast fGSEA on your ranked, averaged effect sizes:
GSEA_Results<-fgsea(BrainGMT, DEResults_Log2FC_forGSEA_Ranked, nperm=10000, minSize = 10, maxSize = 1000)
head(GSEA_Results)


#Pull out the names for the genes that are driving the enrichment of differential expression in each gene set:
GSEA_Results$leadingEdge<-vapply(GSEA_Results$leadingEdge, paste, collapse= ",", character(1L))

#Write out the results:
#write.csv(GSEA_Results, "GSEA_Results.csv")

#Collapse pathways#----------------------------
#duplicate table; filtering for padj < 0.05
GSEA_filtered2 <- GSEA_Results[GSEA_Results$padj<0.05, ]
str(GSEA_Results)

collapsedPathways <- collapsePathways(GSEA_filtered[order(pval)], BrainGMT, DEResults_Log2FC_forGSEA_Ranked, pval.threshold = 0.001, gseaParam=1) #takes forever to run
#unclear how they select the name of a cluster of pathways, difficult to identify relevance

str(collapsedPathways)
# List of 2
# $ mainPathways  : chr [1:143] "GOBP_PROTEIN_CONTAINING_COMPLEX_SUBUNIT_ORGANIZATION" "GOBP_INTRACELLULAR_TRANSPORT" "BLALOCK_ALZHEIMERS_DISEASE_UP" "GOCC_MITOCHONDRION" ...
# $ parentPathways: Named chr [1:280] NA NA NA NA ...
# ..- attr(*, "names")= chr [1:280] "GOBP_PROTEIN_CONTAINING_COMPLEX_SUBUNIT_ORGANIZATION" "GOBP_INTRACELLULAR_TRANSPORT" "BLALOCK_ALZHEIMERS_DISEASE_UP" "GOCC_MITOCHONDRION" ...

write.csv(collapsedPathways[[2]], "collapsedPathways.csv")

mainPathways <- GSEA_Results[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]

str(mainPathways)
mainPathways
write.csv(mainPathways, "mainPathways.csv")

# [1] "DESCARTES_MAIN_FETAL_ENS_GLIA"                                                                      
# [2] "HU_FETAL_RETINA_BLOOD"                                                                              
# [3] "GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX"                                                      
# [4] "MEISSNER_BRAIN_HCP_WITH_H3K4ME3_AND_H3K27ME3"                                                       
# [5] "GOBP_CELL_CELL_SIGNALING"                                                                           
# [6] "GOCC_INTRINSIC_COMPONENT_OF_PLASMA_MEMBRANE"                                                        
# [7] "GOBP_PROTEIN_MODIFICATION_BY_SMALL_PROTEIN_CONJUGATION_OR_REMOVAL"                                  
# [8] "BLALOCK_ALZHEIMERS_DISEASE_UP"                                                                      
# [9] "HP_ABNORMALITY_OF_THE_LIVER"                                                                        
# [10] "Gemma_GSE10415_7_Brain_Regions_in_20_Inbred_Strains_of_Laboratory_Mice_StrainOrLine_129S1.SvImJ_Up" 
# [11] "GOBP_CELLULAR_AMIDE_METABOLIC_PROCESS"                                                              
# [12] "LASTOWSKA_NEUROBLASTOMA_COPY_NUMBER_DN"                                                             
# [13] "GOCC_MITOCHONDRION"                                                                                 
# [14] "GOBP_LYMPHOCYTE_ACTIVATION"                                                                         
# [15] "GOBP_CHROMOSOME_ORGANIZATION"                                                                       
# [16] "GOBP_DNA_METABOLIC_PROCESS"                                                                         
# [17] "GOBP_REGULATION_OF_RESPONSE_TO_BIOTIC_STIMULUS"                                                     
# [18] "Gemma_GSE95449_Experience-dependent_Epigenomic_Reorganization_in_the_Hippocampus_timepoint_1.h_Down"
# [19] "GOCC_NUCLEAR_PROTEIN_CONTAINING_COMPLEX"                                                            
# [20] "JOHNSTONE_PARVB_TARGETS_3_DN"                                                                       
# [21] "GOBP_PEPTIDYL_LYSINE_MODIFICATION"                                                                  
# [22] "GOCC_RIBONUCLEOPROTEIN_COMPLEX"                                                                     
# [23] "GOMF_RNA_BINDING"                                                                                   
# [24] "FAN_EMBRYONIC_CTX_BRAIN_B_CELL"                                                                     
# [25] "GOBP_RNA_PROCESSING"                                                                                
# [26] "HP_ECZEMA"                                                                                          
# [27] "GOBP_DNA_REPLICATION" 




