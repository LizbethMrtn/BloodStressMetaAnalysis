### EIF Meta-Analysis using Gemma Datasets that measure gene-expression (microarray, Agilent) in blood from mice exposed to stress (+/- other variables)
                            #GSE68076
                            #GSE72262
                            #GSE84185

## Following along with https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GemmaDatasets%20/2023_AligningLimmaResults_forMetaAnalysis.R 
#             (from https://github.com/hagenaue/BrainDataAlchemy)


# Aligning mouse datasets ##==============================READ IN THE FUNCTION===============================================================================
AligningMouseDatasets<-function(ListOfMouseDEResults){
  Mouse_MetaAnalysis_FoldChange_Dfs<-list()
  for(i in c(1:length(ListOfMouseDEResults))){
    Mouse_MetaAnalysis_FoldChange_Dfs[[i]]<-data.frame(Mouse_EntrezGene.ID=row.names(ListOfMouseDEResults[[i]][[1]]),ListOfMouseDEResults[[i]][[1]], stringsAsFactors=FALSE)
  }
  print("Mouse_MetaAnalysis_FoldChange_Dfs:")
  print(str(Mouse_MetaAnalysis_FoldChange_Dfs))
  Mouse_MetaAnalysis_FoldChanges<<-join_all(Mouse_MetaAnalysis_FoldChange_Dfs, by="Mouse_EntrezGene.ID", type="full")
  #This function could be join_all (if there are more than 2 datasets) or merge/merge_all (if the plyr package isn't working)
  print("Mouse_MetaAnalysis_FoldChanges:")
  print(str(Mouse_MetaAnalysis_FoldChanges))
  Mouse_MetaAnalysis_SV_Dfs<-list()
  for(i in c(1:length(ListOfMouseDEResults))){
    Mouse_MetaAnalysis_SV_Dfs[[i]]<-data.frame(Mouse_EntrezGene.ID=row.names(ListOfMouseDEResults[[i]][[4]]),ListOfMouseDEResults[[i]][[4]], stringsAsFactors=FALSE)
  }
  print("Mouse_MetaAnalysis_SV_Dfs:")
  print(str(Mouse_MetaAnalysis_SV_Dfs))
  Mouse_MetaAnalysis_SV<<-join_all(Mouse_MetaAnalysis_SV_Dfs, by="Mouse_EntrezGene.ID", type="full")
  #This function could be join_all (if there are more than 2 datasets) or merge/merge_all (if the plyr package isn't working)
  print("Mouse_MetaAnalysis_SV:")
  print(str(Mouse_MetaAnalysis_SV))
  rm(Mouse_MetaAnalysis_SV_Dfs, Mouse_MetaAnalysis_FoldChange_Dfs)
}



####################### USE FUNCTION TO ALIGN SPECIFIC DATASETS
# #Example Usage;
#         ListOfMouseDEResults<-list(DEResults_GSE85136, DEResults_GSE92718)
#         AligningMouseDatasets(ListOfMouseDEResults)

ListOfMouseDEResults<-list(DEResults_GSE68076, DEResults_GSE72262, DEResults_GSE84185)
AligningMouseDatasets(ListOfMouseDEResults)

########################Code for aligning the rat and mice results ===============================================================================================================
#We have the ortholog database that we downloaded from Jackson Labs on July 27, 2023
#This database was trimmed and formatted using the code "FormattingRatMouseOrthologDatabase_20230727.R"
MouseVsRat_NCBI_Entrez<-read.csv("MouseVsRat_NCBI_Entrez_JacksonLab_20230727.csv", header=TRUE, stringsAsFactors = FALSE, row.names=1, colClasses=c("character", "character", "character"))

Mouse_MetaAnalysis_FoldChanges_wOrthologs<-join(MouseVsRat_NCBI_Entrez, Mouse_MetaAnalysis_FoldChanges, by="Mouse_EntrezGene.ID", type="full")
str(Mouse_MetaAnalysis_FoldChanges_wOrthologs)
# 'data.frame':	21799 obs. of  6 variables:
#   $ Rat_EntrezGene.ID        : chr  "498097" "114521" "360576" "24628" ...
# $ Mouse_EntrezGene.ID      : chr  "68980" "21665" "237858" "18591" ...
# $ MouseVsRat_EntrezGene.ID : chr  "68980_498097" "21665_114521" "237858_360576" "18591_24628" ...
# $ GSE68076_SSvsNS          : num  0.0767 0.2614 NA 0.1963 NA ...
# $ GSE72262_StressVsNoStress: num  -0.2923 -0.3025 NA -0.0189 -0.2862 ...
# $ GSE84185_StressVsNoStress: num  0.1406 -0.0805 0.2544 -0.1339 NA ...

#If there aren't any rat datasets:
MetaAnalysis_FoldChanges<-Mouse_MetaAnalysis_FoldChanges_wOrthologs
str(MetaAnalysis_FoldChanges)


Mouse_MetaAnalysis_SV_wOrthologs<-join(MouseVsRat_NCBI_Entrez, Mouse_MetaAnalysis_SV, by="Mouse_EntrezGene.ID", type="full")
str(Mouse_MetaAnalysis_SV_wOrthologs)

#If there aren't any rat datasets:
MetaAnalysis_SV<-Mouse_MetaAnalysis_SV_wOrthologs

#For simplicity's sake, I'm going to replace that Mouse-Rat Entrez annotation
#Because it is missing entries for any genes in the datasets that *don't* have orthologs
MetaAnalysis_FoldChanges$MouseVsRat_EntrezGene.ID<-paste(MetaAnalysis_FoldChanges$Mouse_EntrezGene.ID, MetaAnalysis_FoldChanges$Rat_EntrezGene.ID, sep="_")
MetaAnalysis_SV$MouseVsRat_EntrezGene.ID<-paste(MetaAnalysis_SV$Mouse_EntrezGene.ID, MetaAnalysis_SV$Rat_EntrezGene.ID, sep="_")

#Comparing Log2FC across datasets
#Simple scatterplot... 
colnames(MetaAnalysis_FoldChanges)
#[1] "Rat_EntrezGene.ID"         "Mouse_EntrezGene.ID"       "MouseVsRat_EntrezGene.ID"  "GSE68076_SSvsNS"           "GSE72262_StressVsNoStress" "GSE84185_StressVsNoStress"
#plot(MetaAnalysis_FoldChangesINSERTVARIABLEHERE~MetaAnalysis_FoldChangesINSERTVARIABLEHERE)

plot(MetaAnalysis_FoldChanges$GSE68076_SSvsNS
       ~MetaAnalysis_FoldChanges$GSE72262_StressVsNoStress)
plot(MetaAnalysis_FoldChanges$GSE68076_SSvsNS
     ~MetaAnalysis_FoldChanges$GSE84185_StressVsNoStress)
plot(MetaAnalysis_FoldChanges$GSE72262_StressVsNoStress
     ~MetaAnalysis_FoldChanges$GSE84185_StressVsNoStress)

cor(as.matrix(MetaAnalysis_FoldChanges[,-c(1:3)]), use="pairwise.complete.obs", method="spearman")
#                                GSE68076_SSvsNS       GSE72262_StressVsNoStress    GSE84185_StressVsNoStress
# GSE68076_SSvsNS                1.00000000                0.03046632               -0.03229777
# GSE72262_StressVsNoStress      0.03046632                1.00000000               -0.15717105
# GSE84185_StressVsNoStress     -0.03229777               -0.15717105                1.00000000

#An illustration of the correlation matrix using a hierarchically clustered heatmap, although somewhat pathetic:
heatmap(cor(as.matrix(MetaAnalysis_FoldChanges[,-c(1:3)]), use="pairwise.complete.obs", method="spearman"))

