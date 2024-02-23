#Simple Meta-Analysis Version
### EIF Meta-Analysis using Gemma Datasets that measure gene-expression (microarray, Agilent) in blood from mice exposed to stress (+/- other variables)
#GSE68076
#GSE72262
#GSE84185

## Following along with https://github.com/hagenaue/BrainDataAlchemy/blob/main/MetaAnalysis_GemmaDatasets%20/2023_MetaAnalysisCode.R
#             (from https://github.com/hagenaue/BrainDataAlchemy)
#8) Run a meta-analysis using all of the effect sizes for each gene that has data from at least 2 studies.
#Since I know that the differential expression results from the same study (dataset) are artificially correlated, 
##            I would actually prefer that there are results from more than one dataset.
#How many genes satisfy this criteria?

#This code caculates the number of NAs in each row:
MetaAnalysis_FoldChanges_NAsPerRow<-apply(MetaAnalysis_FoldChanges[,-c(1:3)], 1, function(y) sum(is.na(y)))

#I'm going to make a histogram and a table of the results because I'm curious to see how they are distributed
hist(MetaAnalysis_FoldChanges_NAsPerRow)
table(MetaAnalysis_FoldChanges_NAsPerRow)
# MetaAnalysis_FoldChanges_NAsPerRow
# 0    1    2    3 
# 9221 8004 4024  550 

#Let's try running a meta-analysis using genes that were found in at least 2 sets of differential expression results
#Since there are 3 sets of differential expression results, that means that the genes that we are including need to have 1 or fewer NAs in their results
#I set this conservatively, because there are so few studies in this meta-analysis.
#2 NA is too many

#Installing and loading relevant code packages:
#install.packages("metafor")
library(metafor)
library(plyr)

## READ IN THE FUNCTION ##====================================================
RunBasicMetaAnalysis<-function(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV){
  
  MetaAnalysis_FoldChanges_NAsPerRow<-apply(MetaAnalysis_FoldChanges[,-c(1:3)], 1, function(y) sum(is.na(y)))
  
  print("Table of # of NAs per Row (Gene):")
  print(table(MetaAnalysis_FoldChanges_NAsPerRow))
  
  MetaAnalysis_FoldChanges_ForMeta<<-MetaAnalysis_FoldChanges[MetaAnalysis_FoldChanges_NAsPerRow<CutOffForNAs,]
  MetaAnalysis_SV_ForMeta<<-MetaAnalysis_SV[MetaAnalysis_FoldChanges_NAsPerRow<CutOffForNAs,]
  
  print("MetaAnalysis_FoldChanges_ForMeta:")
  print(str(MetaAnalysis_FoldChanges_ForMeta))
  
  #I'm going to make an empty matrix to store the results of my meta-analysis:
  metaOutput<-matrix(NA, nrow(MetaAnalysis_FoldChanges_ForMeta), 6)
  
  #And then run a loop that run's a meta-analysis on the differential expression results (columns 2-10) for each gene (row):
  for(i in c(1:nrow(MetaAnalysis_FoldChanges_ForMeta))){
    
    effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[i,-c(1:3)])
    var<-as.numeric(MetaAnalysis_SV_ForMeta[i,-c(1:3)])
    
    #I added a function tryCatch that double-checks that the meta-analysis function (rma) doesn't produce errors (which breaks the loop):
    skip_to_next <- FALSE
    tryCatch(TempMeta<-rma(effect, var), error = function(e) {skip_to_next <<- TRUE})
    
    if(skip_to_next){}else{
      TempMeta<-rma(effect, var)
      metaOutput[i, 1]<-TempMeta$b #gives estimate Log2FC
      metaOutput[i, 2]<-TempMeta$se #gives standard error
      metaOutput[i, 3]<-TempMeta$pval #gives pval
      metaOutput[i, 4]<-TempMeta$ci.lb #gives confidence interval lower bound
      metaOutput[i, 5]<-TempMeta$ci.ub #gives confidence interval upper bound
      metaOutput[i, 6]<-NumberOfComparisons-sum(is.na(effect))#Number of comparisons with data
      rm(TempMeta)
    }
    rm(effect, var)
  }
  
  colnames(metaOutput)<-c("Log2FC_estimate", "SE", "pval", "CI_lb", "CI_ub", "Number_Of_Comparisons")
  row.names(metaOutput)<-MetaAnalysis_FoldChanges_ForMeta[,3]
  
  metaOutput<<-metaOutput
  MetaAnalysis_Annotation<<-MetaAnalysis_FoldChanges_ForMeta[,c(1:3)]
  return(metaOutput)
  return(MetaAnalysis_Annotation)
  
  print("metaOutput:")
  print(str(metaOutput))
  
  print("Top of metaOutput:")
  print(head(metaOutput))
  
  print("Bottom of metaOutput")
  print(tail(metaOutput))
  
}

#### ============================================================END OF FUNCTION ============================================================

### ============================================================USE THE FUNCTION============================================================
#Example Usage:
NumberOfComparisons=3
CutOffForNAs=1
#It is the number of group comparisons, so in your case it will be three 
# (study 1: stress vs. no stress, study 2: stress vs. no stress, study 3: stress vs. no stress).
# If you hypothetically had more than one group of interest in a study, it might be more than three 
# (e.g., Study 1: 10 days of stress vs. no stress, Study 1: 5 days of stress vs. no stress...)
#For number of NAs, 2 will be too many for you (your code will crash if you only have one study of results available for a gene...), 
#but you could also say 1 NA is too many and require that a gene have results in all three studies to be included in your meta-analysis.

metaOutput<-RunBasicMetaAnalysis(NumberOfComparisons, CutOffForNAs, MetaAnalysis_FoldChanges, MetaAnalysis_SV)
#Note: this function can take a while to run, especially if you have a lot of data  
#Plug in your computer, take a break, grab some coffee...


# [1] "Table of # of NAs per Row (Gene):"
# MetaAnalysis_FoldChanges_NAsPerRow
# 0    1    2    3 
# 9221 8004 4024  550 
# [1] "MetaAnalysis_FoldChanges_ForMeta:"
# 'data.frame':	9221 obs. of  6 variables:
#   $ Rat_EntrezGene.ID        : chr  "498097" "114521" "24628" "312258" ...
# $ Mouse_EntrezGene.ID      : chr  "68980" "21665" "18591" "57869" ...
# $ MouseVsRat_EntrezGene.ID : chr  "68980_498097" "21665_114521" "18591_24628" "57869_312258" ...
# $ GSE68076_SSvsNS          : num  0.0767 0.2614 0.1963 -0.4161 -0.1684 ...
# $ GSE72262_StressVsNoStress: num  -0.2923 -0.3025 -0.0189 -0.3469 -0.1316 ...
# $ GSE84185_StressVsNoStress: num  0.1406 -0.0805 -0.1339 0.0334 -0.0277 ...
# NULL
# There were 12 warnings (use warnings() to see them)
# > 

warnings()
# Warning messages:
#   1: Fisher scoring algorithm may have gotten stuck at a local maximum.
# Setting tau^2 = 0. Check the profile likelihood plot with profile().
# 2: Fisher scoring algorithm may have gotten stuck at a local maximum.
# Setting tau^2 = 0. Check the profile likelihood plot with profile().
# 3: Fisher scoring algorithm may have gotten stuck at a local maximum.
# Setting tau^2 = 0. Check the profile likelihood plot with profile().
# 4: Fisher scoring algorithm may have gotten stuck at a local maximum.
# Setting tau^2 = 0. Check the profile likelihood plot with profile().
# 5: Fisher scoring algorithm may have gotten stuck at a local maximum.
# Setting tau^2 = 0. Check the profile likelihood plot with profile().
# 6: Fisher scoring algorithm may have gotten stuck at a local maximum.
# Setting tau^2 = 0. Check the profile likelihood plot with profile().
# 7: Fisher scoring algorithm may have gotten stuck at a local maximum.
# Setting tau^2 = 0. Check the profile likelihood plot with profile().
# 8: Fisher scoring algorithm may have gotten stuck at a local maximum.
# Setting tau^2 = 0. Check the profile likelihood plot with profile().
# 9: Fisher scoring algorithm may have gotten stuck at a local maximum.
# Setting tau^2 = 0. Check the profile likelihood plot with profile().
# 10: Fisher scoring algorithm may have gotten stuck at a local maximum.
# Setting tau^2 = 0. Check the profile likelihood plot with profile().
# 11: Fisher scoring algorithm may have gotten stuck at a local maximum.
# Setting tau^2 = 0. Check the profile likelihood plot with profile().
# 12: Fisher scoring algorithm may have gotten stuck at a local maximum.
# Setting tau^2 = 0. Check the profile likelihood plot with profile().
# > 

str(metaOutput)
# #num [1:9221, 1:6] -0.06666 -0.00217 -0.04322 -0.22825 -0.07535 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:9221] "68980_498097" "21665_114521" "18591_24628" "57869_312258" ...
# ..$ : chr [1:6] "Log2FC_estimate" "SE" "pval" "CI_lb" ...

head(metaOutput)
#                 Log2FC_estimate       SE      pval      CI_lb      CI_ub            Number_Of_Comparisons
# 68980_498097     -0.066663135   0.14850066  0.6534980   -0.3577191  0.22439280                     3
# 21665_114521     -0.002172171   0.13661628  0.9873143   -0.2699352  0.26559082                     3
# 18591_24628      -0.043217996   0.09712981  0.6563552   -0.2335889  0.14715293                     3
# 57869_312258     -0.228248440   0.14073850  0.1048481   -0.5040908  0.04759395                     3
# 58240_313950     -0.075346790   0.07835606  0.3362531   -0.2289218  0.07822826                     3
# 104732_690422    -0.102069158   0.24596835  0.6781648   -0.5841583  0.38001996                     3

tail(metaOutput)
#                 Log2FC_estimate         SE       pval      CI_lb       CI_ub        Number_Of_Comparisons
# 13089_NA         0.104962213      0.16596382 0.52709920 -0.2203209 0.430245330                     3
# 100503946_NA     0.029836521      0.12541668 0.81195932 -0.2159757 0.275648696                     3
# 282663_NA        0.117705288      0.16831069 0.48434369 -0.2121776 0.447588177                     3
# 58242_NA        -0.066318309      0.03734082 0.07572848 -0.1395050 0.006868362                     3
# 171486_NA       -0.045576475      0.07793135 0.55866302 -0.1983191 0.107166158                     3
# 26446_NA        -0.006748082      0.13050693 0.95876240 -0.2625370 0.249040802                     3

########################################

## Multiple Comparison corrections
#The following code applies two different types of multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli) 
#Meta-analysis output with adjusted p-values is then outputted along with effect size information.

#We can add some additional gene annotation at this point too:
HOM_MouseVsRat <- read.csv("HOM_MouseVsRat.csv", header = TRUE, row.names = 1)
HOM_MouseVsRat$Mouse_EntrezGene.ID <- as.character(HOM_MouseVsRat$Mouse_EntrezGene.ID)
HOM_MouseVsRat$Rat_EntrezGene.ID <- as.character(HOM_MouseVsRat$Rat_EntrezGene.ID)

#9) Correct the meta-analysis output to take into account the fact that we are running the statistical calculations 
##          many times and therefore have a heightened risk of false discovery (false discovery rate correction) 

#Installing and loading relevant code packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
# Bioconductor version 3.17 (BiocManager 1.30.22), R 4.3.2 (2023-10-31)
# Bioconductor version '3.17' is out-of-date; the current release version '3.18' is available with R version '4.3'; see https://bioconductor.org/install

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")

library(BiocManager)
# Bioconductor version 3.18 (BiocManager 1.30.22), R 4.3.2 (2023-10-31)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("multtest") 

library(multtest)

#Let's functionalize it! =======================================================START OF FUNCTION ======================================================
FalseDiscoveryCorrection<-function(metaOutput, HOM_MouseVsRat, MetaAnalysis_Annotation){
  
  tempPvalAdjMeta<-mt.rawp2adjp(metaOutput[,3], proc=c("BH"))
  
  metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]
  
  metaOutputFDR<-cbind(metaOutput, metaPvalAdj[,2])
  
  colnames(metaOutputFDR)[7]<-"FDR"
  
  metaOutputFDR<<-metaOutputFDR
  
  print("metaOutputFDR:")
  print(str(metaOutputFDR))
  
  TempDF<-cbind.data.frame(metaOutputFDR, MetaAnalysis_Annotation)
  
  TempDF2<-join(TempDF, HOM_MouseVsRat[,c(4:5,9:11)], by="Mouse_EntrezGene.ID", type="left", match="first")
  
  TempDF3<-join(TempDF2, HOM_MouseVsRat[,c(15:16,20:22)], by="Rat_EntrezGene.ID", type="left", match="first")
  
  metaOutputFDR_annotated<-TempDF3
  metaOutputFDR_annotated<<-metaOutputFDR_annotated
  
  write.csv(metaOutputFDR_annotated, "metaOutputFDR_annotated.csv")
  
  #a version of the output in order by p-value:
  metaOutputFDR_OrderbyPval<<-metaOutputFDR_annotated[order(metaOutputFDR_annotated[,5]),]
  
  #Let's write out a version of the output in order by p-value:
  write.csv(metaOutputFDR_OrderbyPval, "metaOutputFDR_orderedByPval.csv")
  
  print("Do we have any genes that are statistically significant following false discovery rate correction?")
  print(sum(metaOutputFDR_annotated[,9]<0.10, na.rm=TRUE))
  
  print("What are the top results?")
  print(head(metaOutputFDR_annotated[order(metaOutputFDR_annotated[,5]),]))
  
  rm(tempPvalAdjMeta, metaPvalAdj)
  
}

## ====================================================================END OF FUNCTION ==========================================================================

#Example usage:
FalseDiscoveryCorrection(metaOutput, HOM_MouseVsRat, MetaAnalysis_Annotation)
# [1] "metaOutputFDR:"
# num [1:9221, 1:7] -0.06666 -0.00217 -0.04322 -0.22825 -0.07535 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:9221] "68980_498097" "21665_114521" "18591_24628" "57869_312258" ...
# ..$ : chr [1:7] "Log2FC_estimate" "SE" "pval" "CI_lb" ...
# NULL
# [1] "Do we have any genes that are statistically significant following false discovery rate correction?"
# [1] 60
# [1] "What are the top results?"
# Rat_EntrezGene.ID Mouse_EntrezGene.ID Log2FC_estimate         SE         pval       CI_lb      CI_ub Number_Of_Comparisons         FDR
# 66836_690285             690285               66836      -0.3129452 0.05839908 8.380645e-08 -0.42740529 -0.1984851                     3 0.000718361
# 20454_83505               83505               20454       0.5552179 0.10658646 1.897747e-07  0.34631225  0.7641235                     3 0.000718361
# 20466_363067             363067               20466      -0.5620947 0.10871600 2.337147e-07 -0.77517414 -0.3490152                     3 0.000718361
# 77531_314721             314721               77531       0.1329329 0.02690335 7.767260e-07  0.08020326  0.1856625                     3 0.001790548
# 260301_246044            246044              260301       0.5563502 0.11612277 1.659093e-06  0.32875375  0.7839466                     3 0.002576818
# 227154_501146            501146              227154       0.3665100 0.07653272 1.676706e-06  0.21650867  0.5165114                     3 0.002576818
# MouseVsRat_EntrezGene.ID Mouse_Symbol Mouse_Genetic.Location Mouse_Genome.Coordinates..mouse..GRCm39.human..GRCh38.p7.
# 66836_690285              66836_690285      Tmem223              Chr19  cM                                  Chr19:8748360-8749839(+)
# 20454_83505                20454_83505      St3gal5               Chr6  cM                                 Chr6:72074576-72131555(+)
# 20466_363067              20466_363067        Sin3a               Chr9  cM                                 Chr9:56979324-57035650(+)
# 77531_314721              77531_314721       Anks1b              Chr10  cM                                Chr10:89709371-90809162(+)
# 260301_246044            260301_246044         Otos               Chr1  cM                                 Chr1:92571940-92576563(-)
# 227154_501146            227154_501146       Stradb               Chr1  cM                                 Chr1:59012681-59034281(+)
# Mouse_Name Rat_Symbol Rat_Genetic.Location Rat_Genome.Coordinates..mouse..GRCm39.human..GRCh38.p7.
# 66836_690285                                    transmembrane protein 223    Tmem223             Chr1 q43                                                      NA
# 20454_83505            ST3 beta-galactoside alpha-2,3-sialyltransferase 5    St3gal5             Chr4 q31                                                      NA
# 20466_363067                     transcriptional regulator, SIN3A (yeast)      Sin3a             Chr8 q24                                                      NA
# 77531_314721  ankyrin repeat and sterile alpha motif domain containing 1B     Anks1b             Chr7 q13                                                      NA
# 260301_246044                                                 otospiralin       Otos             Chr9 q36                                                      NA
# 227154_501146                           STE20-related kinase adaptor beta     Stradb             Chr9 q31                                                      NA
# Rat_Name
# 66836_690285                                    transmembrane protein 223
# 20454_83505            ST3 beta-galactoside alpha-2,3-sialyltransferase 5
# 20466_363067                 SIN3 transcription regulator family member A
# 77531_314721  ankyrin repeat and sterile alpha motif domain containing 1B
# 260301_246044                                                 otospiralin
# 227154_501146                                  STE20 related adaptor beta

############################

#10) Determine which are the top differentially expressed genes and create forest plots to visualize the effect sizes for those top differentially expressed genes across the different studies. 

metaOutputFDR_OrderbyPval$Mouse_Symbol[c(1:17)]
# [1] "Tmem223" "St3gal5" "Sin3a"   "Anks1b"  "Otos"    "Stradb"  "Tecr"    NA        "Wipi2"   "Cd37"    "Mtss1"   "Wdr24"   "Gne"     "Npepl1"  "Psmd7"   "Vars1"  
# [17] "Guca2a"

metaOutputFDR_OrderbyPval$Mouse_EntrezGene.ID[c(1:17)]
#[1] "66836"  "20454"  "20466"  "77531"  "260301" "227154" "106529" "233544" "74781"  "12493"  "211401" "268933" "50798"  "228961" "17463"  "22321"  "14915" 

#Let's plot some of those top results!
#Quickly looking at the range of Log2FC values to figure out the limits for the x-axis for the forest plots:
hist(metaOutputFDR[,1], breaks=40)
#Range is mostly -.5 to .5

###################################### ANOTHER NEW FUNCTION =====================================================
MakeForestPlots<-function(EntrezIDAsCharacter){
  
  pdf(paste("ForestPlot_", EntrezIDAsCharacter, ".pdf", sep=""), height=5, width=8)
  
  effect<-as.numeric(MetaAnalysis_FoldChanges_ForMeta[MetaAnalysis_FoldChanges_ForMeta$Mouse_EntrezGene.ID==EntrezIDAsCharacter,-c(1:3)])
  var<-as.numeric(MetaAnalysis_SV_ForMeta[MetaAnalysis_FoldChanges_ForMeta$Mouse_EntrezGene.ID==EntrezIDAsCharacter,-c(1:3)])
  
  forest.rma(rma(effect, var),slab=colnames(MetaAnalysis_FoldChanges_ForMeta)[-c(1:3)],  xlim=c(-3, 3))
  
  mtext(paste(EntrezIDAsCharacter), line=-1.5, cex=2)
  dev.off()
}
##=================================END OF FUNCTION ===============================
library(ggplot2)
# =============================================================USE THE FUNCTION =============================================================
MakeForestPlots("66836")
MakeForestPlots("20454")
MakeForestPlots("20466")
MakeForestPlots("77531")
MakeForestPlots("260301")
MakeForestPlots("227154")
MakeForestPlots("106529")
MakeForestPlots("233544")

MakeForestPlots("74781")
MakeForestPlots("12493")
MakeForestPlots("211401")
MakeForestPlots("268933")
MakeForestPlots("50798")
MakeForestPlots("228961")
MakeForestPlots("17463")
MakeForestPlots("22321")
MakeForestPlots("14915")
