## 11.9.23 Following along with 2023_Formatting_limmaResults_forMetaAnalysis.R

#Start with an empty workspace
rm(list=ls())

#Set working directory
setwd("/Users/flandree/Documents/2023 2024 Sabbatical/01 UMich Project BrainAlchemy/3FileMeta")

#Load Libraries

#Clean annotation for results
TempResultsJoined<-read.csv("GSE68076Limma_Results.csv", header=TRUE,stringsAsFactors = FALSE)
str(TempResultsJoined)
# 'data.frame':	15575 obs. of  33 variables:
#   $ X                                                                                     : chr  "A_51_P100034" "A_51_P100174" "A_51_P100218" "A_51_P100238" ...
# $ Probe                                                                                   : chr  "A_51_P100034" "A_51_P100174" "A_51_P100218" "A_51_P100238" ...
# $ GeneSymbol                                                                              : chr  "Mif4gd" "Mns1" "Vmn1r234" "Or11l3" ...
# $ GeneName                                                                                : chr  "MIF4G domain containing" "meiosis-specific nuclear structural protein 1" "vomeronasal 1 receptor 234" "olfactory receptor family 11 subfamily L member 3" ...
# $ NCBIid                                                                                  : int  69674 17427 171232 258373 20807 21354 102633419 75163 232717 26406 ...
# $ X.1                                                                                     : chr  "A_51_P100034" "A_51_P100174" "A_51_P100218" "A_51_P100238" ...
# $ AveExpr                                                                                 : num  -1.542 0.424 4.371 2.676 -0.523 ...
# $ Coef..Intercept.                                                                        : num  -1.769 0.244 4.229 2.257 -0.642 ...
# $ Coef.GSE68076.treatment_factorsocial.stress..SS.                                        : num  0.5412 0.3279 0.5236 -0.0812 0.2375 ...
# $ Coef.GSE68076.StressDuration_Factor10                                                   : num  0.9538 0.6177 -0.1957 -0.4749 0.0764 ...
# $ Coef.GSE68076.DelayToSample_Numeric                                                     : num  -0.03085 -0.01295 -0.01445 0.01148 -0.00279 ...
# $ Coef.GSE68076.treatment_factorsocial.stress..SS..GSE68076.StressDuration_Factor10       : num  -1.468 -0.695 0.265 1.088 -0.143 ...
# $ Coef.GSE68076.treatment_factorsocial.stress..SS..GSE68076.DelayToSample_Numeric         : num  0.0418 0.00984 0.00935 0.04607 0.00493 ...
# $ t..Intercept.                                                                           : num  -7.92 1.52 14.98 3.68 -5.25 ...
# $ t.GSE68076.treatment_factorsocial.stress..SS.                                           : num  1.7194 1.4505 1.3158 -0.0941 1.3775 ...
# $ t.GSE68076.StressDuration_Factor10                                                      : num  2.327 2.098 -0.378 -0.422 0.34 ...
# $ t.GSE68076.DelayToSample_Numeric                                                        : num  -2.504 -1.463 -0.928 0.339 -0.414 ...
# $ t.GSE68076.treatment_factorsocial.stress..SS..GSE68076.StressDuration_Factor10          : num  -2.729 -1.799 0.391 0.737 -0.486 ...
# $ t.GSE68076.treatment_factorsocial.stress..SS..GSE68076.DelayToSample_Numeric            : num  2.541 0.833 0.45 1.02 0.548 ...
# $ P.value..Intercept.                                                                     : num  2.90e-09 1.38e-01 1.22e-16 7.86e-04 7.87e-06 ...
# $ P.value.GSE68076.treatment_factorsocial.stress..SS.                                     : num  0.0945 0.156 0.1969 0.9256 0.1772 ...
# $ P.value.GSE68076.StressDuration_Factor10                                                : num  0.026 0.0433 0.708 0.6755 0.7357 ...
# $ P.value.GSE68076.DelayToSample_Numeric                                                  : num  0.0172 0.1525 0.36 0.7363 0.6814 ...
# $ P.value.GSE68076.treatment_factorsocial.stress..SS..GSE68076.StressDuration_Factor10    : num  0.00993 0.08075 0.69857 0.46589 0.63033 ...
# $ P.value.GSE68076.treatment_factorsocial.stress..SS..GSE68076.DelayToSample_Numeric      : num  0.0157 0.4107 0.6558 0.3146 0.5875 ...
# $ P.value.adj..Intercept.                                                                 : num  7.94e-09 1.63e-01 9.30e-16 1.21e-03 1.51e-05 ...
# $ P.value.adj.GSE68076.treatment_factorsocial.stress..SS.                                 : num  0.515 0.602 0.644 0.984 0.62 ...
# $ P.value.adj.GSE68076.StressDuration_Factor10                                            : num  0.284 0.335 0.89 0.874 0.901 ...
# $ P.value.adj.GSE68076.DelayToSample_Numeric                                              : num  0.36 0.67 0.839 0.965 0.954 ...
# $ P.value.adj.GSE68076.treatment_factorsocial.stress..SS..GSE68076.StressDuration_Factor10: num  0.423 0.542 0.919 0.83 0.899 ...
# $ P.value.adj.GSE68076.treatment_factorsocial.stress..SS..GSE68076.DelayToSample_Numeric  : num  0.402 0.825 0.908 0.785 0.889 ...
# $ F                                                                                       : num  32.65 5.46 154.83 13.62 12.3 ...
# $ F.p.value                                                                               : num  8.09e-13 4.70e-04 2.10e-23 7.41e-08 2.35e-07

#Reading in the function.  =====================================================================================================================
#Print the number of rows in results, the # of rows with missing NCBI annotation, # of rows with "NA" 
##                        NCBI annotation, # of rows with missing gene symbol annotation, # of rows mapped to multiple NCBI_IDs, # of rows mapped to multiple Gene Symbols)

FilteringDEResults_GoodAnnotation <-function(TempResultsJoined){
  print(nrow(TempResultsJoined))
#15575
  print(sum(TempResultsJoined$NCBIid==""|TempResultsJoined$NCBIid=="null"))
  print(sum(is.na(TempResultsJoined$NCBIid)))
  #7
  print(sum(TempResultsJoined$GeneSymbol==""|TempResultsJoined$GeneSymbol=="null"))
  #7
  print(length(grep('\\|', TempResultsJoined$NCBIid)))
  #0
  print(length(grep('\\|', TempResultsJoined$GeneSymbol)))
  #0
  
  #include only subset of data with rows that do not contain an ncbi entrezid of ""
  TempResultsJoined_NoNA <- TempResultsJoined[(TempResultsJoined$NCBIid==""|TempResultsJoined$NCBIid=="nul")==FALSE & 
                                                is.na(TempResultsJoined$NCBIid)==FALSE,]
  if(length(grep('\\|', TempResultsJoined_NoNA$NCBIid))==0){
    TempResultsJoined_NoNA_NoMultimapped<<-TempResultsJoined_NoNA
  }else{
    TempResultsJoined_NoNA_NoMultimapped<<-TempResultsJoined_NoNA[-(grep('\\|', TempResultsJoined_NoNA$NCBIid)),]
  }
  
  ##Print # of rows with good annotation
  print(nrow(TempResultsJoined_NoNA_NoMultimapped))
  }
  
#Use the function with the first dataset ==========================================================================================================
FilteringDEResults_GoodAnnotation(TempResultsJoined)
# [1] 15575 rows left after filtering
# [1] NA 
# [1] 7 
# [1] 7
# [1] 0
# [1] 0
# [1] 15568

colnames(TempResultsJoined_NoNA_NoMultimapped)

# TempResultsJoined_NoNA_NoMultimapped_FoldChanges<-cbind(TempResultsJoined_NoNA_NoMultimapped$ ##Coef.SumExp_Subset_noBad_GSE92718_Filtered.treatment_factorcuff.operation##)
## update the part of the code between the two ## to be the column I want to actually use
TempResultsJoined_NoNA_NoMultimapped_FoldChanges<-cbind(TempResultsJoined_NoNA_NoMultimapped$Coef.GSE68076.treatment_factorsocial.stress..SS.)
####  [9] "Coef.GSE68076.treatment_factorsocial.stress..SS
###   use autofill to make sure you're not missing anything

str(TempResultsJoined_NoNA_NoMultimapped_FoldChanges)
#num [1:15568, 1] 0.5412 0.3279 0.5236 -0.0812 0.2375 ...

row.names(TempResultsJoined_NoNA_NoMultimapped_FoldChanges)<-TempResultsJoined_NoNA_NoMultimapped$NCBIid

##TempResultsJoined_NoNA_NoMultimapped_Tstats<-cbind(TempResultsJoined_NoNA_NoMultimapped$***t.SumExp_Subset_noBad_GSE92718_Filtered.treatment_factorcuff.operation***)
# swap name between the ***
TempResultsJoined_NoNA_NoMultimapped_Tstats<-cbind(TempResultsJoined_NoNA_NoMultimapped$t.GSE68076.treatment_factorsocial.stress..SS.)
str(TempResultsJoined_NoNA_NoMultimapped_Tstats)
#num [1:15568, 1] 1.7194 1.4505 1.3158 -0.0941 1.3775 ...

row.names(TempResultsJoined_NoNA_NoMultimapped_Tstats)<-TempResultsJoined_NoNA_NoMultimapped$NCBIid

## rename columns
ComparisonsOfInterest<-c("GSE68076_SSvsNS")
colnames(TempResultsJoined_NoNA_NoMultimapped_FoldChanges)<-ComparisonsOfInterest
colnames(TempResultsJoined_NoNA_NoMultimapped_Tstats)<-ComparisonsOfInterest

#Reading in the Extracting Results function:===================================================================================================================
ExtractingDEResults<-function(GSE_ID, TempResultsJoined_NoNA_NoMultimapped_FoldChanges, TempResultsJoined_NoNA_NoMultimapped_Tstats){
  #We calculate the standard error by dividing the log2FC by the tstat
  TempResultsJoined_NoNA_NoMultimapped_SE<-TempResultsJoined_NoNA_NoMultimapped_FoldChanges/TempResultsJoined_NoNA_NoMultimapped_Tstats
  str(TempResultsJoined_NoNA_NoMultimapped_SE)
  #For running our meta-analysis, we are actually going to need the sampling variance instead of the standard error
  #The sampling variance is just the standard error squared.
  TempResultsJoined_NoNA_NoMultimapped_SV<-(TempResultsJoined_NoNA_NoMultimapped_SE)^2
  str(TempResultsJoined_NoNA_NoMultimapped_SV)
  TempMasterResults<-list(Log2FC=TempResultsJoined_NoNA_NoMultimapped_FoldChanges, Tstat=TempResultsJoined_NoNA_NoMultimapped_Tstats, 
                          SE=TempResultsJoined_NoNA_NoMultimapped_SE, SV=TempResultsJoined_NoNA_NoMultimapped_SV)
  assign(paste("DEResults", GSE_ID, sep="_"), TempMasterResults, envir = as.environment(1))
  print(paste("Output: Named DEResults", GSE_ID, sep="_"))
  rm(TempMasterResults, TempResultsJoined_NoNA_NoMultimapped_SV, TempResultsJoined_NoNA_NoMultimapped_SE, 
     TempResultsJoined_NoNA_NoMultimapped_FoldChanges, TempResultsJoined_NoNA_NoMultimapped_Tstats)
} 
### END OF FUNCTION; read in the whole thing


#Use the extracting results function with the first dataset ============================DATA SET #1 ====================================================================
ExtractingDEResults(GSE_ID="GSE68076", TempResultsJoined_NoNA_NoMultimapped_FoldChanges, TempResultsJoined_NoNA_NoMultimapped_Tstats)
# num [1:15568, 1] 0.315 0.226 0.398 0.864 0.172 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:15568] "69674" "17427" "171232" "258373" ...
# ..$ : chr "GSE68076_SSvsNS"
# num [1:15568, 1] 0.0991 0.0511 0.1583 0.7462 0.0297 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:15568] "69674" "17427" "171232" "258373" ...
# ..$ : chr "GSE68076_SSvsNS"
# [1] "Output: Named DEResults_GSE68076"

str(DEResults_GSE68076)
#List of 4
# $ Log2FC: num [1:15568, 1] 0.5412 0.3279 0.5236 -0.0812 0.2375 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:15568] "69674" "17427" "171232" "258373" ...
# .. ..$ : chr "GSE68076_SSvsNS"
# $ Tstat : num [1:15568, 1] 1.7194 1.4505 1.3158 -0.0941 1.3775 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:15568] "69674" "17427" "171232" "258373" ...
# .. ..$ : chr "GSE68076_SSvsNS"
# $ SE    : num [1:15568, 1] 0.315 0.226 0.398 0.864 0.172 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:15568] "69674" "17427" "171232" "258373" ...
# .. ..$ : chr "GSE68076_SSvsNS"
# $ SV    : num [1:15568, 1] 0.0991 0.0511 0.1583 0.7462 0.0297 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:15568] "69674" "17427" "171232" "258373" ...
# .. ..$ : chr "GSE68076_SSvsNS"

#Sanity Check:
#Calculating the SE from the Log2FC and Tstat for the first gene(row)
# Log2FC/Tstat = SE pull numbers from the str output above.  THe first part is the Log2FC values, the second part is the Tstat, 
####  and then see if Log2FC / Tstat for a given example matches with the SE value listed in the chart.
0.5412/1.7194
#[1] 0.314761

#####   START NEXT DATASET HERE ===================================================================================== Dataset #2 ===============================================
#first clean up the environment
rm(TempResultsJoined, TempResultsJoined_NoNA_NoMultimapped, TempResultsJoined_NoNA_NoMultimapped_FoldChanges, TempResultsJoined_NoNA_NoMultimapped_Tstats, ComparisonsOfInterest)

#Running the next dataset: This part of the code will be dataset specific:
#Setting the working directory to where the limma results for the dataset are located:
setwd("~/Documents/2023 2024 Sabbatical/01 UMich Project BrainAlchemy/3FileMeta")

#Depending on file format, you may need to use read.delim instead of read.csv:
TempResultsJoined<-read.delim("Limma_results72262.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)

### Read in annotation file since 72262 didn't have one file with annotation plus actual data
TempResultsJoined_Annotation<-read.csv("Annotation_LimmaResultsGSE72262.csv", header=TRUE, stringsAsFactors = FALSE)

### combine the limma results with the annotation
library(plyr)
temp<-join(TempResultsJoined_Annotation, TempResultsJoined, by="X", type="right")
str(temp)
TempResultsJoined<-temp

str(TempResultsJoined)
# 'data.frame':	17718 obs. of  16 variables:
#   $ X                                                          : chr  "A_30_P01017530" "A_30_P01017570" "A_30_P01018255" "A_30_P01018819" ...
# $ AveExpr                                                    : num  -5.81 -5.78 -4.37 -6.07 -4.55 ...
# $ Coef..Intercept.                                           : num  -6.38 -5.76 -4.6 -5.4 -3.99 ...
# $ Coef.SummarizedExperiment_Subset.Stress_factorStress       : num  0.3855 0.0856 -0.0927 -0.5981 -0.6991 ...
# $ Coef.SummarizedExperiment_Subset.OvX_factorOvX             : num  0.749 -0.129 0.55 -0.746 -0.425 ...
# $ t..Intercept.                                              : num  -21.5 -18 -18.7 -17.7 -12.3 ...
# $ t.SummarizedExperiment_Subset.Stress_factorStress          : num  1.128 0.232 -0.326 -1.699 -1.861 ...
# $ t.SummarizedExperiment_Subset.OvX_factorOvX                : num  2.19 -0.35 1.93 -2.12 -1.13 ...
# $ P.value..Intercept.                                        : num  1.10e-17 7.23e-16 3.32e-16 1.12e-15 4.42e-12 ...
# $ P.value.SummarizedExperiment_Subset.Stress_factorStress    : num  0.27 0.8181 0.7473 0.1018 0.0745 ...
# $ P.value.SummarizedExperiment_Subset.OvX_factorOvX          : num  0.038 0.7289 0.0646 0.0442 0.2685 ...
# $ P.value.adj..Intercept.                                    : num  3.69e-17 1.76e-15 8.50e-16 2.63e-15 6.76e-12 ...
# $ P.value.adj.SummarizedExperiment_Subset.Stress_factorStress: num  0.427 0.887 0.838 0.213 0.17 ...
# $ P.value.adj.SummarizedExperiment_Subset.OvX_factorOvX      : num  0.0954 0.8106 0.142 0.1068 0.4041 ...
# $ F                                                          : num  387 328 315 399 197 ...
# $ F.p.value                                                  : num  4.20e-21 3.19e-20 5.15e-20 2.92e-21 1.52e-17 ...


#This part is not dataset specific (STILL ON DATASET #2):=======================================================================================
FilteringDEResults_GoodAnnotation(TempResultsJoined)
# [1] 17718  
# [1] NA       "# of rows with missing NCBI annotation:"
# [1] 1       "# of rows with missing Gene Symbol annotation:"
# [1] 1       "# of rows mapped to multiple NCBI_IDs:"
# [1] 0       
# [1] 0       "# of rows mapped to multiple Gene Symbols:"
# [1] 17717   "# of rows with good annotation"

colnames(TempResultsJoined_NoNA_NoMultimapped)
# [1] "X"                                                           "Probe"                                                      
# [3] "GeneSymbol"                                                  "GeneName"                                                   
# [5] "NCBIid"                                                      "AveExpr"                                                    
# [7] "Coef..Intercept."                                            "Coef.SummarizedExperiment_Subset.Stress_factorStress"       
# [9] "Coef.SummarizedExperiment_Subset.OvX_factorOvX"              "t..Intercept."                                              
# [11] "t.SummarizedExperiment_Subset.Stress_factorStress"           "t.SummarizedExperiment_Subset.OvX_factorOvX"                
# [13] "P.value..Intercept."                                         "P.value.SummarizedExperiment_Subset.Stress_factorStress"    
# [15] "P.value.SummarizedExperiment_Subset.OvX_factorOvX"           "P.value.adj..Intercept."                                    
# [17] "P.value.adj.SummarizedExperiment_Subset.Stress_factorStress" "P.value.adj.SummarizedExperiment_Subset.OvX_factorOvX"      
# [19] "F"                                                           "F.p.value"                                  

######===============================
#This part is (somewhat) dataset specific STILL ON DATASET #2
#We just need to replace the column names for the Log2FC ("Coef") and T-stats ("t.") With the appropriate names for this dataset:
#TempResultsJoined_NoNA_NoMultimapped_FoldChanges<-cbind(TempResultsJoined_NoNA_NoMultimapped$ FILL IN THE COEF NAME HERE)
TempResultsJoined_NoNA_NoMultimapped_FoldChanges<-cbind(TempResultsJoined_NoNA_NoMultimapped$Coef.SummarizedExperiment_Subset.Stress_factorStress)

str(TempResultsJoined_NoNA_NoMultimapped_FoldChanges)
#num [1:17717, 1] 0.3855 0.0856 -0.0927 -0.5981 -0.6991 ...

row.names(TempResultsJoined_NoNA_NoMultimapped_FoldChanges)<-TempResultsJoined_NoNA_NoMultimapped$NCBIid

#We also need to make a matrix of the Tstat values for the comparisons of interest
#These columns start with "t." and are specific to your comparisons of interest
#TempResultsJoined_NoNA_NoMultimapped_Tstats<-cbind(TempResultsJoined_NoNA_NoMultimapped$FIULL IN THE t. HERE)
TempResultsJoined_NoNA_NoMultimapped_Tstats<-cbind(TempResultsJoined_NoNA_NoMultimapped$t.SummarizedExperiment_Subset.Stress_factorStress)
str(TempResultsJoined_NoNA_NoMultimapped_Tstats)
# num [1:17717, 1] 1.128 0.232 -0.326 -1.699 -1.861 ...

row.names(TempResultsJoined_NoNA_NoMultimapped_Tstats)<-TempResultsJoined_NoNA_NoMultimapped$NCBIid

#Let's rename our columns to something nicer describing the effect of interest:
#Note - we later discovered that this name needs to include the dataset identifier (GSEID#) for later joining and plotting purposes
ComparisonsOfInterest<-c("GSE72262_StressVsNoStress")

colnames(TempResultsJoined_NoNA_NoMultimapped_FoldChanges)<-ComparisonsOfInterest
colnames(TempResultsJoined_NoNA_NoMultimapped_Tstats)<-ComparisonsOfInterest

#This is dataset specific, only because we need to provide the GSE#:
ExtractingDEResults(GSE_ID="GSE72262", TempResultsJoined_NoNA_NoMultimapped_FoldChanges, TempResultsJoined_NoNA_NoMultimapped_Tstats)
# # num [1:17717, 1] 0.342 0.368 0.284 0.352 0.376 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:17717] "751866" "69729" "102639188" "105246162" ...
# ..$ : chr "GSE72262_StressVsNoStress"
# num [1:17717, 1] 0.1168 0.1357 0.0809 0.124 0.1411 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:17717] "751866" "69729" "102639188" "105246162" ...
# ..$ : chr "GSE72262_StressVsNoStress"
# [1] "Output: Named DEResults_GSE72262"

str(DEResults_GSE72262)
# list of 4
# $ Log2FC: num [1:17717, 1] 0.3855 0.0856 -0.0927 -0.5981 -0.6991 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:17717] "751866" "69729" "102639188" "105246162" ...
# .. ..$ : chr "GSE72262_StressVsNoStress"
# $ Tstat : num [1:17717, 1] 1.128 0.232 -0.326 -1.699 -1.861 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:17717] "751866" "69729" "102639188" "105246162" ...
# .. ..$ : chr "GSE72262_StressVsNoStress"
# $ SE    : num [1:17717, 1] 0.342 0.368 0.284 0.352 0.376 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:17717] "751866" "69729" "102639188" "105246162" ...
# .. ..$ : chr "GSE72262_StressVsNoStress"
# $ SV    : num [1:17717, 1] 0.1168 0.1357 0.0809 0.124 0.1411 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:17717] "751866" "69729" "102639188" "105246162" ...
# .. ..$ : chr "GSE72262_StressVsNoStress"


####START THE NEXT DATASET HERE ========================================================================== DATASET #3==========================================================
#Clean up environment (save here as 11.16.23 format for Meta.R and 11.16.23 workspace)
rm(TempResultsJoined, TempResultsJoined_NoNA_NoMultimapped, TempResultsJoined_NoNA_NoMultimapped_FoldChanges, TempResultsJoined_NoNA_NoMultimapped_Tstats, ComparisonsOfInterest)

#Running the next dataset: This part of the code will be dataset specific:
#Setting the working directory to where the limma results for the dataset are located:
setwd("~/Documents/2023 2024 Sabbatical/01 UMich Project BrainAlchemy/3FileMeta")

#Depending on file format, you may need to use read.delim instead of read.csv:
TempResultsJoined<-read.delim("Limma_results84185.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)

### Read in annotation file since 72262 didn't have one file with annotation plus actual data
TempResultsJoined_Annotation<-read.csv("Annotation_LimmaResultsGSE84185.csv", header=TRUE, stringsAsFactors = FALSE)

### combine the limma results with the annotation
library(plyr)
temp<-join(TempResultsJoined_Annotation, TempResultsJoined, by="X", type="right")
str(temp)

## overwrite temp to be TempResultsJoined since that's the object name in the remaining code.
TempResultsJoined<-temp

str(TempResultsJoined)
# 'data.frame':	14414 obs. of  24 variables:
#   $ X                                                                                                                                                     : chr  "10003" "10011" "10015" "10016" ...
# $ Probe                                                                                                                                                 : chr  "10003" "10011" "10015" "10016" ...
# $ GeneSymbol                                                                                                                                            : chr  "Sirt3" "Folr1" "Mkln1" "Ppp1r36dn" ...
# $ GeneName                                                                                                                                              : chr  "sirtuin 3" "folate receptor alpha" "muskelin 1, intracellular mediator containing kelch motifs" "Ppp1r36 downstream neighbor" ...
# $ NCBIid                                                                                                                                                : int  64384 14275 27418 100041694 80987 13607 20845 102162 231130 74044 ...
# $ AveExpr                                                                                                                                               : num  7.85 4.69 4.68 4.67 4.67 ...
# $ Coef..Intercept.                                                                                                                                      : num  7.89 4.66 4.66 4.66 4.67 ...
# $ Coef.SummarizedExperiment_Filtered..1...Stress_factorunpredictable.chronic.mild.stress                                                                : num  0.052275 -0.005788 -0.006475 -0.000463 -0.00265 ...
# $ Coef.SummarizedExperiment_Filtered..1...Drug_factorfluoxetine                                                                                         : num  -0.0161 0.0253 0.054 0.0352 -0.0159 ...
# $ Coef.SummarizedExperiment_Filtered..1...Stress_factorunpredictable.chronic.mild.stress.SummarizedExperiment_Filtered..1...Drug_factorfluoxetine       : num  -0.23819 0.075 -0.03187 -0.03665 0.00805 ...
# $ t..Intercept.                                                                                                                                         : num  106 101 231 278 456 ...
# $ t.SummarizedExperiment_Filtered..1...Stress_factorunpredictable.chronic.mild.stress                                                                   : num  0.4978 -0.0883 -0.2273 -0.0195 -0.183 ...
# $ t.SummarizedExperiment_Filtered..1...Drug_factorfluoxetine                                                                                            : num  -0.153 0.387 1.894 1.489 -1.097 ...
# $ t.SummarizedExperiment_Filtered..1...Stress_factorunpredictable.chronic.mild.stress.SummarizedExperiment_Filtered..1...Drug_factorfluoxetine          : num  -1.604 0.81 -0.791 -1.095 0.393 ...
# $ P.value..Intercept.                                                                                                                                   : num  1.72e-40 8.97e-40 1.07e-50 4.01e-53 1.29e-59 ...
# $ P.value.SummarizedExperiment_Filtered..1...Stress_factorunpredictable.chronic.mild.stress                                                             : num  0.622 0.93 0.822 0.985 0.856 ...
# $ P.value.SummarizedExperiment_Filtered..1...Drug_factorfluoxetine                                                                                      : num  0.8793 0.7016 0.0678 0.147 0.2813 ...
# $ P.value.SummarizedExperiment_Filtered..1...Stress_factorunpredictable.chronic.mild.stress.SummarizedExperiment_Filtered..1...Drug_factorfluoxetine    : num  0.119 0.425 0.435 0.282 0.697 ...
# $ P.value.adj..Intercept.                                                                                                                               : num  3.45e-40 1.70e-39 4.85e-50 2.29e-52 1.79e-58 ...
# $ P.value.adj.SummarizedExperiment_Filtered..1...Stress_factorunpredictable.chronic.mild.stress                                                         : num  0.886 0.996 0.969 1 0.976 ...
# $ P.value.adj.SummarizedExperiment_Filtered..1...Drug_factorfluoxetine                                                                                  : num  0.967 0.855 0.227 0.344 0.469 ...
# $ P.value.adj.SummarizedExperiment_Filtered..1...Stress_factorunpredictable.chronic.mild.stress.SummarizedExperiment_Filtered..1...Drug_factorfluoxetine: num  0.784 0.787 0.789 0.784 0.896 ...
# $ F                                                                                                                                                     : num  11182 10253 53931 77768 207579 ...
# $ F.p.value                                                                                                                                             : num  1.82e-47 6.75e-47 8.50e-58 3.35e-60 1.19e-66 ...
# > 60 1.19e-66 ...


#This part is not dataset specific:=======================================================================================
FilteringDEResults_GoodAnnotation(TempResultsJoined)
# [1] 14414
# [1] NA
# [1] 1
# [1] 1
# [1] 0
# [1] 0
# [1] 14413

colnames(TempResultsJoined_NoNA_NoMultimapped)
# [1] "X"                                                                                                                                                     
# [2] "Probe"                                                                                                                                                 
# [3] "GeneSymbol"                                                                                                                                            
# [4] "GeneName"                                                                                                                                              
# [5] "NCBIid"                                                                                                                                                
# [6] "AveExpr"                                                                                                                                               
# [7] "Coef..Intercept."                                                                                                                                      
# [8] "Coef.SummarizedExperiment_Filtered..1...Stress_factorunpredictable.chronic.mild.stress"                                                                
# [9] "Coef.SummarizedExperiment_Filtered..1...Drug_factorfluoxetine"                                                                                         
# [10] "Coef.SummarizedExperiment_Filtered..1...Stress_factorunpredictable.chronic.mild.stress.SummarizedExperiment_Filtered..1...Drug_factorfluoxetine"       
# [11] "t..Intercept."                                                                                                                                         
# [12] "t.SummarizedExperiment_Filtered..1...Stress_factorunpredictable.chronic.mild.stress"                                                                   
# [13] "t.SummarizedExperiment_Filtered..1...Drug_factorfluoxetine"                                                                                            
# [14] "t.SummarizedExperiment_Filtered..1...Stress_factorunpredictable.chronic.mild.stress.SummarizedExperiment_Filtered..1...Drug_factorfluoxetine"          
# [15] "P.value..Intercept."                                                                                                                                   
# [16] "P.value.SummarizedExperiment_Filtered..1...Stress_factorunpredictable.chronic.mild.stress"                                                             
# [17] "P.value.SummarizedExperiment_Filtered..1...Drug_factorfluoxetine"                                                                                      
# [18] "P.value.SummarizedExperiment_Filtered..1...Stress_factorunpredictable.chronic.mild.stress.SummarizedExperiment_Filtered..1...Drug_factorfluoxetine"    
# [19] "P.value.adj..Intercept."                                                                                                                               
# [20] "P.value.adj.SummarizedExperiment_Filtered..1...Stress_factorunpredictable.chronic.mild.stress"                                                         
# [21] "P.value.adj.SummarizedExperiment_Filtered..1...Drug_factorfluoxetine"                                                                                  
# [22] "P.value.adj.SummarizedExperiment_Filtered..1...Stress_factorunpredictable.chronic.mild.stress.SummarizedExperiment_Filtered..1...Drug_factorfluoxetine"
# [23] "F"                                                                                                                                                     
# [24] "F.p.value"                                 

######===============================
#This part is (somewhat) dataset specific STILL ON DATASET #3
#We just need to replace the column names for the Log2FC ("Coef") and T-stats ("t.") With the appropriate names for this dataset:
#TempResultsJoined_NoNA_NoMultimapped_FoldChanges<-cbind(TempResultsJoined_NoNA_NoMultimapped$ FILL IN THE COEF NAME HERE)
TempResultsJoined_NoNA_NoMultimapped_FoldChanges<-cbind(TempResultsJoined_NoNA_NoMultimapped$Coef.SummarizedExperiment_Filtered..1...Stress_factorunpredictable.chronic.mild.stress)

str(TempResultsJoined_NoNA_NoMultimapped_FoldChanges)
##### num [1:14413, 1] 0.052275 -0.005788 -0.006475 -0.000463 -0.00265 ...

row.names(TempResultsJoined_NoNA_NoMultimapped_FoldChanges)<-TempResultsJoined_NoNA_NoMultimapped$NCBIid

#We also need to make a matrix of the Tstat values for the comparisons of interest
#These columns start with "t." and are specific to your comparisons of interest
#TempResultsJoined_NoNA_NoMultimapped_Tstats<-cbind(TempResultsJoined_NoNA_NoMultimapped$FIULL IN THE t. HERE)
TempResultsJoined_NoNA_NoMultimapped_Tstats<-cbind(TempResultsJoined_NoNA_NoMultimapped$t.SummarizedExperiment_Filtered..1...Stress_factorunpredictable.chronic.mild.stress)
str(TempResultsJoined_NoNA_NoMultimapped_Tstats)
####  num [1:14413, 1] 0.4978 -0.0883 -0.2273 -0.0195 -0.183 ...

row.names(TempResultsJoined_NoNA_NoMultimapped_Tstats)<-TempResultsJoined_NoNA_NoMultimapped$NCBIid

#Let's rename our columns to something nicer describing the effect of interest:
#Note - we later discovered that this name needs to include the dataset identifier (GSEID#) for later joining and plotting purposes
ComparisonsOfInterest<-c("GSE84185_StressVsNoStress")

colnames(TempResultsJoined_NoNA_NoMultimapped_FoldChanges)<-ComparisonsOfInterest
colnames(TempResultsJoined_NoNA_NoMultimapped_Tstats)<-ComparisonsOfInterest

#This is dataset specific, only because we need to provide the GSE# STILL ON DATASET #3:
ExtractingDEResults(GSE_ID="GSE84185", TempResultsJoined_NoNA_NoMultimapped_FoldChanges, TempResultsJoined_NoNA_NoMultimapped_Tstats)
# num [1:14413, 1] 0.105 0.0655 0.0285 0.0237 0.0145 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:14413] "64384" "14275" "27418" "100041694" ...
# ..$ : chr "GSE84185_StressVsNoStress"
# num [1:14413, 1] 0.011029 0.004291 0.000811 0.00056 0.00021 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:14413] "64384" "14275" "27418" "100041694" ...
# ..$ : chr "GSE84185_StressVsNoStress"
# [1] "Output: Named DEResults_GSE84185"

str(DEResults_GSE84185)
#List of 4
# $ Log2FC: num [1:14413, 1] 0.052275 -0.005788 -0.006475 -0.000463 -0.00265 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:14413] "64384" "14275" "27418" "100041694" ...
# .. ..$ : chr "GSE84185_StressVsNoStress"
# $ Tstat : num [1:14413, 1] 0.4978 -0.0883 -0.2273 -0.0195 -0.183 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:14413] "64384" "14275" "27418" "100041694" ...
# .. ..$ : chr "GSE84185_StressVsNoStress"
# $ SE    : num [1:14413, 1] 0.105 0.0655 0.0285 0.0237 0.0145 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:14413] "64384" "14275" "27418" "100041694" ...
# .. ..$ : chr "GSE84185_StressVsNoStress"
# $ SV    : num [1:14413, 1] 0.011029 0.004291 0.000811 0.00056 0.00021 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : chr [1:14413] "64384" "14275" "27418" "100041694" ...
# .. ..$ : chr "GSE84185_StressVsNoStress"

rm(TempResultsJoined, TempResultsJoined_NoNA_NoMultimapped, TempResultsJoined_NoNA_NoMultimapped_FoldChanges, TempResultsJoined_NoNA_NoMultimapped_Tstats, ComparisonsOfInterest)