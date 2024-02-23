#Megan Hagenauer Sept 21, 2023
#Flandreau Datset GSE72262
#This code document is set up to provide an example of analyzing a single Gemma dataset
################### Overall Design of this study
# Female mice were subjected to chronic ultra-mild stress, a bilateral ovariectomy, or both. 
#Sham-operated mice without stress were used as the control. 
#Medial prefrontal cortex and whole blood were obtained from the same individual (n = 6 in each group), 
#and analyzed using an Agilent SurePrint G3 Mouse GE 8×60K Microarray (Design ID: 028005)

########################### highlight text then use command shift c on mac to add # to each line

##############################
#Setting the working directory

setwd ("/Users/flandree/Documents/2023 2024 Sabbatical/01 UMich Project BrainAlchemy/GSE72262 Analysis")

#######################Gemma applies a standardized filter across all platforms (RNA-Seq/microarray) to remove genes with minimal variance
#These genes are likely to suffer from severe floor effects (very low expression/unmeasurable)
#Or not actually be expressed at all and the measurements that we have are just noise.
#"There are two stages to our variance filter. 
# First is just to remove genes with zero variance. 
# The second is to remove genes that have too many identical values (within a small tolerance), 
# where "too many" is currently defined as <70% distinct values, but we have tweaked this over time (slightly - we're never perfectly happy with anything)."
#I call this a variance filter as a shorthand, of course it's not that exactly.
#But the latter filter was originally intended to address cases where microarray data had been clipped by the submitter.
#Excluding low level expressed genes within an individual dataset reduces false positives due to noise within low-variance data
#It also decreases the severity of the multiple comparisons correction (because there are fewer genes in the dataset)
#Within a meta-analysis, we do not need to worry about effects caused by noise in individual datasets quite as much

###############load useful code packages (must do this for each new workspace; every time you create a new environment of objects)
library(SummarizedExperiment)
library(gemma.R)
library(plyr)
library(tidyr)

#Input the summarized experiment object for your dataset: GSE72262 (FYI This takes a really long time and then shows up in the "data" area)
SummarizedExperiment_Filtered<-gemma.R::get_dataset_object("GSE72262", type = 'se', filter=TRUE, consolidate="average")
SummarizedExperiment_Filtered
#### Copied from console (comment / uncomment lines under the code tab)
# $`9654`
# class: SummarizedExperiment 
# dim: 17718 48 
# metadata(8): title abstract ... GemmaSuitabilityScore taxon
# assays(1): counts
# rownames(17718): A_30_P01017530 A_30_P01017570 ... Averaged from A_55_P2084248 A_66_P140742 Averaged from A_55_P2056463
# A_66_P140794
# rowData names(4): Probe GeneSymbol GeneName NCBIid
# colnames(48): WB_OVX_5 WB_control_6 ... mPFC_OVX_1 mPFC_control_2
# colData names(3): factorValues treatment organism part

ExpressionData_Filtered<-assay(SummarizedExperiment_Filtered[[1]])
str(ExpressionData_Filtered)
# num [1:17718, 1:48] -5.67 -6.54 -3.75 -4.04 -5.45 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:17718] "A_30_P01017530" "A_30_P01017570" "A_30_P01018255" "A_30_P01018819" ...
# ..$ : chr [1:48] "WB_OVX_5" "WB_control_6" "mPFC_OVX+stress_1" "WB_stress_2" ...

hist(ExpressionData_Filtered)
####Each row = gene or probe; different expression values for each probe / row.  Gemma does log 2 regression on datasets it
# can access.  Some it can't access (agilent) then it uses "external data" values that aren't standardized because we don't know what the "external analyst"
# did with the data before it got to Gemma.  That can result in massive differences in expression data OR a duplication in the log transformation, which
# creates an artificially low level of variability, which is also bad.  The histogram helps tell us if we need to update the transformation to have the 'correct'
# degree of variability.). CURRENT histogram -5 to +10 looks like log 2 transformed values (typical = -5 to +12; if i see a range from -1 to +2.5 it's probably
# double log transformed.  If the range is 100s or 1000s it likely hasn't been log 2 transformed at all)

min(ExpressionData_Filtered)
# [1] -6.754

#Let's look at the mean vs. variance curve (corresponds to line 114 - 116 on GSE68076 file).
#These two lines create vectors for the mean per gene and the standard deviation per gene:
ExpressionData_Filtered_MeanPerGene<-apply(ExpressionData_Filtered, 1, mean)
ExpressionData_Filtered_SDPerGene<-apply(ExpressionData_Filtered, 1, sd)

#This line generates a graph of the mean vs. the standard deviation.  If there is equal variance, it will basically be a horizontal line.
plot(ExpressionData_Filtered_SDPerGene~ExpressionData_Filtered_MeanPerGene)
#heteroskedasticity may be a problem in this dataset
#The trend/voom function will help this data still be useable in regression equations <- this will come up later

#How many genes have zero variance?
sum(ExpressionData_Filtered_SDPerGene==0)
### the answer is zero, so that's good.

###########################Dataset Subsetting:
#Before we do much more with the dataset, let's subset down to the samples that we actually plan to use:
#First, we need to know what we have:
#How to access different parts of the Summarized Experiment object:

colData(SummarizedExperiment_Filtered[[1]])
#DataFrame with 48 rows and 3 columns
# factorValues
# <list>
#   WB_OVX_5                               bilateral ovariectomy:NA:NA:...,reference subject role:http://purl.obolibra..:NA:...,blood:http://purl.obolibra..:NA:...,...
# WB_control_6      reference subject role:http://purl.obolibra..:NA:...,reference subject role:http://purl.obolibra..:NA:...,blood:http://purl.obolibra..:NA:...,...
# mPFC_OVX+stress_1                                      bilateral ovariectomy:NA:NA:...,chronic:http://purl.obolibra..:NA:...,mild:http://purl.obolibra..:NA:...,...
# WB_stress_2                       chronic:http://purl.obolibra..:NA:...,mild:http://purl.obolibra..:NA:...,stimulus or stress d..:http://www.ebi.ac.uk..:NA:...,...
# mPFC_OVX+stress_3                                      bilateral ovariectomy:NA:NA:...,chronic:http://purl.obolibra..:NA:...,mild:http://purl.obolibra..:NA:...,...
# ...                                                                                                                                                             ...
# mPFC_OVX+stress_4                                      bilateral ovariectomy:NA:NA:...,chronic:http://purl.obolibra..:NA:...,mild:http://purl.obolibra..:NA:...,...
# mPFC_stress_5                     chronic:http://purl.obolibra..:NA:...,mild:http://purl.obolibra..:NA:...,stimulus or stress d..:http://www.ebi.ac.uk..:NA:...,...
# mPFC_OVX+stress_6                                      bilateral ovariectomy:NA:NA:...,chronic:http://purl.obolibra..:NA:...,mild:http://purl.obolibra..:NA:...,...
# mPFC_OVX_1                              bilateral ovariectomy:NA:NA:...,medial to:http://purl.obolibra..:NA:...,prefrontal cortex:http://purl.obolibra..:NA:...,...
# mPFC_control_2     reference subject role:http://purl.obolibra..:NA:...,medial to:http://purl.obolibra..:NA:...,prefrontal cortex:http://purl.obolibra..:NA:...,...
# treatment          organism part
# <character>            <character>
#   WB_OVX_5          bilateral ovariectom..                  blood
# WB_control_6      reference subject ro..                  blood
# mPFC_OVX+stress_1 bilateral ovariectom.. medial to,prefrontal..
# WB_stress_2       chronic,mild,stimulu..                  blood
# mPFC_OVX+stress_3 bilateral ovariectom.. medial to,prefrontal..
# ...                                  ...                    ...
# mPFC_OVX+stress_4 bilateral ovariectom.. medial to,prefrontal..
# mPFC_stress_5     chronic,mild,stimulu.. medial to,prefrontal..
# mPFC_OVX+stress_6 bilateral ovariectom.. medial to,prefrontal..
# mPFC_OVX_1        bilateral ovariectom.. medial to,prefrontal..
# mPFC_control_2    reference subject ro.. medial to,prefrontal..

head(colData(SummarizedExperiment_Filtered[[1]]))
# DataFrame with 6 rows and 3 columns
# factorValues
# <list>
#   WB_OVX_5                               bilateral ovariectomy:NA:NA:...,reference subject role:http://purl.obolibra..:NA:...,blood:http://purl.obolibra..:NA:...,...
# WB_control_6      reference subject role:http://purl.obolibra..:NA:...,reference subject role:http://purl.obolibra..:NA:...,blood:http://purl.obolibra..:NA:...,...
# mPFC_OVX+stress_1                                      bilateral ovariectomy:NA:NA:...,chronic:http://purl.obolibra..:NA:...,mild:http://purl.obolibra..:NA:...,...
# WB_stress_2                       chronic:http://purl.obolibra..:NA:...,mild:http://purl.obolibra..:NA:...,stimulus or stress d..:http://www.ebi.ac.uk..:NA:...,...
# mPFC_OVX+stress_3                                      bilateral ovariectomy:NA:NA:...,chronic:http://purl.obolibra..:NA:...,mild:http://purl.obolibra..:NA:...,...
# mPFC_stress_4                     chronic:http://purl.obolibra..:NA:...,mild:http://purl.obolibra..:NA:...,stimulus or stress d..:http://www.ebi.ac.uk..:NA:...,...
# treatment          organism part
# <character>            <character>
#   WB_OVX_5          bilateral ovariectom..                  blood
# WB_control_6      reference subject ro..                  blood
# mPFC_OVX+stress_1 bilateral ovariectom.. medial to,prefrontal..
# WB_stress_2       chronic,mild,stimulu..                  blood
# mPFC_OVX+stress_3 bilateral ovariectom.. medial to,prefrontal..
# mPFC_stress_4     chronic,mild,stimulu.. medial to,prefrontal..

rowData(SummarizedExperiment_Filtered[[1]])
# DataFrame with 17718 rows and 4 columns
# Probe    GeneSymbol               GeneName      NCBIid
# <character>   <character>            <character> <character>
#   A_30_P01017530                                   A_30_P01017530          Kis2 Kaplan integration s..      751866
# A_30_P01017570                                   A_30_P01017570 2410003L11Rik RIKEN cDNA 2410003L1..       69729
# A_30_P01018255                                   A_30_P01018255       Gm27008  predicted gene, 27008   102639188
# A_30_P01018819                                   A_30_P01018819       Gm41496  predicted gene, 41496   105246162
# A_30_P01018893                                   A_30_P01018893         Ndc80 NDC80 kinetochore co..       67052
# ...                                                         ...           ...                    ...         ...
# Averaged from A_55_P2106255 A_66_P140335 Averaged from A_55_P..         Mbnl3 muscleblind like spl..      171170
# Averaged from A_55_P2157966 A_66_P140452 Averaged from A_55_P..         Map1a microtubule-associat..       17754
# Averaged from A_51_P286946 A_66_P140507  Averaged from A_51_P..          Lhpp phospholysine phosph..       76429
# Averaged from A_55_P2084248 A_66_P140742 Averaged from A_55_P..         Azin2   antizyme inhibitor 2      242669
# Averaged from A_55_P2056463 A_66_P140794 Averaged from A_55_P..       Gm10584   predicted gene 10584   100043682

head(rowData(SummarizedExperiment_Filtered[[1]]))
######### Columns = Probe, GeneSymbol, GeneName, and NCBIid; Rows = Individual samples? <-- NOT SURE HERE
# DataFrame with 6 rows and 4 columns
# Probe    GeneSymbol               GeneName      NCBIid
# <character>   <character>            <character> <character>
#   A_30_P01017530 A_30_P01017530          Kis2 Kaplan integration s..      751866
# A_30_P01017570 A_30_P01017570 2410003L11Rik RIKEN cDNA 2410003L1..       69729
# A_30_P01018255 A_30_P01018255       Gm27008  predicted gene, 27008   102639188
# A_30_P01018819 A_30_P01018819       Gm41496  predicted gene, 41496   105246162
# A_30_P01018893 A_30_P01018893         Ndc80 NDC80 kinetochore co..       67052
# A_30_P01018927 A_30_P01018927       Gm36480  predicted gene, 36480   102640415


##################The distribution of samples from different brain regions or blood
table(SummarizedExperiment_Filtered[[1]]$`organism part`)
# blood medial to,prefrontal cortex 
# 24                          24 

################ #The distribution of the different phenotypes in this experiment
table(SummarizedExperiment_Filtered[[1]]$treatment)
# bilateral ovariectomy,chronic,mild,stimulus or stress design                  bilateral ovariectomy,reference subject role 
# 12                                                            12 
# chronic,mild,stimulus or stress design,reference subject role                 reference subject role,reference subject role 
# 12                                                            12

SampleFilter<-
  SummarizedExperiment_Filtered[[1]]$`organism part`=="blood"

#Subsetting the data to have only blood:
SummarizedExperiment_Subset<-SummarizedExperiment_Filtered[[1]][,SampleFilter]

SummarizedExperiment_Subset
# class: SummarizedExperiment 
# dim: 17718 24 
# metadata(8): title abstract ... GemmaSuitabilityScore taxon
# assays(1): counts
# rownames(17718): A_30_P01017530 A_30_P01017570 ... Averaged from A_55_P2084248 A_66_P140742 Averaged from A_55_P2056463 A_66_P140794
# rowData names(4): Probe GeneSymbol GeneName NCBIid
# colnames(24): WB_OVX_5 WB_control_6 ... WB_OVX+stress_2 WB_stress_3
# colData names(3): factorValues treatment organism part

table(SummarizedExperiment_Subset$`organism part`)
#blood
#24

table(SummarizedExperiment_Subset$treatment)
# bilateral ovariectomy,chronic,mild,stimulus or stress design                  bilateral ovariectomy,reference subject role 
# 6                                                             6 
# chronic,mild,stimulus or stress design,reference subject role                 reference subject role,reference subject role 
# 6                                                             6 

#### So the first one is OVX, Stress; the second on the same line is OVX, Non_Stress; The third is Sham, Stress; the fourth is Sham, No Stress
#I will probably need to separate out the OVX vs. Sham variable from the Stress vs. No-Stress variable similar to separating hte 
#stress variable and time variable in the dataset for GSE68076

#Pulling out a matrix of gene expression data for this subset (to use in functions that require matrices).  Correspond to line 199 of GSE68076
ExpressionData_Subset<-assay(SummarizedExperiment_Subset)
str(ExpressionData_Subset)
# num [1:17718, 1:24] -5.67 -6.54 -3.75 -4.04 -5.45 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:17718] "A_30_P01017530" "A_30_P01017570" "A_30_P01018255" "A_30_P01018819" ...
# ..$ : chr [1:24] "WB_OVX_5" "WB_control_6" "WB_stress_2" "WB_control_4" ...

#Outlier Removal:
#This would be a good time to check for outliers and remove them if they are present.
#Creating an example sample-sample correlation scatterplot (data for all genes for 1 sample versus the data for all genes for the second sample)
plot(ExpressionData_Subset[,1]~ExpressionData_Subset[,2])

########################10.24.23################################################
#Creating a matrix showing how each sample correlates with every other sample:
CorMatrix<-cor(ExpressionData_Subset)
heatmap(CorMatrix)
boxplot(CorMatrix)
ExpressionData_Subset<-assay(SummarizedExperiment_Subset)
str(ExpressionData_Subset)
# num [1:17718, 1:24] -5.67 -6.54 -3.75 -4.04 -5.45 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:17718] "A_30_P01017530" "A_30_P01017570" "A_30_P01018255" "A_30_P01018819" ...
# ..$ : chr [1:24] "WB_OVX_5" "WB_control_6" "WB_stress_2" "WB_control_4" ...

ExpressionData_Subset_SDperGene<-apply(ExpressionData_Subset, 1, sd)
min(ExpressionData_Subset_SDperGene)
##[1] 0.02966256
sum(ExpressionData_Subset_SDperGene==0)
#[1] 0

tapply(ExpressionData_Subset[1,], SummarizedExperiment_Subset$treatment, sd)
# bilateral ovariectomy,chronic,mild,stimulus or stress design                  bilateral ovariectomy,reference subject role 
#                      0.7467044                                                     1.1757341 
# chronic,mild,stimulus or stress design,reference subject role                 reference subject role,reference subject role 
#                       0.6464559                                                     0.8318001 

#We want the minimum sd for all of our groups to not be zero:
min(tapply(ExpressionData_Subset[1,], SummarizedExperiment_Subset$treatment, sd))!=0
# TRUE

#... and we need to know that for all rows (this next line of code takes a long time to run)
GenesToFilter<-apply(ExpressionData_Subset, 1, function(y) (min(tapply(y, SummarizedExperiment_Subset$treatment, sd))!=0))
head(GenesToFilter)
# A_30_P01017530 A_30_P01017570 A_30_P01018255 A_30_P01018819 A_30_P01018893 A_30_P01018927 
#      TRUE           TRUE           TRUE           TRUE           TRUE           TRUE

#How many genes we'll end up keeping
sum(GenesToFilter)
#[1] 17718 ##there were no genes with zero variability so we have the same number as before running this.
SummarizedExperiment_Subset
# class: SummarizedExperiment 
# dim: 17718 24 
# metadata(8): title abstract ... GemmaSuitabilityScore taxon
# assays(1): counts
# rownames(17718): A_30_P01017530 A_30_P01017570 ... Averaged from A_55_P2084248 A_66_P140742 Averaged from A_55_P2056463 A_66_P140794
# rowData names(4): Probe GeneSymbol GeneName NCBIid
# colnames(24): WB_OVX_5 WB_control_6 ... WB_OVX+stress_2 WB_stress_3
# colData names(3): factorValues treatment organism part

str(ExpressionData_Subset)
# num [1:17718, 1:24] -5.67 -6.54 -3.75 -4.04 -5.45 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:17718] "A_30_P01017530" "A_30_P01017570" "A_30_P01018255" "A_30_P01018819" ...
# ..$ : chr [1:24] "WB_OVX_5" "WB_control_6" "WB_stress_2" "WB_control_4" ...

#Checking for batch confound:
#Since we are down to the subset of samples that we plan to use, this would be a good time to check for confounds
#Unfortunately, processing batches are often unbalanced in regards to variables of interest
#For this dataset, Gemma has all of the processing batch information lumped into one variable

table(SummarizedExperiment_Subset$treatment)
# bilateral ovariectomy,chronic,mild,stimulus or stress design                  bilateral ovariectomy,reference subject role 
#                     6                                                             6 
# chronic,mild,stimulus or stress design,reference subject role                 reference subject role,reference subject role 
#                     6                                                             6

#### =======================================Separating out the stress variable from the OvX variable
SummarizedExperiment_Subset$Stress<- SummarizedExperiment_Subset$treatment
SummarizedExperiment_Subset$Stress[SummarizedExperiment_Subset$treatment=="bilateral ovariectomy,chronic,mild,stimulus or stress design"]<-"Stress"
SummarizedExperiment_Subset$Stress[SummarizedExperiment_Subset$treatment=="chronic,mild,stimulus or stress design,reference subject role"]<-"Stress"
SummarizedExperiment_Subset$Stress[SummarizedExperiment_Subset$treatment=="bilateral ovariectomy,reference subject role"]<-"NS"
SummarizedExperiment_Subset$Stress[SummarizedExperiment_Subset$treatment=="reference subject role,reference subject role"]<-"NS"
table(SummarizedExperiment_Subset$Stress)
# NS Stress 
# 12     12 

SummarizedExperiment_Subset$OvX<- SummarizedExperiment_Subset$treatment
SummarizedExperiment_Subset$OvX[SummarizedExperiment_Subset$treatment=="bilateral ovariectomy,chronic,mild,stimulus or stress design"]<-"OvX"
SummarizedExperiment_Subset$OvX[SummarizedExperiment_Subset$treatment=="bilateral ovariectomy,reference subject role"]<-"OvX"
SummarizedExperiment_Subset$OvX[SummarizedExperiment_Subset$treatment=="chronic,mild,stimulus or stress design,reference subject role"]<-"OvControl"
SummarizedExperiment_Subset$OvX[SummarizedExperiment_Subset$treatment=="reference subject role,reference subject role"]<-"OvControl"
table(SummarizedExperiment_Subset$OvX)
# OvControl       OvX 
# 12        12

###================== Ok now let's go back and do principal component analysis
#Principal components analysis:
#Let's see what the largest sources of variation are in the dataset
#We are particularly interested in determining which variables we should include as co-variates in our model
#E.g., if a technical variable (e.g., batch) is the main source of variation in the data, we should probably control for it

library(stats)
pca_output<-prcomp(t(ExpressionData_Subset), scale=TRUE)

PCeigenvectors<-pca_output$rotation[ ,c(1:4)]
PCeigenvectors2<-cbind(PCeigenvectors, rowData(SummarizedExperiment_Subset))
write.csv(PCeigenvectors2, "PCeigenvectors.csv")

PC1<-pca_output$x[,1]
PC2<-pca_output$x[,2]
PC3<-pca_output$x[,3]
PC4<-pca_output$x[,4]

#Output a scree plot for the PCA:
#This plot illustrates the proportion of variance explained by each principal component (PC):
png("10 PCA Scree Plot1.png")
plot(summary(pca_output)$importance[2,]~(c(1:length(summary(pca_output)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Proportion of Variance Explained", col=2)
dev.off()
##### These plots show up in the folder I assigned.  The def.off() makes sure the file is not overwritten by the next figure.

#Output a scatterplot illustrating the relationship between Principal components 1 & 2 (PC1 & PC2):
png("10 PC1 vs PC2.png")
plot(PC1~PC2, main="Principal Components Analysis")
dev.off()

########================= SAVING AS HERE ================================================
##You can color these plots in using different variables in the dataset 
#This can help explain the main sources of variation in the data
#Technical variables (brain region, batch) tend to be the main culprits
#plot(PC1~PC2, main="Principal Components Analysis", col=as.factor(BatchVariables$Run))
#PC1 and 2, don't seem to be related to our main batch variable <-- this is empty for me because I don't have "run" as a variable

plot(PC1~PC2, main="Principal Components Analysis", col=as.factor(SummarizedExperiment_Subset$Stress))
#Here we see PC1 and PC2 for treatment (stress vs. non stress).  
plot(PC1~PC2, main="Principal Components Analysis", col=as.factor(SummarizedExperiment_Subset$OvX))
#Here we see PC1 and PC2 for treatment (OvX).  


plot(PC3~PC4, main="Principal Components Analysis", col=as.factor(SummarizedExperiment_Subset$Stress))
plot(PC3~PC4, main="Principal Components Analysis", col=as.factor(SummarizedExperiment_Subset$OvX))

#If we want to zoom in on the relationship between PC1 and treatment we can make a boxplot
boxplot(PC1~SummarizedExperiment_Subset$Stress, las=2, xlab="")
boxplot(PC2~SummarizedExperiment_Subset$Stress, las=2, xlab="")
boxplot(PC3~SummarizedExperiment_Subset$Stress, las=2, xlab="")
boxplot(PC4~SummarizedExperiment_Subset$Stress, las=2, xlab="")

#### Let's do the same for OvX
boxplot(PC1~SummarizedExperiment_Subset$OvX, las=2, xlab="")
boxplot(PC2~SummarizedExperiment_Subset$OvX, las=2, xlab="")
boxplot(PC3~SummarizedExperiment_Subset$OvX, las=2, xlab="")
boxplot(PC4~SummarizedExperiment_Subset$OvX, las=2, xlab="")
###PC1 definitely looks like it could be OvX related

str(SummarizedExperiment_Subset)

boxplot(PC1~SummarizedExperiment_Subset$Stress*SummarizedExperiment_Subset$OvX)

#Getting the expression data for a particular gene (Xist specific to females; Ddx3y specific to males):
Xist<-assay(SummarizedExperiment_Subset)[rowData(SummarizedExperiment_Subset)$GeneSymbol=="Xist",]
Ddx3y<-assay(SummarizedExperiment_Subset)[rowData(SummarizedExperiment_Subset)$GeneSymbol=="Ddx3y",]
plot(Ddx3y~Xist)
###### Error in (function (formula, data = NULL, subset = NULL, na.action = na.fail,  : 
######variable lengths differ (found for 'Xist')

#We can add jittered data points to our graph so that we can see the values for the individual samples in each group.
#To do this, we use the function stripchart and the exact same y~x formula that we used for the boxplot
#The parameter pch determines the size of the data points.
#The parameter "add" places the points on top of the boxplot that we already created.
boxplot(PC1~SummarizedExperiment_Subset$Stress*SummarizedExperiment_Subset$OvX, xlab="Stress and OvX", main="GSE72262", col=2)
stripchart(PC1~SummarizedExperiment_Subset$Stress*SummarizedExperiment_Subset$OvX, pch = 19, method = "jitter", jitter = 0.2, vertical = TRUE, add=TRUE)

###################### <-- EIF resume here
#Now: Running a differential expression analysis for an entire dataset:
#Making a design matrix:

table(SummarizedExperiment_Subset$Stress)
# NS Stress 
# 12     12
table(SummarizedExperiment_Subset$OvX)
# OvControl       OvX 
# 12        12 

str(SummarizedExperiment_Subset$Stress)
#chr [1:24] "NS" "NS" "Stress" "NS" "NS" "Stress" "NS" "Stress" "Stress" "Stress" "Stress" "NS" "Stress" "NS" "NS" "NS" "NS" ...
### this tells us it's a character vector

###============== 
library(limma)
# Attaching package: ‘limma’
# The following object is masked from ‘package:BiocGenerics’:
#   plotMA

str(colData(SummarizedExperiment_Subset))

#Converting our (categorical) variables of interest to a factor
SummarizedExperiment_Subset$Stress_factor<-as.factor(SummarizedExperiment_Subset$Stress)
levels(SummarizedExperiment_Subset$Stress_factor)
#[1] "NS"     "Stress"
SummarizedExperiment_Subset$OvX_factor<-as.factor(SummarizedExperiment_Subset$OvX)
levels(SummarizedExperiment_Subset$OvX_factor)
#[1] "OvControl" "OvX" 

str(colData(SummarizedExperiment_Subset))
# Formal class 'DFrame' [package "S4Vectors"] with 6 slots
# ..@ rownames       : chr [1:24] "WB_OVX_5" "WB_control_6" "WB_stress_2" "WB_control_4" ...
# ..@ nrows          : int 24
# ..@ elementType    : chr "ANY"
# ..@ elementMetadata: NULL
# ..@ metadata       : list()
# ..@ listData       :List of 7
# .. ..$ factorValues :List of 24

###I don't see a variable for block or batch.
table(SummarizedExperiment_Subset$Stress_factor,SummarizedExperiment_Subset$OvX_factor)
#           OvControl OvX
# NS             6   6
# Stress         6   6


#### SAVE AS 1025

#You can add co-variates to the design matrix using an addition symbol, e.g. ~treatment_factor+block
### This model matrix includes the main factors of stress and OvX but also the interaction
### Can set things so there's an interaction term but stress is still the default is "contrast treatment" so the reference level is important if there's an interaction
##### maybe not have an interaction term for this one... avoid underpowered problems
##======================================================= 10.30.23 UPDATE LATER 
design <- model.matrix(~SummarizedExperiment_Subset$Stress_factor+SummarizedExperiment_Subset$OvX_factor+
                         SummarizedExperiment_Subset$Stress_factor*SummarizedExperiment_Subset$OvX_factor
                      , data=colData(SummarizedExperiment_Subset))
design
# (Intercept) SummarizedExperiment_Subset$Stress_factorStress SummarizedExperiment_Subset$OvX_factorOvX
# WB_OVX_5                  1                                               0                                         1
# WB_control_6              1                                               0                                         0
# WB_stress_2               1                                               1                                         0
# WB_control_4              1                                               0                                         0
# WB_OVX_3                  1                                               0                                         1
# WB_stress_6               1                                               1                                         0
# WB_control_1              1                                               0                                         0
# WB_OVX+stress_5           1                                               1                                         1
# WB_OVX+stress_1           1                                               1                                         1
# WB_OVX+stress_3           1                                               1                                         1
# WB_stress_4               1                                               1                                         0
# WB_OVX_6                  1                                               0                                         1
# WB_stress_1               1                                               1                                         0
# WB_OVX_2                  1                                               0                                         1
# WB_control_3              1                                               0                                         0
# WB_OVX_4                  1                                               0                                         1
# WB_control_5              1                                               0                                         0
# WB_OVX+stress_4           1                                               1                                         1
# WB_stress_5               1                                               1                                         0
# WB_OVX+stress_6           1                                               1                                         1
# WB_OVX_1                  1                                               0                                         1
# WB_control_2              1                                               0                                         0
# WB_OVX+stress_2           1                                               1                                         1
# WB_stress_3               1                                               1                                         0
# SummarizedExperiment_Subset$Stress_factorStress:SummarizedExperiment_Subset$OvX_factorOvX
# WB_OVX_5                                                                                                0
# WB_control_6                                                                                            0
# WB_stress_2                                                                                             0
# WB_control_4                                                                                            0
# WB_OVX_3                                                                                                0
# WB_stress_6                                                                                             0
# WB_control_1                                                                                            0
# WB_OVX+stress_5                                                                                         1
# WB_OVX+stress_1                                                                                         1
# WB_OVX+stress_3                                                                                         1
# WB_stress_4                                                                                             0
# WB_OVX_6                                                                                                0
# WB_stress_1                                                                                             0
# WB_OVX_2                                                                                                0
# WB_control_3                                                                                            0
# WB_OVX_4                                                                                                0
# WB_control_5                                                                                            0
# WB_OVX+stress_4                                                                                         1
# WB_stress_5                                                                                             0
# WB_OVX+stress_6                                                                                         1
# WB_OVX_1                                                                                                0
# WB_control_2                                                                                            0
# WB_OVX+stress_2                                                                                         1
# WB_stress_3                                                                                             0
# attr(,"assign")
# [1] 0 1 2 3
# attr(,"contrasts")
# attr(,"contrasts")$`SummarizedExperiment_Subset$Stress_factor`
# [1] "contr.treatment"
# 
# attr(,"contrasts")$`SummarizedExperiment_Subset$OvX_factor`
# [1] "contr.treatment"

######################SAMPLE CODE FROM MH
# #You can add co-variates to the design matrix using an addition symbol, e.g. ~treatment_factor+ timepoint (If I want to use all time points)
#           design <- model.matrix(~treatment_factor + timepoint, data=colData(SummarizedExperiment_Subset))
#           design
# #testing it out with the data for 1 gene to make sure it doesn't crash: 
# #we have to tell the lm function to not add an intercept
# summary.lm(lm(ExpressionData_Subset_noBad_Filtered[1,]~0+design))
#Applying the model to all genes using limma:

fit <- lmFit(ExpressionData_Subset, design)

#Adding an eBayes correction to help reduce the influence of outliers/small sample size on estimates
efit <- eBayes(fit, trend=TRUE)
dt<-decideTests(efit)
summary(decideTests(efit))
#              (Intercept) SummarizedExperiment_Subset$Stress_factorStress SummarizedExperiment_Subset$OvX_factorOvX
# Down         11821                                            1853                                       944
# NotSig         985                                           15037                                     15917
# Up            4912                                             828                                       857
# SummarizedExperiment_Subset$Stress_factorStress:SummarizedExperiment_Subset$OvX_factorOvX
# Down                                                                                           0
# NotSig                                                                                     17718
# Up                                                                                             0

str(efit)

write.fit(efit, adjust="BH", file="Limma_results.txt")
write.csv(rowData(SummarizedExperiment_Subset), "Annotation_LimmaResultsGSE72262.csv")


#########===========================10.30.23 GSE72262_1030.R
> setwd ("/Users/flandree/Documents/2023 2024 Sabbatical/01 UMich Project BrainAlchemy/GSE72262 Analysis")
> load("~/Documents/2023 2024 Sabbatical/01 UMich Project BrainAlchemy/GSE72262 Analysis/GSE72262 Workspace.RData")

design <- model.matrix(~SummarizedExperiment_Subset$Stress_factor+SummarizedExperiment_Subset$OvX_factor
                       , data=colData(SummarizedExperiment_Subset))
design
# (Intercept) SummarizedExperiment_Subset$Stress_factorStress SummarizedExperiment_Subset$OvX_factorOvX
# WB_OVX_5                  1                                               0                                         1
# WB_control_6              1                                               0                                         0
# WB_stress_2               1                                               1                                         0
# WB_control_4              1                                               0                                         0
# WB_OVX_3                  1                                               0                                         1
# WB_stress_6               1                                               1                                         0
# WB_control_1              1                                               0                                         0
# WB_OVX+stress_5           1                                               1                                         1
# WB_OVX+stress_1           1                                               1                                         1
# WB_OVX+stress_3           1                                               1                                         1
# WB_stress_4               1                                               1                                         0
# WB_OVX_6                  1                                               0                                         1
# WB_stress_1               1                                               1                                         0
# WB_OVX_2                  1                                               0                                         1
# WB_control_3              1                                               0                                         0
# WB_OVX_4                  1                                               0                                         1
# WB_control_5              1                                               0                                         0
# WB_OVX+stress_4           1                                               1                                         1
# WB_stress_5               1                                               1                                         0
# WB_OVX+stress_6           1                                               1                                         1
# WB_OVX_1                  1                                               0                                         1
# WB_control_2              1                                               0                                         0
# WB_OVX+stress_2           1                                               1                                         1
# WB_stress_3               1                                               1                                         0
# attr(,"assign")
# [1] 0 1 2
# attr(,"contrasts")
# attr(,"contrasts")$`SummarizedExperiment_Subset$Stress_factor`
# [1] "contr.treatment"
# 
# attr(,"contrasts")$`SummarizedExperiment_Subset$OvX_factor`
# [1] "contr.treatment"



fit <- lmFit(ExpressionData_Subset, design)

#Adding an eBayes correction to help reduce the influence of outliers/small sample size on estimates
efit <- eBayes(fit, trend=TRUE)
dt<-decideTests(efit)
summary(decideTests(efit))
# (Intercept) SummarizedExperiment_Subset$Stress_factorStress SummarizedExperiment_Subset$OvX_factorOvX
# Down         11918                                            3602                                      2097
# NotSig         850                                           12750                                     12218
# Up            4950                                            1366                                      3403
str(efit)

write.fit(efit, adjust="BH", file="Limma_results72262.txt")
write.csv(rowData(SummarizedExperiment_Subset), "Annotation_LimmaResultsGSE72262.csv")
