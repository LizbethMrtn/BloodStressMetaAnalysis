#Megan Hagenauer July 13, 2023
#Flandreau Datset GSE84185 10.25.23

##################################
#This code document is set up to provide an example of analyzing a single Gemma dataset
#To keep things neat and tidy, I'm only analyzing one dataset in this document & workspace


#Setting the working directory
setwd("/Users/flandree/Documents/2023 2024 Sabbatical/01 UMich Project BrainAlchemy/GSE84185 Analysis")

#load useful code packages (only need to download once, but need to load in each new workspace)
library(SummarizedExperiment)
library(gemma.R)
library(plyr)
library(tidyr)

#Input the summarized experiment object for your dataset: GSE84185 (FYI This takes a really long time and then shows up in the "data" area)
SummarizedExperiment_Filtered<-gemma.R::get_dataset_object("GSE84185", type = 'se', filter=TRUE, consolidate="average")
SummarizedExperiment_Filtered
# $`18398`
# class: SummarizedExperiment 
# dim: 14414 32 
# metadata(8): title abstract ... GemmaSuitabilityScore taxon
# assays(1): counts
# rownames(14414): 10003 10011 ... Averaged from 6027 998 Averaged from 34179 999
# rowData names(4): Probe GeneSymbol GeneName NCBIid
# colnames(32): Blood-S-T_32 Blood-NS-FLX_46 ... Blood-NS-T_52 Blood-S-T_25
# colData names(3): factorValues environmental stress treatment

ExpressionData_Filtered<-assay(SummarizedExperiment_Filtered[[1]])
str(ExpressionData_Filtered)
# num [1:14414, 1:32] 7.69 4.67 4.65 4.67 4.66 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:14414] "10003" "10011" "10015" "10016" ...
# ..$ : chr [1:32] "Blood-S-T_32" "Blood-NS-FLX_46" "Blood-NS-T_53" "Blood-NS-T_64" ...

hist(ExpressionData_Filtered)
#Produces a histogram of expression data (MH 8.31.23 Frequency = number of values of expression data matrix.  
# Each row = gene or probe; different expression values for each probe / row.  Gemma does log 2 regression on datasets it can access.  
# Some it can't access (agilent) then it uses "external data" values that aren't standardized because we don't know what the "external analyst"
# did with the data before it got to Gemma.  
# That can result in massive differences in expression data OR a duplication in the log transformation, which
# creates an artificially low level of variability, which is also bad.  
# The histogram helps tell us if we need to update the transformation to have the 'correct' degree of variability.). 
# EIF histogram -5 to +7 looks like log 2 transformed values (typical = -5 to +12; if i see a range from -1 to +2.5 it's probably double log transformed.  
# If the range is 100s or 1000s it likely hasn't been log 2 transformed at all)
##### Range is ~+5 to +15, I think that's ok?


#MH 8.23 another way to judge range of values to determine if it's correctly log transformed.  
#The smallest value in an RNA seq dataset is probably artificial so that they don't have to worry about zeros)
#This is probably an artificial floor
#A small positive value (e.g., 0.5) is often added to each datapoint before log transformation Because you can't log transform zeroes
min(ExpressionData_Filtered)
#[1] 4.6443

#Let's look at the mean vs. variance curve:
ExpressionData_Filtered_MeanPerGene<-apply(ExpressionData_Filtered, 1, mean)
ExpressionData_Filtered_SDPerGene<-apply(ExpressionData_Filtered, 1, sd)
plot(ExpressionData_Filtered_SDPerGene~ExpressionData_Filtered_MeanPerGene)
#heteroskedasticity = unequal variance.  Ideally this will be horizontal line.  The scattered dots could be a problem.
#The trend/voom function will help this data still be useable in regression equations <- this will come up later

#How many genes have zero variance?
sum(ExpressionData_Filtered_SDPerGene==0)
# [1] 0 <-- this is good
sum(c(TRUE, FALSE, TRUE, FALSE)) ##<-- not sure why we do this or if it's necessary
#[1] 2

#Dataset Subsetting: Bfore we do much more with the dataset, let's subset down to the samples that we actually plan to use
#First, we need to know what we have:
#How to access different parts of the Summarized Experiment object:

colData(SummarizedExperiment_Filtered[[1]])
#DataFrame with 32 rows and 3 columns

head(colData(SummarizedExperiment_Filtered[[1]]))
# DataFrame with 6 rows and 3 columns
## Columns = factorValues (as list); environmental stress (as character), treatment (as character)

rowData(SummarizedExperiment_Filtered[[1]])
# DataFrame with 14414 rows and 4 columns
# Probe    GeneSymbol               GeneName      NCBIid

head(rowData(SummarizedExperiment_Filtered[[1]]))
# DataFrame with 6 rows and 4 columns
# Probe  GeneSymbol               GeneName      NCBIid
# <character> <character>            <character> <character>

#The distribution of samples from strss or non stress
table(SummarizedExperiment_Filtered[[1]]$`environmental stress`)
# reference subject role unpredictable chronic mild stress 
#        16                                16

#The distribution of the different drug treatments
table(SummarizedExperiment_Filtered[[1]]$treatment)
# fluoxetine       reference substance role <-- does this mean the reference is backwards?<<<<<<<<<<<<<<<<<<<<<<<<
#      16                       16 

##Checking to see what is the reference treatment
levels(SummarizedExperiment_Filtered[[1]]$treatment) ### Null because this is in as character not factor, will need to change later.

#### I don't think I need to remove anything because it's all blood.STOPPED HERE 10.25.23 due to error with next code.  
### MH meeting 10am 10.26.23

#Outlier Removal:
#This would be a good time to check for outliers and remove them if they are present.
#Creating an example sample-sample correlation scatterplot (data for all genes for 1 sample versus the data for all genes for the second sample)
plot(ExpressionData_Filtered[,1]~ExpressionData_Filtered[,2])
CorMatrix<-cor(ExpressionData_Filtered)
CorMatrix

#Writing that matrix out to your working directory to save it:
write.csv(CorMatrix, "CorMatrix.csv")
#### open this file (from working directory in excel to find the specific ones that are outliers
#### conditional formatting in excel...)

#Creating a hierarchically clustered heatmap illustrating the sample-sample correlation matrix:
heatmap(CorMatrix)
##### white may be outliers.

#Creating a boxplot illustrating the sample-sample correlations for each sample. Outliers should be obvious at this point.
boxplot(CorMatrix)
#I don't see any outliers. 
#######SAVED 

#This function calculates the StDev for each treatment group for a particular gene (row of data) because still can't div by zero and we 
##don't want to crash:
tapply(ExpressionData_Filtered[1,], SummarizedExperiment_Filtered[[1]]$treatment, sd)
# fluoxetine reference substance role 
# 0.2240557                0.2072331 
tapply(ExpressionData_Filtered[1,], SummarizedExperiment_Filtered[[1]]$`environmental stress`, sd)
# reference subject role unpredictable chronic mild stress 
# 0.1847168                         0.2576648 

#... and we need to know that for all rows (this next line of code takes a long time to run)
GenesToFilter<-apply(ExpressionData_Filtered, 1, function(y) (min(tapply(y, SummarizedExperiment_Filtered[[1]]$treatment, sd))!=0))
head(GenesToFilter)
# 10003 10011 10015 10016 10018 10024 
# TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
sum(GenesToFilter)
#####[1] 14414 <-- I think this means there are no genes with sd=0

GenesToFilter2<-apply(ExpressionData_Filtered, 1, function(y) (min(tapply(y, SummarizedExperiment_Filtered[[1]]$`environmental stress`, sd))!=0))
head(GenesToFilter2)
# 10003 10011 10015 10016 10018 10024 
# TRUE  TRUE  TRUE  TRUE  TRUE  TRUE 
sum(GenesToFilter2)
# [1] 14414 <-- I think this means there are no genes with sd=0

ExpressionData_Filtered<-assay(SummarizedExperiment_Filtered[[1]])
str(ExpressionData_Filtered)
# num [1:14414, 1:32] 7.69 4.67 4.65 4.67 4.66 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:14414] "10003" "10011" "10015" "10016" ...
# ..$ : chr [1:32] "Blood-S-T_32" "Blood-NS-FLX_46" "Blood-NS-T_53" "Blood-NS-T_64" ...

###### The next part of the code is looking for batch but we don't have 'batch' as a variable here.  Moving on to PCA
#Principal components analysis:
#Let's see what the largest sources of variation are in the dataset
#We are particularly interested in determining which variables we should include as co-variates in our model
#E.g., if a technical variable (e.g., batch) is the main source of variation in the data, we should probably control for it

library(stats)
pca_output<-prcomp(t(ExpressionData_Filtered), scale=TRUE)

PCeigenvectors<-pca_output$rotation[ ,c(1:4)]
PCeigenvectors2<-cbind(PCeigenvectors, rowData(SummarizedExperiment_Filtered[[1]]))
write.csv(PCeigenvectors2, "PCeigenvectors.csv")

PC1<-pca_output$x[,1]
PC2<-pca_output$x[,2]
PC3<-pca_output$x[,3]
PC4<-pca_output$x[,4]

#Output a scree plot for the PCA (no outliers):
#This plot illustrates the proportion of variance explained by each principal component (PC):
png("10 PCA Scree Plot1.png")
plot(summary(pca_output)$importance[2,]~(c(1:length(summary(pca_output)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Proportion of Variance Explained", col=2)
dev.off()

#Output a scatterplot illustrating the relationship between Principal components 1 & 2 (PC1 & PC2):
png("10 PC1 vs PC2.png")
plot(PC1~PC2, main="Principal Components Analysis")
dev.off()

plot(PC1~PC2, main="Principal Components Analysis", col=as.factor(SummarizedExperiment_Filtered[[1]]$`environmental stress`))
#Here we see PC1 and PC2 for treatment (stress vs. non stress).  These don't seem to be responsible for PC1 and 2
plot(PC1~PC2, main="Principal Components Analysis", col=as.factor(SummarizedExperiment_Filtered[[1]]$treatment))
#Here we see PC1 and PC2 for treatment (fluoxetine). HO LY SHIT clearly fluoxetine has a massive impact. 

plot(PC3~PC4, main="Principal Components Analysis", col=as.factor(SummarizedExperiment_Filtered[[1]]$`environmental stress`))
### WTAF, this looks like there's two massive outliers
plot(PC3~PC4, main="Principal Components Analysis", col=as.factor(SummarizedExperiment_Filtered[[1]]$treatment))
### This is crazy.  Fluoxetine's impact is massive.

#If we want to zoom in on the relationship between PC1 and treatment we can make a boxplot
boxplot(PC1~SummarizedExperiment_Filtered[[1]]$`environmental stress`, las=2, xlab="")
boxplot(PC2~SummarizedExperiment_Filtered[[1]]$`environmental stress`, las=2, xlab="")
boxplot(PC3~SummarizedExperiment_Filtered[[1]]$`environmental stress`, las=2, xlab="")
boxplot(PC4~SummarizedExperiment_Filtered[[1]]$`environmental stress`, las=2, xlab="")
#### how do I find out who this outlier is???

boxplot(PC1~SummarizedExperiment_Filtered[[1]]$treatment, las=2, xlab="")
boxplot(PC2~SummarizedExperiment_Filtered[[1]]$treatment, las=2, xlab="")
boxplot(PC3~SummarizedExperiment_Filtered[[1]]$treatment, las=2, xlab="")
boxplot(PC4~SummarizedExperiment_Filtered[[1]]$treatment, las=2, xlab="")
### OK so PC1 and PC2 are clearly carried by drug

str(SummarizedExperiment_Filtered[[1]])

boxplot(PC1~SummarizedExperiment_Filtered[[1]]$treatment*SummarizedExperiment_Filtered[[1]]$`environmental stress`)

#We can add jittered data points to our graph so that we can see the values for the individual samples in each group.
#To do this, we use the function stripchart and the exact same y~x formula that we used for the boxplot
#The parameter pch determines the size of the data points.
#The parameter "add" places the points on top of the boxplot that we already created.
boxplot(PC1~SummarizedExperiment_Filtered[[1]]$treatment*SummarizedExperiment_Filtered[[1]]$`environmental stress`, xlab="Stress and Drug", main="GSE84185", col=2)

boxplot(PC1~SummarizedExperiment_Filtered[[1]]$`environmental stress`*SummarizedExperiment_Filtered[[1]]$treatment, xlab="Stress and Drug", main="GSE84185", col=2)
stripchart(PC1~SummarizedExperiment_Filtered[[1]]$`environmental stress`*SummarizedExperiment_Filtered[[1]]$treatment, pch = 19, method = "jitter", jitter = 0.2, vertical = TRUE, add=TRUE)

Xist<-assay(SummarizedExperiment_Filtered[[1]])[rowData(SummarizedExperiment_Filtered[[1]])$GeneSymbol=="Xist",]
Ddx3y<-assay(SummarizedExperiment_Filtered[[1]])[rowData(SummarizedExperiment_Filtered[[1]])$GeneSymbol=="Ddx3y",]
plot(Ddx3y~Xist). ###<-- I got an error again here as well

###============== 
library(limma)
# Attaching package: ‘limma’
# The following object is masked from ‘package:BiocGenerics’:
#   plotMA

str(colData(SummarizedExperiment_Filtered[[1]]))

#Converting our (categorical) variables of interest to a factor
SummarizedExperiment_Filtered[[1]]$Stress_factor<-as.factor(SummarizedExperiment_Filtered[[1]]$`environmental stress`)
levels(SummarizedExperiment_Filtered[[1]]$Stress_factor)
#[1] "reference subject role"            "unpredictable chronic mild stress"

SummarizedExperiment_Filtered[[1]]$Drug_factor<-as.factor(SummarizedExperiment_Filtered[[1]]$treatment)
levels(SummarizedExperiment_Filtered[[1]]$Drug_factor)
#[1] "fluoxetine"               "reference substance role"

str(colData(SummarizedExperiment_Filtered[[1]]))

table(SummarizedExperiment_Filtered[[1]]$Stress_factor,SummarizedExperiment_Filtered[[1]]$Drug_factor)
#                                     fluoxetine        reference substance role
# reference subject role                     8                        8
# unpredictable chronic mild stress          8                        8

####################SAVE HERE GSE84185 Script1026c and workspace1026c
### Ok now trying to update the reference level for the Drug_factor 

SummarizedExperiment_Filtered[[1]]$Drug_factor <-relevel(SummarizedExperiment_Filtered[[1]]$Drug_factor, ref="reference substance role")
table(SummarizedExperiment_Filtered[[1]]$Stress_factor,SummarizedExperiment_Filtered[[1]]$Drug_factor)
#### OMG I DID IT
#                                      reference substance role    fluoxetine
# reference subject role                                   8          8
# unpredictable chronic mild stress                        8          8


#Create model matrix with two main factors: stress and drug
design <- model.matrix(~SummarizedExperiment_Filtered[[1]]$Stress_factor+SummarizedExperiment_Filtered[[1]]$Drug_factor 
                       , data=colData(SummarizedExperiment_Filtered[[1]]))


##Create model matrix with two main factors plus interaction (Just replaced the "+" with "*" and R will expand to the full model)
###Contrast Sum = more traditional ANOVA output (coefficients are not interpretable, in hypothetical space between groups...)
###Contrast Treatment = preference in social science (?) coefficients are interpretable in the model.  e.g. it would be the reference sample fluoxetine in co-efficient
### MH prefers contrast treatment and avoiding interaction terms to streamline interpretation of the model in small n studies 
###############(balance with strong biological rationale for including interaction)
design <- model.matrix(~SummarizedExperiment_Filtered[[1]]$Stress_factor*SummarizedExperiment_Filtered[[1]]$Drug_factor 
                       , data=colData(SummarizedExperiment_Filtered[[1]]))

design
# (Intercept) SummarizedExperiment_Filtered[[1]]$Stress_factorunpredictable chronic mild stress SummarizedExperiment_Filtered[[1]]$Drug_factorfluoxetine
# Blood-S-T_32              1                                                                                 1                                                        0
# Blood-NS-FLX_46           1                                                                                 0                                                        1
# Blood-NS-T_53             1                                                                                 0                                                        0
# Blood-NS-T_64             1                                                                                 0                                                        0...

str(SummarizedExperiment_Filtered)

library(limma)

fit <- lmFit(ExpressionData_Filtered, design)
#Adding an eBayes correction to help reduce the influence of outliers/small sample size on estimates
efit <- eBayes(fit, trend=TRUE)
dt<-decideTests(efit)
summary(decideTests(efit))

# (Intercept) SummarizedExperiment_Filtered[[1]]$Stress_factorunpredictable chronic mild stress SummarizedExperiment_Filtered[[1]]$Drug_factorfluoxetine
# Down             0                                                                                 2                                                     1213
# NotSig           0                                                                             14411                                                    12334
# Up           14414                                                                                 1                                                      867
# SummarizedExperiment_Filtered[[1]]$Stress_factorunpredictable chronic mild stress:SummarizedExperiment_Filtered[[1]]$Drug_factorfluoxetine
# Down                                                                                                                                            0
# NotSig                                                                                                                                      14412
# Up                                                                                                                                              2

str(efit)
# Formal class 'MArrayLM' [package "limma"] with 1 slot
# ..@ .Data:List of 23
# .. ..$ : num [1:14414, 1:4] 7.89 4.66 4.66 4.66 4.67 ...
# .. .. ..- attr(*, "dimnames")=List of 2
# .. .. .. ..$ : chr [1:14414] "10003" "10011" "10015" "10016" ...
# .. .. .. ..$ : chr [1:4] "(Intercept)" "SummarizedExperiment_Filtered[[1]]$Stress_factorunpredictable chronic mild stress" "SummarizedExperiment_Filtered[[1]]$Drug_factorfluoxetine" "SummarizedExperiment_Filtered[[1]]$Stress_factorunpredictable chronic mild stress:SummarizedExperiment_Filtered"| __truncated__
# .. ..$ : int 4
# .. ..$ : int [1:4] 0 1 2 3
# .. ..$ :List of 5

write.fit(efit, adjust="BH", file="Limma_results84185.txt")
write.csv(rowData(SummarizedExperiment_Filtered[[1]]), "Annotation_LimmaResultsGSE84185.csv")