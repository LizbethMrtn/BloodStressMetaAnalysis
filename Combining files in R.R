TempResultsJoined<-read.delim("Limma_results72262.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
### Read in annotation file since 72262 didn't have one file with annotation plus actual data
TempResultsJoined_Annotation<-read.csv("Annotation_LimmaResultsGSE72262.csv", header=TRUE, stringsAsFactors = FALSE)
### combine the limma results with the annotation
library(plyr)
temp<-join(TempResultsJoined_Annotation, TempResultsJoined, by="X", type="right")
str(temp)
TempResultsJoined<-temp