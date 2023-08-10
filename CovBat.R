#library(devtools)
#devtools::install_github("andy1764/CovBat_Harmonization/R")
#if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("genefilter")
#install_github("jfortin1/neuroCombat_Rpackage")
library(CovBat)
library(neuroCombat)
library(genefilter)
Dataset<-read.csv("/Users/sakiko/Desktop/19_Komagino_MRS_MAGet/CIVET_202111/CovBat/data_all_merge_combat.csv",header=T,sep=",",fileEncoding="utf-8")
Dataset$sex<-as.factor(Dataset$sex)
Dataset$group<-as.factor(Dataset$group)
Mod <- model.matrix(~Dataset$age+Dataset$sex+Dataset$group)
colnames<-names(Dataset)[startsWith(names(Dataset),"DKT")]
dat<-as.matrix(t(Dataset[,colnames]))
bat<-as.factor(t(Dataset$Study))
CovBatData<-covbat(dat, bat, mod = Mod,percent.var = 0.95, n.pc = NULL,
                   train = NULL, mean.only = FALSE, std.var = TRUE, 
                   resid = FALSE, eb = TRUE, parametric = TRUE,
                   score.eb = FALSE, score.parametric = TRUE, verbose = FALSE)
Output<-cbind(Dataset[,c("group","age","sex")],t(CovBatData$dat.covbat),Dataset[,"Study"])
names(Output)<-names(Dataset)
write.csv(Output,"/Users/sakiko/Desktop/19_Komagino_MRS_MAGet/CIVET_202111/CovBat/data_all_merge_harmonized.csv",row.names = F)                   

