#################load packages and data#################
library(psych)
library(ggplot2)
library(qgraph)
library(gplots)
library(corrplot)
library(Hmisc)
library(reshape2)
library(plyr)
library(boot)
library(RColorBrewer)
library(stringr)
options(digits = 3)
theme_set(theme_bw(base_size = 20))


gf<-read.csv("/Users/sakiko/Desktop/19_Komagino_MRS_MAGet/CIVET_202111/ANOVA/residual_harmonized_data2.csv",header=T,sep=",",fileEncoding="utf-8")


for (i in 1:nrow(gf)){
  if (gf[i,"group"]==0){
    gf[i,"group"]<-"HC"
  }else if (gf[i,"group"]==1){
    gf[i,"group"]<-"TnRS"
  }else if (gf[i,"group"]==2){
    gf[i,"group"]<-"TRS"
  }
}

corticalnames<-names(gf)[startsWith(names(gf),"residual_DKT_CT_")]
subset_gf <-gf[c("ID", "group", corticalnames)]

TRSstudy_characteristics<-read.csv("/Users/sakiko/Desktop/19_Komagino_MRS_MAGet/CIVET_202111/Dataset/TRSstudy_characteristics_data_20220607.csv")





###################################################split by groups##########################
splitgroups<-split(subset_gf, subset_gf$group, drop=FALSE)

hconlydata<-splitgroups[["HC"]][1:ncol(subset_gf)]
trsonlydata<-splitgroups[["TRS"]][1:ncol(subset_gf)]
ntrsonlydata<-splitgroups[["TnRS"]][1:ncol(subset_gf)]

hconlydatasubset<-na.omit(hconlydata[3:ncol(subset_gf)])
trsonlydatasubset<-na.omit(trsonlydata[3:ncol(subset_gf)])
ntrsonlydatasubset<-na.omit(ntrsonlydata[3:ncol(subset_gf)])

###HC###
patientdata<-hconlydatasubset
patientdata2<-hconlydata
Newdataset<-hconlydata
PCC_matrix<-cor(Newdataset[3:ncol(subset_gf)])

dataN<-length(corticalnames)

#Pearson correlation coefficience
difZ_result<-data.frame(matrix(rep(NA, nrow(patientdata)*dataN*(dataN-1)/2), ncol=dataN*(dataN-1)/2))

matrixnames<-str_sub(rownames(PCC_matrix),17,-1)

difZ_result_name_matrix<-data.frame(matrix(rep(NA,dataN*dataN),dataN,dataN))
for (i in 1:dataN){
  for (k in  1:dataN){
    difZ_result_name_matrix[i,k]<-paste(matrixnames[k],"_",matrixnames[i],sep="")
  }}
difZ_result_name<-difZ_result_name_matrix[lower.tri(difZ_result_name_matrix)]

names(difZ_result)<-difZ_result_name

for (subject in 1:nrow(patientdata)){
REFDATA<-hconlydatasubset[-subject,]
REF_PCC_matrix<-cor(REFDATA)
difZ_PCC_matrix<-data.frame(matrix(rep(NA,dataN*dataN),dataN,dataN))

for (i in 1:dataN){
  for (j in 1:dataN){
    difZ_PCC_matrix[i,j]<-(PCC_matrix[i,j]-REF_PCC_matrix[i,j])*(nrow(REFDATA)-1)/(1-REF_PCC_matrix[i,j]^2)
  }
}
difZ_result[subject,]<-difZ_PCC_matrix[lower.tri(difZ_PCC_matrix)]
}

difZ_result$ID<-patientdata2$ID
difZ_result$group<-patientdata2$group
HC_difZ_result<-difZ_result



###TRS###
patientdata<-trsonlydatasubset
patientdata2<-trsonlydata
Newdataset<-trsonlydata
PCC_matrix<-cor(Newdataset[3:ncol(subset_gf)])

dataN<-length(corticalnames)
#Pearson correlation coefficience
difZ_result<-data.frame(matrix(rep(NA, nrow(patientdata)*dataN*(dataN-1)/2), ncol=dataN*(dataN-1)/2))

matrixnames<-str_sub(rownames(PCC_matrix),17,-1)

difZ_result_name_matrix<-data.frame(matrix(rep(NA,dataN*dataN),dataN,dataN))
for (i in 1:dataN){
  for (k in  1:dataN){
    difZ_result_name_matrix[i,k]<-paste(matrixnames[k],"_",matrixnames[i],sep="")
  }}
difZ_result_name<-difZ_result_name_matrix[lower.tri(difZ_result_name_matrix)]

names(difZ_result)<-difZ_result_name

for (subject in 1:nrow(patientdata)){
REFDATA<-trsonlydatasubset[-subject,]
REF_PCC_matrix<-cor(REFDATA)
difZ_PCC_matrix<-data.frame(matrix(rep(NA,dataN*dataN),dataN,dataN))

for (i in 1:dataN){
  for (j in 1:dataN){
    difZ_PCC_matrix[i,j]<-(PCC_matrix[i,j]-REF_PCC_matrix[i,j])*(nrow(REFDATA)-1)/(1-REF_PCC_matrix[i,j]^2)
  }
}
difZ_result[subject,]<-difZ_PCC_matrix[lower.tri(difZ_PCC_matrix)]
}

difZ_result$ID<-patientdata2$ID
difZ_result$group<-patientdata2$group
TRS_difZ_result<-difZ_result


###NTRS###
patientdata<-ntrsonlydatasubset
patientdata2<-ntrsonlydata
Newdataset<-ntrsonlydata
PCC_matrix<-cor(Newdataset[3:ncol(subset_gf)])

dataN<-length(corticalnames)
#Pearson correlation coefficience
difZ_result<-data.frame(matrix(rep(NA, nrow(patientdata)*dataN*(dataN-1)/2), ncol=dataN*(dataN-1)/2))

matrixnames<-str_sub(rownames(PCC_matrix),17,-1)

difZ_result_name_matrix<-data.frame(matrix(rep(NA,dataN*dataN),dataN,dataN))
for (i in 1:dataN){
  for (k in  1:dataN){
    difZ_result_name_matrix[i,k]<-paste(matrixnames[k],"_",matrixnames[i],sep="")
  }}
difZ_result_name<-difZ_result_name_matrix[lower.tri(difZ_result_name_matrix)]

names(difZ_result)<-difZ_result_name

for (subject in 1:nrow(patientdata)){
REFDATA<-ntrsonlydatasubset[-subject,]
REF_PCC_matrix<-cor(REFDATA)
difZ_PCC_matrix<-data.frame(matrix(rep(NA,dataN*dataN),dataN,dataN))

for (i in 1:dataN){
  for (j in 1:dataN){
    difZ_PCC_matrix[i,j]<-(PCC_matrix[i,j]-REF_PCC_matrix[i,j])*(nrow(REFDATA)-1)/(1-REF_PCC_matrix[i,j]^2)
  }
}
difZ_result[subject,]<-difZ_PCC_matrix[lower.tri(difZ_PCC_matrix)]
}

difZ_result$ID<-patientdata2$ID
difZ_result$group<-patientdata2$group
NTRS_difZ_result<-difZ_result
Merged_difZ_result<-rbind(HC_difZ_result,TRS_difZ_result,NTRS_difZ_result)

write.csv(Merged_difZ_result,"/Users/sakiko/Desktop/19_Komagino_MRS_MAGet/CIVET_202111/IndividualSCN/Ingroup_Covbat_IDSC_data.csv",row.names = F)
