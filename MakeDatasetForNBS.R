
combatData<-read.csv("/Users/sakiko/Desktop/19_Komagino_MRS_MAGet/CIVET_202111/Dataset/data_all_merge_harmonized.csv",header=T,sep=",")

names(combatData)<-c("group","age","sex",corticalnames,"Study")
combatData<-combatData[,c("group","sex","age",corticalnames)]
ATRS0HC1<-subset(combatData,combatData$group==0|combatData$group==2)
for(i in 1:nrow(ATRS0HC1)){
  if (ATRS0HC1$group[i]==0){
    ATRS0HC1$group[i]<-1
  }else if (ATRS0HC1$group[i]==2){
    ATRS0HC1$group[i]<-0
  }
}
ATnRS0HC1<-subset(combatData,combatData$group==0|combatData$group==1)
for(i in 1:nrow(ATnRS0HC1)){
  if (ATnRS0HC1$group[i]==0){
    ATnRS0HC1$group[i]<-1
  }else if (ATnRS0HC1$group[i]==1){
    ATnRS0HC1$group[i]<-0
  }
}

ATRS0TnRS1<-subset(combatData,combatData$group==1|combatData$group==2)
for(i in 1:nrow(ATRS0TnRS1)){
  if (ATRS0TnRS1$group[i]==1){
    ATRS0TnRS1$group[i]<-1
  }else if (ATRS0TnRS1$group[i]==2){
    ATRS0TnRS1$group[i]<-0
  }
}

write.csv(ATRS0HC1,"/Users/sakiko/Desktop/19_Komagino_MRS_MAGet/CIVET_202111/Dataset/All_TRS0HC1.csv",row.names = F)
write.csv(ATnRS0HC1,"/Users/sakiko/Desktop/19_Komagino_MRS_MAGet/CIVET_202111/Dataset/All_TnRS0HC1.csv",row.names = F)
write.csv(ATRS0TnRS1,"/Users/sakiko/Desktop/19_Komagino_MRS_MAGet/CIVET_202111/Dataset/All_TRS0TnRS1.csv",row.names = F)

