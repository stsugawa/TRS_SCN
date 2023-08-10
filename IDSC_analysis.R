library(ggplot2)
library(Hmisc)
library(reshape2)
library(plyr)
library(metafor)


parseLm<- function (model) 
{
  modelcoef <- summary(model)[[4]]
  rsquared <- summary(model)[[9]]
  fstatistic <- as.numeric(summary(model)[[10]][1])
  df<- model$df.residual
  rows <- nrow(modelcoef)
  cols <- ncol(modelcoef)
  
  betacolnames= paste("beta_", rownames(modelcoef)[2:rows], sep="")
  secolnames= paste("se_", rownames(modelcoef)[2:rows], sep="")
  tcolnames= paste("t_", rownames(modelcoef)[2:rows], sep="")
  pcolnames= paste("p_", rownames(modelcoef)[2:rows], sep="")
  
  output=cbind(fstatistic, df, rsquared)
  
  a<-c(1,2,3,4)
  for (i in a){
    for (j in 2:rows) {
      output= cbind(output, modelcoef[j,i])
    }
  }
  
  colnames(output) <- c("F_statistic","df","R_squared", betacolnames, secolnames, tcolnames, pcolnames)
  return(as.data.frame(output))
}

TRSstudy_characteristics<-read.csv("/Users/sakiko/Desktop/19_Komagino_MRS_MAGet/CIVET_202111/Dataset/TRSstudy_characteristics_data_20220607.csv",header=T,sep=",",fileEncoding="utf-8")
gf<-read.csv("/Users/sakiko/Desktop/19_Komagino_MRS_MAGet/CIVET_202111/IndividualSCN/Ingroup_Covbat_IDSC_data.csv",header=T,sep=",",fileEncoding="utf-8")
combatdf<-read.csv("/Users/sakiko/Desktop/19_Komagino_MRS_MAGet/CIVET_202111/ANOVA/residual_harmonized_data2.csv",header=T,sep=",",fileEncoding="utf-8")
GF<-merge(TRSstudy_characteristics[,names(TRSstudy_characteristics)!="group"],gf,by="ID",all.y=T)
GF<-merge(GF,combatdf[,names(combatdf)!="group"],by="ID",all.x=T)

##DOIとmeanIDSC

NBSfinalresult<-read.csv("/Users/sakiko/Desktop/19_Komagino_MRS_MAGet/CIVET_202111/SCN_result/Combat_SCN_finalresult_202111.csv")

TRSHC<-subset(NBSfinalresult,!is.na(NBSfinalresult$TRSHC))
TnRSHC<-subset(NBSfinalresult,!is.na(NBSfinalresult$TnRSHC))
TRSHC_NBSedge<-rep(NA,nrow(TRSHC))
TnRSHC_NBSedge<-rep(NA,nrow(TnRSHC))
for (i in 1:nrow(TRSHC)){
  TRSHC_NBSedge[i]<-paste0(TRSHC$ROI1[i],"_",TRSHC$ROI2[i])
}
for (i in 1:nrow(TnRSHC)){
  TnRSHC_NBSedge[i]<-paste0(TnRSHC$ROI1[i],"_",TnRSHC$ROI2[i])
}
TnRSHC_NBSedge

###Calculate average IDSC value of the significantly different network between TnRS and HC
GF$Mean_TnRSHC_IDSC<-apply(GF[,TnRSHC_NBSedge], 1, mean)

chara_col<-c("PANSS_total","PANSS_positive","PANSS_general","PANSS_negative","CPZ","duration_of_disease")

###Make correlation plot
for (i in chara_col){ 
  for (j in c("Mean_TnRSHC_IDSC")){
    g<-ggplot(GF,aes_string(x=i,y=j,col="group"))
    g <- g + geom_point()
    g <- g + geom_smooth(method = "lm")
    pdf(paste0("/Users/sakiko/Desktop/19_Komagino_MRS_MAGet/CIVET_202111/IndividualSCN/Ingroup_Correlation_",j,"_",i,".pdf"),onefile=F,width = 400/72,height = 300/72)
    plot(g)
    dev.off()
    
  }
}

#### Regression analysis
datalist<-list()
datalist[[1]]<-GF
datalist[[2]]<-subset(GF,GF$group=="TnRS")
datalist[[3]]<-subset(GF,GF$group=="TRS")

datalistname<-c("IDSC","TnRS","TRS")
xname<-chara_col
yname<-c("Mean_TnRSHC_IDSC") 

all_regression_result<-data.frame()
for (k in 1:length(datalist)){#
  regression_result<-data.frame()
  gf<-datalist[[k]]
  for (j in yname) {
    final.data_1<- data.frame()
    for (i in 1:length(xname)){
      if (all(is.na(gf[,xname[i]])) | all(is.na(gf[,j]))){
        output<-as.data.frame(t(c(j,rep(NA,(3+3*2)))))
        output[2:(1+3+3*2)]<-as.numeric(output[2:(1+3+3*2)])
      }else{
        model <-lm(as.formula(paste(j," ~",xname[i])), data=gf)
        ###Store coefficients
        coefs<- parseLm(model)
        output=cbind(xname[i],round(coefs[,1:3],2),round(coefs[,startsWith(names(coefs),"t")],2),coefs[,startsWith(names(coefs),"p")])
        names(output)<-c("xname","F_statistic","GF","R_squared", "t_xname" ,"p_xname")
      }
      final.data_1 <- rbind.fill(final.data_1, output)
      
    }
    names(final.data_1)<-paste(j,names(final.data_1))
    if (j==yname[1]){
      regression_result<-final.data_1
    }else{
      regression_result<-cbind(regression_result,final.data_1)
    }
  }
  if (k==1){
    all_regression_result<-regression_result
  }else{
    all_regression_result<-rbind(all_regression_result,regression_result)
  }
}

rownames(all_regression_result)<-as.character(paste(rep(datalistname,each=length(xname)),1:length(xname)))
for (i in names(all_regression_result)[endsWith(names(all_regression_result),"p_xname")]){
  all_regression_result[,paste(i,"FDR")]<-round(p.adjust(all_regression_result[,i], method = "BH"),3)
}

Patient_all_regression_result<-  all_regression_result

write.csv(Patient_all_regression_result,paste0("/Users/sakiko/Desktop/19_Komagino_MRS_MAGet/CIVET_202111/IndividualSCN/Ingroup_Correlation_IDSC_clinicalvalue.csv"),row.names = T)

