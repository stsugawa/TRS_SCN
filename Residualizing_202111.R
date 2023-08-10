#未完成

library(corrplot)
library(ggplot2)
library(Hmisc)
library(reshape2)
library(plyr)

options(digits=3)
theme_set(theme_bw(base_size = 20))


####residuals of cortical thickness####

df<-read.csv("/Users/sakiko/Desktop/19_Komagino_MRS_MAGet/CIVET_202111/Dataset/data_all_merge_harmonized.csv",header=T,sep=",",stringsAsFactors=FALSE)

f <- function(name) {
  linmod <- lm(as.formula(paste(name, "~ age + sex")), data = df)
  residu_res<-c()
  for (i in 1:nrow(df)){
    if (is.na(df[i,name])){
      residu_res[i]<-NA
    }
    else{
      residu_res[i]<-df[i,name]- linmod$coefficients[[1]] 
      for (k in 2:length(linmod$coefficients)){
        residu_res[i]<-residu_res[i]-linmod$coefficients[[k]]*df[i,names(linmod$coefficients)[k]]
      }
    }
  }
  return (residu_res) }


results = sapply(thickness_colnames, f)
results<-as.data.frame(results)
for (j in names(results)){
  df[,paste0("residual_",j)]<-results[,j]
}
write.csv(df, file="/Users/sakiko/Desktop/19_Komagino_MRS_MAGet/CIVET_202111/ANOVA/residual_harmonized_data.csv",row.names = F)