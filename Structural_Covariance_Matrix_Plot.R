#script made by Eric Plitman and Sakiko Tsugawa

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
options(digits = 3)
theme_set(theme_bw(base_size = 20))

cortical_colnames<-c(
  #frontal
  "SuperiorFrontalGyrus.L",
  "RostralMiddleFrontal.L",	
  "CaudalMiddleFrontal.L",
  "LateralFrontalOpercularis.L",
  "LateralFrontalTriangularis.L", 
  "LateralFrontalOrbitalis.L",
  "LateralOrbitofrontal.L",	
  "MedialOrbitofrontal.L",
  "PrecentralGyrus.L",
  "ParacentralGyrus.L",
  
  #parietal
  "SuperiorParietal.L",
  "InferiorParietal.L",
  "SupramarginalGyrus.L",
  "PostcentralGyrus.L",
  "Precuneus.L",
  
  "Insula.L",
  
  #temporal
  "SuperiorTemporal.L",	
  "MiddleTemporal.L",	
  "InferiorTemporal.L",
  "FusiformGyrus.L",
  "TransverseTemporal.L",
  "EntorhinalCortex.L",
  "Parahippocampal.L",
  
  #occipital
  "InferiorOccipitalCortex.L",
  "LingualGyrus.L",
  "Cuneus.L",	
  "Pericalcarine.L",
  
  #cingulate
  "RostralAnteriorCingulate.L",
  "CaudalAnteriorCingulate.L",
  "PosteriorCingulate.L",	
  "IsthmusCingulateGyrus.L",
  
  #right
  #frontal
  "SuperiorFrontalGyrus.R",
  "RostralMiddleFrontal.R",	
  "CaudalMiddleFrontal.R",
  "LateralFrontalOpercularis.R",
  "LateralFrontalTriangularis.R", 
  "LateralFrontalOrbitalis.R",
  "LateralOrbitofrontal.R",	
  "MedialOrbitofrontal.R",
  "PrecentralGyrus.R",
  "ParacentralGyrus.R",
  
  #parietal
  "SuperiorParietal.R",
  "InferiorParietal.R",
  "SupramarginalGyrus.R",
  "PostcentralGyrus.R",
  "Precuneus.R",
  
  "Insula.R",
  
  #temporal
  "SuperiorTemporal.R",	
  "MiddleTemporal.R",	
  "InferiorTemporal.R",
  "FusiformGyrus.R",
  "TransverseTemporal.R",
  "EntorhinalCortex.R",
  "Parahippocampal.R",
  
  #occipital
  "InferiorOccipitalCortex.R",
  "LingualGyrus.R",
  "Cuneus.R",	
  "Pericalcarine.R",
  
  #cingulate
  "RostralAnteriorCingulate.R",
  "CaudalAnteriorCingulate.R",
  "PosteriorCingulate.R",	
  "IsthmusCingulateGyrus.R"
)





####Komagino dataについて

groupname<-c("Komagino","Toronto","Shimofusa","Combat")
eachstudy<-read.csv("/Users/sakiko/Desktop/19_Komagino_MRS_MAGet/CIVET_202111/ANOVA/residual_data_all_merge.csv",header=T,sep=",",fileEncoding="utf-8")
combatdata<-read.csv("/Users/sakiko/Desktop/19_Komagino_MRS_MAGet/CIVET_202111/ANOVA/residual_harmonized_data2.csv",header=T,sep=",",stringsAsFactors=FALSE)
eachstudy<-na.omit(eachstudy)
GN<-4  ###ここを変える
if (GN ==4){
  gf<-combatdata
  corticalnames<-paste0("DKT_CT_",cortical_colnames)
  for (row in 1:nrow(gf)){
    if (  gf$group[row]==0){
      gf$group[row]<-"HC"
    }else    if (  gf$group[row]==1){
      gf$group[row]<-"TnRS"
    }else    if (  gf$group[row]==2){
      gf$group[row]<-"TRS"
    }
  }

}else {
 gf<-subset(eachstudy,eachstudy$Study==GN)
 corticalnames<-paste0("residual_DKT_CT_",cortical_colnames)
}

subset_gf <-gf[c("ID", "group", corticalnames)]

###################################################split by groups##########################
splitgroups<-split(subset_gf, subset_gf$group, drop=FALSE)
trsonlydata<-splitgroups[["TRS"]][1:ncol(subset_gf)]
tnrsonlydata<-splitgroups[["TnRS"]][1:ncol(subset_gf)]
hconlydata<-splitgroups[["HC"]][1:ncol(subset_gf)]
trsonlydatasubset<-na.omit(trsonlydata[3:ncol(subset_gf)])
tnrsonlydatasubset<-na.omit(tnrsonlydata[3:ncol(subset_gf)])
hconlydatasubset<-na.omit(hconlydata[3:ncol(subset_gf)])

dataN<-length(corticalnames)

SC_HC<-cor(hconlydatasubset)
SC_TnRS<-cor(tnrsonlydatasubset)
SC_TRS<-cor(trsonlydatasubset)

####################################################ST code####################

setwd("/Users/sakiko/Desktop/19_Komagino_MRS_MAGet/CIVET_202111/Visualization/SCN")
pdf(paste0(groupname[GN],"_cor_HC_CT.pdf"))
heatmap.2(SC_HC,Rowv=FALSE,Colv=FALSE,dendrogram="none", col=rev(brewer.pal(11,"RdYlBu")),
          scale="none", density.info="none", trace="none",margins = c(10,10))
dev.off()

pdf(paste0(groupname[GN],"_cor_nTRS_CT.pdf"))
heatmap.2(SC_TnRS,Rowv=FALSE,Colv=FALSE,dendrogram="none", col=rev(brewer.pal(11,"RdYlBu")),
          scale="none", density.info="none", trace="none",margins = c(10,10))
dev.off()

pdf(paste0(groupname[GN],"_cor_TRS_CT.pdf"))
heatmap.2(SC_TRS,Rowv=FALSE,Colv=FALSE,dendrogram="none", col=rev(brewer.pal(11,"RdYlBu")),
          scale="none", density.info="none", trace="none",margins = c(10,10))
dev.off()

#なぜ自分でZを求めるか？→r.testだとzの絶対値しかでないから。

#Fisher Z 相関比較
ntrsn<-length(tnrsonlydatasubset[,1])
hcn<-length(hconlydatasubset[,1])
trsn<-length(trsonlydatasubset[,1])

hc_trs_scn<-data.frame(matrix(rep(NA, dataN*dataN), nrow=dataN))
hc_ntrs_scn<-data.frame(matrix(rep(NA, dataN*dataN), nrow=dataN))
ntrs_trs_scn<-data.frame(matrix(rep(NA, dataN*dataN), nrow=dataN))


for (i in  1:dataN){
  for (k in  1:dataN){
    hc_trs_scn[k,i]<-(fisherz(SC_TRS[k,i])-fisherz(SC_HC[k,i]))/sqrt(1/(trsn-3)+1/(hcn-3))
  }
}
for (i in  1:dataN){
  for (k in  1:dataN){
    hc_ntrs_scn[k,i]<-(fisherz(SC_TnRS[k,i])-fisherz(SC_HC[k,i]))/sqrt(1/(ntrsn-3)+1/(hcn-3))
  }
}

for (i in  1:dataN){
  for (k in  1:dataN){
    ntrs_trs_scn[k,i]<-(fisherz(SC_TRS[k,i])-fisherz(SC_TnRS[k,i]))/sqrt(1/(trsn-3)+1/(ntrsn-3))
  }
}

hc_trs_scn<-as.matrix(hc_trs_scn)
hc_ntrs_scn<-as.matrix(hc_ntrs_scn)
ntrs_trs_scn<-as.matrix(ntrs_trs_scn)

rownames(hc_trs_scn)<-cortical_colnames
colnames(hc_trs_scn)<-cortical_colnames
rownames(hc_ntrs_scn)<-cortical_colnames
colnames(hc_ntrs_scn)<-cortical_colnames

rownames(ntrs_trs_scn)<-cortical_colnames
colnames(ntrs_trs_scn)<-cortical_colnames

pdf(paste0(groupname[GN],"_hc_ntrs_scn.pdf"),onefile=F,width = 1200/72,height = 1200/72)
heatmap.2(hc_ntrs_scn,Rowv=FALSE,Colv=FALSE,dendrogram="none", col=rev(brewer.pal(11,"RdYlBu")),
          scale="none", density.info="none", trace="none",margins = c(10,10))
dev.off()
pdf(paste0(groupname[GN],"_hc_trs_scn.pdf"),onefile=F,width = 1200/72,height = 1200/72)
heatmap.2(hc_trs_scn,Rowv=FALSE,Colv=FALSE,dendrogram="none", col=rev(brewer.pal(11,"RdYlBu")),
          scale="none", density.info="none", trace="none",margins = c(10,10))
dev.off()
pdf(paste0(groupname[GN],"_ntrs_trs_scn.pdf"),onefile=F,width = 1200/72,height = 1200/72)
heatmap.2(ntrs_trs_scn,Rowv=FALSE,Colv=FALSE,dendrogram="none", col=rev(brewer.pal(11,"RdYlBu")),
          scale="none", density.info="none", trace="none",margins = c(10,10))
dev.off()

###Z valueのヒストグラム作成
histData<-as.data.frame(fisherz(SC_HC)[upper.tri(fisherz(SC_HC))])
names(histData)<-"HC"
histData$TnRS<-fisherz(SC_TnRS)[upper.tri(fisherz(SC_TnRS))]
histData$TRS<-fisherz(SC_TRS)[upper.tri(fisherz(SC_TRS))]
df <- melt(histData) 
g <- ggplot(df, aes(x = value, y = ..density.., fill = variable))
g <- g + geom_histogram(position = "identity", alpha = 0.8,binwidth = 0.025)
g <- g + geom_density(aes(color = variable, alpha = 0.2), show.legend = F)
pdf(paste0(groupname[GN],"_Histgram_structural_covariance.pdf"),onefile=F,width = 1200/72,height = 400/72)
plot(g)
dev.off()

df <- histData
g <- ggplot(df, aes(x = HC, y = ..density..))
g <- g + geom_histogram(position = "identity", alpha = 0.8,binwidth = 0.025)
g <- g + geom_density(aes( alpha = 0.2), show.legend = F)
g1<-g

g <- ggplot(df, aes(x = TnRS, y = ..density..))
g <- g + geom_histogram(position = "identity", alpha = 0.8,binwidth = 0.025)
g <- g + geom_density(aes( alpha = 0.2), show.legend = F)
g2<-g

g <- ggplot(df, aes(x = TRS, y = ..density..))
g <- g + geom_histogram(position = "identity", alpha = 0.8,binwidth = 0.025)
g <- g + geom_density(aes( alpha = 0.2), show.legend = F)
g3<-g

g1 <- ggplotGrob(g1)
g2 <- ggplotGrob(g2)
g3 <- ggplotGrob(g3)
g <- rbind(g1,g2, g3, size = "first")
g$widths = grid::unit.pmax(g1$widths,g2$widths, g3$widths)

GP<-ggplot(stack(histData), aes(x=values)) + geom_histogram(binwidth = 0.025) + facet_wrap(~ind, ncol=1)

pdf(paste0(groupname[GN],"_Each_Histgram_structural_covariance.pdf"),onefile=F,width = 1200/72,height = 800/72)
plot(GP)
dev.off()

# rownames(hc_ntrs_scn)[1]
hc_trs_scn.pvalue<-data.frame(matrix(rep(NA, 4*dataN*(dataN-1)/2), nrow=dataN*(dataN-1)/2))
#x!
hc_ntrs_scn.pvalue<-data.frame(matrix(rep(NA, 4*dataN*(dataN-1)/2), nrow=dataN*(dataN-1)/2))
ntrs_trs_scn.pvalue<-data.frame(matrix(rep(NA, 4*dataN*(dataN-1)/2), nrow=dataN*(dataN-1)/2))

for (i in 1:(dataN-1)){
  for (k in (i+1):dataN){
    hc_trs_scn.pvalue[(k-i+dataN*(dataN-1)/2-(dataN-(i-1))*((dataN-(i-1))-1)/2),1]<-paste(rownames(hc_trs_scn)[i],"_",colnames(hc_trs_scn)[k],sep="")
    hc_trs_scn.pvalue[(k-i+dataN*(dataN-1)/2-(dataN-(i-1))*((dataN-(i-1))-1)/2),2]<-hc_trs_scn[i,k]
    hc_ntrs_scn.pvalue[(k-i+dataN*(dataN-1)/2-(dataN-(i-1))*((dataN-(i-1))-1)/2),1]<-paste(rownames(hc_ntrs_scn)[i],"_",colnames(hc_ntrs_scn)[k],sep="")
    hc_ntrs_scn.pvalue[(k-i+dataN*(dataN-1)/2-(dataN-(i-1))*((dataN-(i-1))-1)/2),2]<-hc_ntrs_scn[i,k]
    ntrs_trs_scn.pvalue[(k-i+dataN*(dataN-1)/2-(dataN-(i-1))*((dataN-(i-1))-1)/2),1]<-paste(rownames(ntrs_trs_scn)[i],"_",colnames(ntrs_trs_scn)[k],sep="")
    ntrs_trs_scn.pvalue[(k-i+dataN*(dataN-1)/2-(dataN-(i-1))*((dataN-(i-1))-1)/2),2]<-ntrs_trs_scn[i,k]    

  }
}
x<-paste(rownames(hc_ntrs_scn)[1],"_",colnames(hc_ntrs_scn)[2],sep="")
names(hc_trs_scn.pvalue)<-c("Region", "HC_TRS_Z", "HC_TRS_P", "HC_TRS_FDR")
names(hc_ntrs_scn.pvalue)<-c("Region", "HC_TnRS_Z", "HC_TnRS_P", "HC_TnRS_FDR")
names(ntrs_trs_scn.pvalue)<-c("Region", "TnRS_TRS_Z", "TnRS_TRS_P", "TnRS_TRS_FDR")

for (i in 1:nrow(hc_trs_scn.pvalue)){
  
  hc_trs_scn.pvalue[i,3]<- 2*pnorm(-abs(hc_trs_scn.pvalue[i,2]))
  hc_ntrs_scn.pvalue[i,3]<- 2*pnorm(-abs(hc_ntrs_scn.pvalue[i,2]))
  
  ntrs_trs_scn.pvalue[i,3]<- 2*pnorm(-abs(ntrs_trs_scn.pvalue[i,2]))
  
}

hc_trs_scn.pvalue[,4]<-p.adjust(hc_trs_scn.pvalue[,3],method = "fdr")
hc_ntrs_scn.pvalue[,4]<-p.adjust(hc_ntrs_scn.pvalue[,3],method = "fdr")
ntrs_trs_scn.pvalue[,4]<-p.adjust(ntrs_trs_scn.pvalue[,3],method = "fdr")

result<-merge(hc_trs_scn.pvalue,hc_ntrs_scn.pvalue,by="Region",all=T)
result<-merge(result,ntrs_trs_scn.pvalue,by="Region")
if (GN!=2){
  write.csv(result,paste0(groupname[GN],"_SCN_results.csv"))
}


map_trshc_ntrshc<-data.frame(matrix(rep(NA, dataN*dataN), nrow=dataN))
map_trshc_ntrshc
for (i in  1:dataN){
  for (k in  1:dataN){
    if (i<k){
      map_trshc_ntrshc[k,i]<-hc_trs_scn[k,i]
    }
    else if (i==k){
      map_trshc_ntrshc[k,i]<-0
    }
    else if (i>k){
      map_trshc_ntrshc[k,i]<-hc_ntrs_scn[k,i]
    }
  }
}
map_trshc_ntrshc
rownames(map_trshc_ntrshc)<-cortical_colnames
colnames(map_trshc_ntrshc)<-cortical_colnames
map_trshc_ntrshc<-as.matrix(map_trshc_ntrshc)



if (GN==2){



utrsonlydata<-splitgroups[["TRS"]][1:ncol(subset_gf)]
utrsonlydatasubset<-na.omit(utrsonlydata[3:ncol(subset_gf)])
SC_UTRS<-cor(trsonlydatasubset)

####################################################ST code####################

setwd("/Users/sakiko/Desktop/19_Komagino_MRS_MAGet/CIVET_202111/Visualization/SCN")

pdf(paste0(groupname[GN],"_cor_UTRS_CT.pdf"))
heatmap.2(SC_UTRS,Rowv=FALSE,Colv=FALSE,dendrogram="none", col=rev(brewer.pal(11,"RdYlBu")),
          scale="none", density.info="none", trace="none",margins = c(10,10))
dev.off()

#なぜ自分でZを求めるか？→r.testだとzの絶対値しかでないから。



#Fisher Z 相関比較

utrsn<-length(trsonlydatasubset[,1])

hc_utrs_scn<-data.frame(matrix(rep(NA, dataN*dataN), nrow=dataN))
ntrs_utrs_scn<-data.frame(matrix(rep(NA, dataN*dataN), nrow=dataN))
trs_utrs_scn<-data.frame(matrix(rep(NA, dataN*dataN), nrow=dataN))


for (i in  1:dataN){
  for (k in  1:dataN){
    hc_utrs_scn[k,i]<-(fisherz(SC_TRS[k,i])-fisherz(SC_HC[k,i]))/sqrt(1/(trsn-3)+1/(hcn-3))
  }
}
for (i in  1:dataN){
  for (k in  1:dataN){
    ntrs_utrs_scn[k,i]<-(fisherz(SC_TnRS[k,i])-fisherz(SC_HC[k,i]))/sqrt(1/(ntrsn-3)+1/(hcn-3))
  }
}

for (i in  1:dataN){
  for (k in  1:dataN){
    trs_utrs_scn[k,i]<-(fisherz(SC_TRS[k,i])-fisherz(SC_TnRS[k,i]))/sqrt(1/(trsn-3)+1/(ntrsn-3))
  }
}

hc_utrs_scn<-as.matrix(hc_utrs_scn)
ntrs_utrs_scn<-as.matrix(ntrs_utrs_scn)
trs_utrs_scn<-as.matrix(trs_utrs_scn)

rownames(hc_utrs_scn)<-cortical_colnames
colnames(hc_utrs_scn)<-cortical_colnames
rownames(ntrs_utrs_scn)<-cortical_colnames
colnames(ntrs_utrs_scn)<-cortical_colnames

rownames(trs_utrs_scn)<-cortical_colnames
colnames(trs_utrs_scn)<-cortical_colnames

pdf(paste0(groupname[GN],"_hc_utrs_scn.pdf"),onefile=F,width = 1200/72,height = 1200/72)
heatmap.2(hc_utrs_scn,Rowv=FALSE,Colv=FALSE,dendrogram="none", col=rev(brewer.pal(11,"RdYlBu")),
          scale="none", density.info="none", trace="none",margins = c(10,10))
dev.off()
pdf(paste0(groupname[GN],"_ntrs_utrs_scn.pdf"),onefile=F,width = 1200/72,height = 1200/72)
heatmap.2(ntrs_utrs_scn,Rowv=FALSE,Colv=FALSE,dendrogram="none", col=rev(brewer.pal(11,"RdYlBu")),
          scale="none", density.info="none", trace="none",margins = c(10,10))
dev.off()
pdf(paste0(groupname[GN],"_trs_utrs_scn.pdf"),onefile=F,width = 1200/72,height = 1200/72)
heatmap.2(trs_utrs_scn,Rowv=FALSE,Colv=FALSE,dendrogram="none", col=rev(brewer.pal(11,"RdYlBu")),
          scale="none", density.info="none", trace="none",margins = c(10,10))
dev.off()


# rownames(hc_ntrs_scn)[1]

hc_utrs_scn.pvalue<-data.frame(matrix(rep(NA, 4*dataN*(dataN-1)/2), nrow=dataN*(dataN-1)/2))
#x!
ntrs_utrs_scn.pvalue<-data.frame(matrix(rep(NA, 4*dataN*(dataN-1)/2), nrow=dataN*(dataN-1)/2))
trs_utrs_scn.pvalue<-data.frame(matrix(rep(NA, 4*dataN*(dataN-1)/2), nrow=dataN*(dataN-1)/2))

for (i in 1:(dataN-1)){
  for (k in (i+1):dataN){#dataN*(i-1)-i*(i+1)/2+k-i
    hc_utrs_scn.pvalue[(k-i+dataN*(dataN-1)/2-(dataN-(i-1))*((dataN-(i-1))-1)/2),1]<-paste(rownames(hc_trs_scn)[i],"_",colnames(hc_trs_scn)[k],sep="")
    hc_utrs_scn.pvalue[(k-i+dataN*(dataN-1)/2-(dataN-(i-1))*((dataN-(i-1))-1)/2),2]<-hc_trs_scn[i,k]
    ntrs_utrs_scn.pvalue[(k-i+dataN*(dataN-1)/2-(dataN-(i-1))*((dataN-(i-1))-1)/2),1]<-paste(rownames(hc_ntrs_scn)[i],"_",colnames(hc_ntrs_scn)[k],sep="")
    ntrs_utrs_scn.pvalue[(k-i+dataN*(dataN-1)/2-(dataN-(i-1))*((dataN-(i-1))-1)/2),2]<-hc_ntrs_scn[i,k]
    trs_utrs_scn.pvalue[(k-i+dataN*(dataN-1)/2-(dataN-(i-1))*((dataN-(i-1))-1)/2),1]<-paste(rownames(ntrs_trs_scn)[i],"_",colnames(ntrs_trs_scn)[k],sep="")
    trs_utrs_scn.pvalue[(k-i+dataN*(dataN-1)/2-(dataN-(i-1))*((dataN-(i-1))-1)/2),2]<-ntrs_trs_scn[i,k]    
    
  }
}
x<-paste(rownames(hc_utrs_scn)[1],"_",colnames(hc_utrs_scn)[2],sep="")

names(hc_utrs_scn.pvalue)<-c("Region", "HC_UTRS_Z", "HC_UTRS_P", "HC_UTRS_FDR")
names(trs_utrs_scn.pvalue)<-c("Region", "TRS_UTRS_Z", "TRS_UTRS_P", "TRS_UTRS_FDR")
names(ntrs_utrs_scn.pvalue)<-c("Region", "TnRS_UTRS_Z", "TnRS_UTRS_P", "TnRS_UTRS_FDR")

for (i in 1:nrow(hc_utrs_scn.pvalue)){
  
  hc_utrs_scn.pvalue[i,3]<- 2*pnorm(-abs(hc_utrs_scn.pvalue[i,2]))
  ntrs_utrs_scn.pvalue[i,3]<- 2*pnorm(-abs(ntrs_utrs_scn.pvalue[i,2]))
  
  trs_utrs_scn.pvalue[i,3]<- 2*pnorm(-abs(trs_utrs_scn.pvalue[i,2]))
  
}

hc_utrs_scn.pvalue[,4]<-p.adjust(hc_utrs_scn.pvalue[,3],method = "fdr")
ntrs_utrs_scn.pvalue[,4]<-p.adjust(ntrs_utrs_scn.pvalue[,3],method = "fdr")
trs_utrs_scn.pvalue[,4]<-p.adjust(trs_utrs_scn.pvalue[,3],method = "fdr")

result<-merge(result,hc_utrs_scn.pvalue,by="Region")
result<-merge(result,ntrs_utrs_scn.pvalue,by="Region")
result<-merge(result,trs_utrs_scn.pvalue,by="Region")
write.csv(result,paste0(groupname[GN],"_SCN_results.csv"))
}
