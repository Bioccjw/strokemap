library(Seurat)


setwd('....')


#####################################################Add data#########################################################################################################################################################################################################################
load('all.merge_region-result.Rdata')
cluster=as.matrix(all.merge@meta.data$seurat_clusters)
rownames(cluster)=rownames(all.merge@meta.data)  
x=unlist(strsplit(rownames(cluster),split='_'))
region_sample=as.matrix(all.merge@meta.data$orig.ident)
region_sample=cbind(rownames(all.merge@meta.data),region_sample)
region_sample[,1]=x[seq(1,length(x),2)]

load('/pt1-weichuli.Rdata')
pt1$Spatial@scale.data=all.merge$SCT@scale.data[,which(region_sample[,2]=="pt1")]
colnames(pt1$Spatial@scale.data)=region_sample[which(region_sample[,2]=="pt1"),1]
pt1.region=all.merge@meta.data$region[which(region_sample[,2]=="pt1")]
pt1@meta.data=cbind(pt1@meta.data,pt1.region)
colnames(pt1@meta.data)[ncol(pt1@meta.data)]='region'

load('pt2-weichuli.Rdata')
pt2$Spatial@scale.data=all.merge$SCT@scale.data[,which(region_sample[,2]=="pt2")]
colnames(pt2$Spatial@scale.data)=region_sample[which(region_sample[,2]=="pt2"),1]
pt2.region=all.merge@meta.data$region[which(region_sample[,2]=="pt2")]
pt2@meta.data=cbind(pt2@meta.data,pt2.region)
colnames(pt2@meta.data)[ncol(pt2@meta.data)]='region'

#####################################################HVG dropout percentage###########################################################################################################################################################################################################
HVG=VariableFeatures(all.merge)
pt1_other_counts=pt1$Spatial@counts[which(rownames(pt1$Spatial@counts)%in%HVG),which(pt1@meta.data$region!='ICA')]  ##
pt2_other_counts=pt2$Spatial@counts[which(rownames(pt2$Spatial@counts)%in%HVG),which(pt2@meta.data$region!='ICA')]  ##

pt1.dropout_percentage=matrix(0,ncol=1,nrow=length(HVG))
rownames(pt1.dropout_percentage)=rownames(pt1_other_counts)
pt2.dropout_percentage=matrix(0,ncol=1,nrow=length(HVG))
rownames(pt2.dropout_percentage)=rownames(pt2_other_counts)
for(i in 1:length(HVG)){
  pt1.dropout_percentage[i]=length(which(pt1_other_counts[i,]==0))/ncol(pt1_other_counts)
  pt2.dropout_percentage[i]=length(which(pt2_other_counts[i,]==0))/ncol(pt2_other_counts)
}
write.table(pt1.dropout_percentage,sep='\t',quote=F,file='./pt1.dropout_percentage.txt')
write.table(pt2.dropout_percentage,sep='\t',quote=F,file='./pt2.dropout_percentage.txt')

#####################################################Min distance####################################################################################################################################################################################################################
object=c(pt1,pt2)
min_distance <-function(pt_sham.num){
  if(pt_sham.num==1){ 
    pt.Coordinate=object[[pt_sham.num]]@images$pt1@coordinates 
  }
  else{
    pt.Coordinate=object[[pt_sham.num]]@images$pt2@coordinates 
  }
  pt_ICA=pt.Coordinate[which(rownames(pt.Coordinate)%in%rownames(object[[pt_sham.num]]@meta.data)[which(object[[pt_sham.num]]@meta.data$region=='ICA')]),] 
  pt_other=pt.Coordinate[which(rownames(pt.Coordinate)%in%rownames(pt_ICA)==F),]  #
  #remove Outliers 
  if(pt_sham.num==1){
    pt_ICA=pt_ICA[which(rownames(pt_ICA)%in%c('AGAAGAGCGCCGTTCC-1','ACTATTTCCGGGCCCA-1')==F),]
  }  
  else{
    pt_ICA=pt_ICA[which(rownames(pt_ICA)%in%c('TTAGAAGAACATGACT-1','GTTTGGGCTTGTGAGC-1','TGGAAGGATAAAGATG-1','AAGCGCAGGGCTTTGA-1','CCAGATAGTTGAGTGA-1','AGAAGAGCGCCGTTCC-1','AGGTATAATTGATAGT-1','GCTTGCAGCACAATTG-1','CCGCTTGCTGACATGG-1','GTCGTGTCTGGTCATC-1')==F),]

  }
  pt.min.distance=matrix(10000,ncol=2,nrow=nrow(pt_other))
  for(i in 1:nrow(pt_other)){
    x=matrix(10000,ncol=2,nrow=nrow(pt_ICA))
    for(j in 1:nrow(pt_ICA)){
      x[j,1]=dist(rbind(pt_other[i,2:3],pt_ICA[j,2:3]),method='euclidean')
	  x[j,2]=dist(rbind(pt_other[i,4:5],pt_ICA[j,4:5]),method='euclidean')
    } 
    pt.min.distance[i,1]=min(x[,1])
    pt.min.distance[i,2]=min(x[,2])
  }
  rownames(pt.min.distance)=rownames(pt_other)
  colnames(pt.min.distance)=c('x_y.min_distance','xiangsu.min_distance')
  return(pt.min.distance)
}
pt1.min_distance=min_distance(pt_sham.num=1)
cor.test(x=pt1.min_distance[,1],y=pt1.min_distance[,2],method="pearson")  
#cor=0.9906336
write.table(pt1.min_distance,sep='\t',quote=F,file='./pt1.min_distance.txt')

pt2.min_distance=min_distance(pt_sham.num=2)
cor.test(x=pt2.min_distance[,1],y=pt2.min_distance[,2],method="pearson")  
#cor=0.9604831 
write.table(pt2.min_distance,sep='\t',quote=F,file='./pt2.min_distance.txt')

#####################################################Correlation between gene expression and distance of PT##########################################################################################################################################################################
object=c(pt1,pt2)
cor_exp.dis <-function(pt_sham.num){
  pt_other_exp=object[[pt_sham.num]]$Spatial@scale.data[,which(object[[pt_sham.num]]@meta.data$region!='ICA')]
  pt_distance=read.table(sep='\t',header=T,file=paste0('./pt',pt_sham.num,'.min_distance.txt'))
  pt_distance=as.matrix(pt_distance)
  pearson.p=matrix(NA,ncol=1,nrow=nrow(pt_other_exp))
  pearson.cor=matrix(NA,ncol=1,nrow=nrow(pt_other_exp))
  if(all(rownames(pt_distance)==colnames(pt_other_exp))){
    for(i in 1:nrow(pt_other_exp)){
     x=cor.test(x=pt_distance[,1],y=pt_other_exp[i,],method="pearson") 
     pearson.p[i]=x$p.value
     pearson.cor[i]=x$estimate
   }
   pt_result=cbind(pearson.p,pearson.cor)            
   colnames(pt_result)=c('p-value','cor-value')
   rownames(pt_result)=rownames(pt_other_exp)
  }
  else{
    print('error')
  }
  return(pt_result)
}
pt1_result=cor_exp.dis(pt_sham.num=1)
write.table(pt1_result,sep='\t',quote=F,file='./pt1.allHVG_cor.txt')

pt1.dropout_percentage=read.table(sep='\t',header=T,file='pt1.dropout_percentage.txt')
gene=rownames(pt1.dropout_percentage)[which(pt1.dropout_percentage<=0.7)]
re=pt1_result[which(rownames(pt1_result)%in%gene),]
write.table(re,file='./pt1.cor_drop out=0.7.txt',quote=F,sep='\t')

pt2_result=cor_exp.dis(pt_sham.num=2)
write.table(pt1_result,sep='\t',quote=F,file='./pt2_all.txt')

pt2.dropout_percentage=read.table(sep='\t',header=T,file='./pt2.dropout_percentage.txt')
gene=rownames(pt2.dropout_percentage)[which(pt2.dropout_percentage<=0.7)]
re=pt2_result[which(rownames(pt2_result)%in%gene),]
write.table(re,file='./pt2.cor_drop out=0.7.txt',quote=F,sep='\t')

#####################################################average distance of brain regions################################################################################################################################################################################################
distance_mean<-function(seurat_object,all_distance){
  allregion=sort(unique(seurat_object@meta.data$region))
  result=matrix(NA,nrow=length(allregion),ncol=1)
  rownames(result)=allregion
  for(i in 1:length(allregion)){
   result[i]=mean(all_distance[which(rownames(all_distance)%in%rownames(seurat_object@meta.data)[which(seurat_object@meta.data$region==allregion[i])]),1])
   }
  return(result)
}
pt1.region=unique(pt1@meta.data$region)
pt1.region=pt1.region[-which(pt1.region=='ICA')]
pt1.min_distance=read.table(sep='\t',header=T,file='pt1.min_distance.txt')
pt1.region_mean_distance=distance_mean(seurat_object=pt1,all_distance=pt1.min_distance)  
pt1.region_mean_distance=pt1.region_mean_distance[-which(rownames(pt1.region_mean_distance)=='ICA'),]
write.table(pt1.region_mean_distance,sep='\t',quote=F,file='./pt1.region_mean_distance.txt')

pt2.region=unique(pt2@meta.data$region)
pt2.region=pt2.region[-which(pt2.region=='ICA')]
pt2.min_distance=read.table(sep='\t',header=T,file='pt2.min_distance.txt')
pt2.region_mean_distance=distance_mean(seurat_object=pt2,all_distance=pt2.min_distance)  
pt2.region_mean_distance=pt2.region_mean_distance[-which(rownames(pt2.region_mean_distance)=='ICA'),]
write.table(pt2.region_mean_distance,sep='\t',quote=F,file='./pt2.region_mean_distance.txt')

#####################################################Self-validation##################################################################################################################################################################################################################
distance_pred=function(seurat_object,pt.region,pt.sham.num,pt.region_mean_distance){
  pt.distance_pred=matrix(0,ncol=1,nrow=length(which(seurat_object@meta.data$region%in%pt.region)))
  rownames(pt.distance_pred)=rownames(seurat_object@meta.data)[which(seurat_object@meta.data$region%in%pt.region)]
  colnames(pt.distance_pred)=paste0('pt',pt.sham.num,'_pred_distance')
  for(i in 1:length(pt.region)){
    pt.region_spot=rownames(seurat_object@meta.data)[which(seurat_object@meta.data$region%in%pt.region[i])] 										   
    pt.distance_pred[which(rownames(pt.distance_pred)%in%pt.region_spot)]=as.numeric(pt.region_mean_distance[which(rownames(pt.region_mean_distance)==pt.region[i]),1])   
  }  
  return(pt.distance_pred)
}

cor_exp.dis=function(pt_other_exp,seurat_object,pt.region,pt.distance_pred){
  pearson.p=matrix(NA,ncol=1,nrow=nrow(pt_other_exp))
  pearson.cor=matrix(NA,ncol=1,nrow=nrow(pt_other_exp))
  pt.region_exp=seurat_object@assays$Spatial@scale.data[,which(seurat_object@meta.data$region%in%pt.region)]

  if(all(rownames(pt.distance_pred)==colnames(pt.region_exp))==T){
   for(i in 1:nrow(pt_other_exp)){
    x=cor.test(x=pt.distance_pred[,1],y=pt.region_exp[which(rownames(pt.region_exp)==rownames(pt_other_exp)[i]),],method="pearson") 
    pearson.p[i]=x$p.value
    pearson.cor[i]=x$estimate
   }
  pt.pearson.result.pred=cbind(pearson.p,pearson.cor)            #
  colnames(pt.pearson.result.pred)=c('p-value','cor-value')
  rownames(pt.pearson.result.pred)=rownames(pt_other_exp)
  return(pt.pearson.result.pred)
  }
  else{
    print('error')
  }
}
###Calculating the predicted distance
pt1.region_mean_distance=read.table(sep='\t',header=T,file='pt1.region_mean_distance.txt')
pt1.region_mean_distance=as.matrix(pt1.region_mean_distance)
pt1.region=rownames(pt1.region_mean_distance)
pt1.distance_pred=distance_pred(seurat_object=pt1,pt.region=pt1.region,pt.sham.num=1,pt.region_mean_distance=pt1.region_mean_distance)
write.table(pt1.distance_pred,sep='\t',quote=F,file='./pt1.distance_pred.txt')

pt2.region_mean_distance=read.table(sep='\t',header=T,file='pt2.region_mean_distance.txt')
pt2.region_mean_distance=as.matrix(pt2.region_mean_distance)
pt2.region=rownames(pt2.region_mean_distance)
pt2.distance_pred=distance_pred(seurat_object=pt2,pt.region=pt2.region,pt.sham.num=1,pt.region_mean_distance=pt2.region_mean_distance)
write.table(pt2.distance_pred,sep='\t',quote=F,file='./pt2.distance_pred.txt')

###Calculating the correlation coefficient
pt1_other_exp=pt1$Spatial@scale.data[,which(pt1@meta.data$region!='ICA')]  
pt1.pearson_result=read.table(sep='\t',header=T,file='pt1.cor_drop out=0.7.txt')
pt1_other_exp=pt1_other_exp[which(rownames(pt1_other_exp)%in%rownames(pt1.pearson_result)),]  
pt1.pearson_result_pred=cor_exp.dis(pt_other_exp=pt1_other_exp,seurat_object=pt1,pt.region=pt1.region,pt.distance_pred=pt1.distance_pred) 
#colnames(pt1.pearson_result_pred) "p.value"   "cor.value" 
#colnames(pt1.pearson_result) "p.value"   "cor.value"  
cor.test(pt1.pearson_result[,2],pt1.pearson_result_pred[,2],method='pearson')  
#cor=0.9960667 

pt2_other_exp=pt2$Spatial@scale.data[,which(pt2@meta.data$region!='ICA')]  
pt2.pearson_result=read.table(sep='\t',header=T,file='pt2.cor_drop out=0.7.txt')
pt2_other_exp=pt2_other_exp[which(rownames(pt2_other_exp)%in%rownames(pt2.pearson_result)),]  
pt2.pearson_result_pred=cor_exp.dis(pt_other_exp=pt2_other_exp,seurat_object=pt2,pt.region=pt2.region,pt.distance_pred=pt2.distance_pred) 
#colnames(pt2.pearson_result_pred) "p.value"   "cor.value" 
#colnames(pt2.pearson_result) "p.value"   "cor.value"  
cor.test(pt2.pearson_result[,2],pt2.pearson_result_pred[,2],method='pearson')  
#cor=0.9918282  

#####################################################Cross-validation################################################################################################################################################################################################################
distance_pred <-function(intersect.region,observe_object,observe.num,pt.region_mean_distance){
  pt.distance_pred=matrix(0,ncol=1,nrow=length(which(observe_object@meta.data$region%in%intersect.region)))
  rownames(pt.distance_pred)=rownames(observe_object@meta.data)[which(observe_object@meta.data$region%in%intersect.region)]
  colnames(pt.distance_pred)=paste0('pt',observe.num,'_pred_distance')
  for(i in 1:length(intersect.region)){
    intersect.region_spot=rownames(observe_object@meta.data)[which(observe_object@meta.data$region%in%intersect.region[i])]     
    pt.distance_pred[which(rownames(pt.distance_pred)%in%intersect.region_spot),1]=pt.region_mean_distance[which(rownames(pt.region_mean_distance)==intersect.region[i]),1]
  }
  return(pt.distance_pred)
}

cor_exp.dis <-function(pt_other_exp,intersect.region,observe_object,pt.distance_pred){
  pearson.p=matrix(NA,ncol=1,nrow=nrow(pt_other_exp))
  pearson.cor=matrix(NA,ncol=1,nrow=nrow(pt_other_exp))
  pt_intersect.region_exp=observe_object@assays$Spatial@scale.data[,which(observe_object@meta.data$region%in%intersect.region)]
  if(all(rownames(pt.distance_pred)==colnames(pt_intersect.region_exp))){
    for(i in 1:nrow(pt_other_exp)){
     x=cor.test(x=as.numeric(pt.distance_pred[,1]),y=pt_intersect.region_exp[which(rownames(pt_intersect.region_exp)==rownames(pt_other_exp)[i]),],method="pearson") 
     pearson.p[i]=x$p.value
     pearson.cor[i]=x$estimate
    }
    pt.pearson.result.pred=cbind(pearson.p,pearson.cor)            #
    colnames(pt.pearson.result.pred)=c('p-value','cor-value')
    rownames(pt.pearson.result.pred)=rownames(pt_other_exp)
    return(pt.pearson.result.pred)
  }
  else{
    print('error')
  }
}
###Calculating the predicted distance
pt1.region_mean_distance=read.table(sep='\t',header=T,file='pt1.region_mean_distance.txt')
pt1.region_mean_distance=as.matrix(pt1.region_mean_distance)
pt1.region=rownames(pt1.region_mean_distance)

pt2.region_mean_distance=read.table(sep='\t',header=T,file='pt2.region_mean_distance.txt')
pt2.region_mean_distance=as.matrix(pt2.region_mean_distance)
pt2.region=rownames(pt2.region_mean_distance)
intersect.region=intersect(pt1.region,pt2.region)

pt1_pt2.distance_pred=distance_pred(intersect.region=intersect.region,observe_object=pt2,observe.num=2,pt.region_mean_distance=pt1.region_mean_distance) 
write.table(pt1_pt2.distance_pred,sep='\t',quote=F,file='./pt1_pt2.distance_pred.txt')
pt2_pt1.distance_pred=distance_pred(intersect.region=intersect.region,observe_object=pt1,observe.num=1,pt.region_mean_distance=pt2.region_mean_distance) 
write.table(pt2_pt1.distance_pred,sep='\t',quote=F,file='./pt2_pt1.distance_pred.txt')

###Calculating the correlation coefficient
pt2_other_exp=pt2$Spatial@scale.data[,which(pt2@meta.data$region!='ICA')]  
pt2.pearson_result=read.table(sep='\t',header=T,file='pt2.cor_drop out=0.7.txt')
pt2_other_exp=pt2_other_exp[which(rownames(pt2_other_exp)%in%rownames(pt2.pearson_result)),]  
pt1_pt2.pearson_result_pred=cor_exp.dis(pt_other_exp=pt2_other_exp,intersect.region=intersect.region,observe_object=pt2,pt.distance_pred=pt1_pt2.distance_pred)
#colnames(pt2.pearson_result):"p.value"   "cor.value"
#colnames(pt1_pt2.pearson_result_pred):"p.value"   "cor.value"   
cor.test(pt2.pearson_result[,2],pt1_pt2.pearson_result_pred[,2],method='pearson')  
#cor=0.9669309

pt1_other_exp=pt1$Spatial@scale.data[,which(pt1@meta.data$region!='ICA')]  
pt1.pearson_result=read.table(sep='\t',header=T,file='pt1.cor_drop out=0.7.txt')
pt1_other_exp=pt1_other_exp[which(rownames(pt1_other_exp)%in%rownames(pt1.pearson_result)),]  
pt2_pt1.pearson_result_pred=cor_exp.dis(pt_other_exp=pt1_other_exp,intersect.region=intersect.region,observe_object=pt1,pt.distance_pred=pt2_pt1.distance_pred)
#colnames(pt1.pearson_result):"p.value"   "cor.value"
#colnames(pt2_pt1.pearson_result_pred):"p.value"   "cor.value"   
cor.test(pt1.pearson_result[,2],pt2_pt1.pearson_result_pred[,2],method='pearson')  
#cor=0.9862425 

#####################################################Correlation between gene expression and distance of Sham########################################################################################################################################################################
load('sham1-weichuli.Rdata')
sham1$Spatial@scale.data=all.merge$SCT@scale.data[,which(region_sample[,2]=="sham1")]
colnames(sham1$Spatial@scale.data)=region_sample[which(region_sample[,2]=="sham1"),1]
sham1.region=all.merge@meta.data$region[which(region_sample[,2]=="sham1")]
sham1@meta.data=cbind(sham1@meta.data,sham1.region)
colnames(sham1@meta.data)[ncol(sham1@meta.data)]='region'

load('sham2-weichuli.Rdata')
sham2$Spatial@scale.data=all.merge$SCT@scale.data[,which(region_sample[,2]=="sham2")]
colnames(sham2$Spatial@scale.data)=region_sample[which(region_sample[,2]=="sham2"),1]
sham2.region=all.merge@meta.data$region[which(region_sample[,2]=="sham2")]
sham2@meta.data=cbind(sham2@meta.data,sham2.region)
colnames(sham2@meta.data)[ncol(sham2@meta.data)]='region'

sham_distance_pred <-function(pt.region,sham,sham.num,pt.region_mean_distance){
  intersect.region=intersect(pt.region,sham@meta.data$region)   
  sham.distance_pred=matrix(0,ncol=1,nrow=length(which(sham@meta.data$region%in%intersect.region)))   
  rownames(sham.distance_pred)=rownames(sham@meta.data)[which(sham@meta.data$region%in%intersect.region)]
  colnames(sham.distance_pred)=paste0('sham',sham.num,'_pred_distance')
  for(i in 1:length(intersect.region)){
   intersect.region_spot=rownames(sham@meta.data)[which(sham@meta.data$region%in%intersect.region[i])]  										   
   sham.distance_pred[which(rownames(sham.distance_pred)%in%intersect.region_spot)]=as.numeric(pt.region_mean_distance[which(rownames(pt.region_mean_distance)==intersect.region[i]),1])
  }
  return(sham.distance_pred)
}

cor_exp.dis <-function(pt.pearson.result,sham,pt.region,sham.distance_pred){
  pearson.p=matrix(NA,ncol=1,nrow=nrow(pt.pearson.result))
  pearson.cor=matrix(NA,ncol=1,nrow=nrow(pt.pearson.result))
  intersect.region=intersect(pt.region,sham@meta.data$region)
  sham_intersect.region_exp=sham$Spatial@scale.data[,which(sham@meta.data$region%in%intersect.region)]
  if(all(rownames(sham.distance_pred)==colnames(sham_intersect.region_exp))==T){
    for(i in 1:nrow(pt.pearson.result)){
     x=cor.test(x=sham.distance_pred[,1],y=sham_intersect.region_exp[which(rownames(sham_intersect.region_exp)==rownames(pt.pearson.result)[i]),],method="pearson") 
     pearson.p[i]=x$p.value
     pearson.cor[i]=x$estimate
    }
    sham.pearson.result.pred=cbind(pearson.p,pearson.cor)            #
    colnames(sham.pearson.result.pred)=c('p-value','cor-value')
    rownames(sham.pearson.result.pred)=rownames(pt.pearson.result)
	return(sham.pearson.result.pred)
  }
  else{
    print('error')
  }
}
###Calculating the predicted distance
pt1.region_mean_distance=read.table(sep='\t',header=T,file='pt1.region_mean_distance.txt')
pt1.region_mean_distance=as.matrix(pt1.region_mean_distance)
pt1.region=rownames(pt1.region_mean_distance)
allregion=pt1.region[-which(pt1.region%in%c('PIA_D','PIA_P'))]  
sham1.distance_pred=sham_distance_pred(pt.region=allregion,sham=sham1,sham.num=1,pt.region_mean_distance=pt1.region_mean_distance)
write.table(sham1.distance_pred,sep='\t',quote=F,file='./sham1.distance_pred.txt')

pt2.region_mean_distance=read.table(sep='\t',header=T,file='pt2.region_mean_distance.txt')
pt2.region_mean_distance=as.matrix(pt2.region_mean_distance)
pt2.region=rownames(pt2.region_mean_distance)
allregion=pt2.region[-which(pt2.region%in%c('PIA_D','PIA_P'))]  
sham2.distance_pred=sham_distance_pred(pt.region=allregion,sham=sham2,sham.num=2,pt.region_mean_distance=pt2.region_mean_distance)
write.table(sham2.distance_pred,sep='\t',quote=F,file='./sham2.distance_pred.txt')

###Calculating the correlation coefficient
allregion=pt1.region[-which(pt1.region%in%c('PIA_D','PIA_P'))] 
pt1.pearson_result=read.table(sep='\t',header=T,file='pt1.cor_drop out=0.7.txt')
sham1.pearson_result_pred=cor_exp.dis(pt.pearson.result=pt1.pearson_result,sham=sham1,pt.region=allregion,sham.distance_pred=sham1.distance_pred)
write.table(sham1.pearson_result_pred,sep='\t',quote=F,file='./sham1.pearson_result_pred.txt')

allregion=pt2.region[-which(pt2.region%in%c('PIA_D','PIA_P'))]  
pt2.pearson_result=read.table(sep='\t',header=T,file='pt2.cor_drop out=0.7.txt')
sham2.pearson_result_pred=cor_exp.dis(pt.pearson.result=pt2.pearson_result,sham=sham2,pt.region=allregion,sham.distance_pred=sham2.distance_pred)
write.table(sham2.pearson_result_pred,sep='\t',quote=F,file='./sham2.pearson_result_pred.txt')




