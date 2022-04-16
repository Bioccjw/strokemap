library(Seurat)


setwd("....")


#####################################################DEGs##################################################################################################################################
load(file='all.merge_cluster-result.Rdata')
naoqu=unique(all.merge@meta.data$region)  
naoqu=naoqu[-which(naoqu%in%c('HY_MEZ','ICA','PIA_P','PIA_D','VL'))]  ##HY_MEZ only appears on sham4
for(i in 1:length(naoqu)){
    markers <- FindMarkers(all.merge,ident.1 ='pt',ident.2 ='sham',group.by = 'stim',subset.ident = naoqu[i],only.pos=F,logfc.threshold=log(1.5))
	if(i%in%grep('/',naoqu,value=F)){
	  write.table(markers,sep='\t',quote=F,file=paste0('./',gsub('/','_',naoqu[i]),'_DEGs.txt'))
	}
	else{
	  write.table(markers,sep='\t',quote=F,file=paste0('./',naoqu[i],'_DEGs.txt'))
	}
}

ISCM=c('ICA','PIA_P','PIA_D')
for(i in 1:length(ISCM)){
  print(ISCM[i])
  all_region=all.merge@meta.data$region
  all_region=as.matrix(all_region)
  rownames(all_region)=rownames(all.merge@meta.data)
  all_region[which(all_region[,1]%in%ISCM&all.merge@meta.data$orig.ident%in%c('sham1','sham2','sham3','sham4')),1]='error'		
  all_region[which(all_region[,1]%in%c('CTX_L1','CTX_L2/3','CTX_L4','CTX_L5','CTX_L6','CTX_L6b')&all.merge@meta.data$orig.ident%in%c('sham1','sham2','sham3','sham4')),1]=ISCM[i]		
  all.merge@meta.data=cbind(all.merge@meta.data,all_region)
  Idents(all.merge)=all_region
  markers=FindMarkers(all.merge,ident.1 ='pt',ident.2 ='sham', group.by = 'stim', subset.ident = ISCM[i],
  logfc.threshold = log(1.5),only.pos=F)  
  write.table(markers,sep='\t',quote=F,file=paste0('./',ISCM[i],'_DEGs.txt'))
}


