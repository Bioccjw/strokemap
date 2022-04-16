library(Seurat)


setwd('....')

#####################################################Cell mapping#######################################################################################
load('all.merge_cluster-result.Rdata')
part.merge=subset(all.merge[,which(all.merge@meta.data$region%in%c('PIA_P','PIA_D'))])
load('Stroke.rename.rds')
Stroke2=SCTransform(Stroke,do.scale=T,do.center=T)
Stroke2=RunPCA(Stroke2,verbose = FALSE)
anterset_SCT <- FindTransferAnchors(reference = Stroke2,query =part.merge ,normalization.method="SCT") 
predictions.assay_SCT <- TransferData(anchorset = anterset_SCT, refdata = Stroke2@active.ident, prediction.assay = TRUE)
save(predictions.assay_SCT,file='Mapping_result.Rdata')


