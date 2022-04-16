library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)
library(RColorBrewer)


setwd("....")


######################################Initial Data##########################################
###Loading data
Sham1.exp=Read10X_h5(filename ="Sham1-1/2.3.h5_files/filtered_feature_bc_matrix.h5")
Sham2.exp=Read10X_h5(filename ="Sham1-2/2.3.h5_files/filtered_feature_bc_matrix.h5")
Sham3.exp=Read10X_h5(filename ="Sham1-3/2.3.h5_files/filtered_feature_bc_matrix.h5")
Sham4.exp=Read10X_h5(filename ="Sham1-4/2.3.h5_files/filtered_feature_bc_matrix.h5")
PT1.exp=Read10X_h5(filename ="PT1-1/2.3.h5_files/filtered_feature_bc_matrix.h5")
PT2.exp=Read10X_h5(filename ="PT1-2/2.3.h5_files/filtered_feature_bc_matrix.h5")
PT3.exp=Read10X_h5(filename ="PT1-3/2.3.h5_files/filtered_feature_bc_matrix.h5")
PT4.exp=Read10X_h5(filename ="PT1-4/2.3.h5_files/filtered_feature_bc_matrix.h5")
dim(Sham1.exp)
#32285  1752
dim(Sham2.exp)
#32285  2031
dim(Sham3.exp)
#32285  2618
dim(Sham4.exp)
#32285  2793
dim(PT1.exp)
#32285  2397
dim(PT2.exp)
#32285  2650
dim(PT3.exp)
#32285  2781
dim(PT4.exp)
#32285  2755
sham1<- CreateSeuratObject(counts = Sham1.exp,project = 'sham1', assay = 'Spatial')						 
sham.img <- Read10X_Image(image.dir = 'Sham1-1/spatial sham1')
DefaultAssay(object =sham.img) <- 'Spatial'
img<-sham.img[colnames(x = sham1)]
sham1[['sham1']] <- img
sham2<- CreateSeuratObject(counts = Sham2.exp,project = 'sham2', assay = 'Spatial')						 
sham.img <- Read10X_Image(image.dir = 'Sham1-2/spatial sham2')
DefaultAssay(object =sham.img) <- 'Spatial'
img <-sham.img[colnames(x = sham2)]
sham2[['sham2']] <- img
sham3<- CreateSeuratObject(counts = Sham3.exp,project = 'sham3', assay = 'Spatial')						 
sham.img <- Read10X_Image(image.dir = 'Sham1-3/spatial sham3')
DefaultAssay(object =sham.img) <- 'Spatial'
img <-sham.img[colnames(x = sham3)]
sham3[['sham3']] <- img
sham4<- CreateSeuratObject(counts = Sham4.exp,project = 'sham4', assay = 'Spatial')						 
sham.img <- Read10X_Image(image.dir = 'Sham1-4/spatial sham4')
DefaultAssay(object =sham.img) <- 'Spatial'
img <-sham.img[colnames(x = sham4)]
sham4[['sham4']] <- img
pt1<- CreateSeuratObject(counts = PT1.exp,project = 'pt1', assay = 'Spatial')						 
pt.img <- Read10X_Image(image.dir = 'PT1-1/spatial PT1')
DefaultAssay(object =pt.img) <- 'Spatial'
img <-pt.img[colnames(x = pt1)]
pt1[['pt1']] <- img
pt2<- CreateSeuratObject(counts = PT2.exp,project = 'pt2', assay = 'Spatial')						 
pt.img <- Read10X_Image(image.dir = 'PT1-2/spatial PT2')
DefaultAssay(object =pt.img) <- 'Spatial'
img <-pt.img[colnames(x = pt2)]
pt2[['pt2']] <- img
pt3<- CreateSeuratObject(counts = PT3.exp,project = 'pt3', assay = 'Spatial')						 
pt.img <- Read10X_Image(image.dir = 'PT1-3/spatial PT3')
DefaultAssay(object =pt.img) <- 'Spatial'
img <-pt.img[colnames(x = pt3)]
pt3[['pt3']] <- img
pt4<- CreateSeuratObject(counts = PT4.exp,project = 'pt4', assay = 'Spatial')						 
pt.img <- Read10X_Image(image.dir = 'PT1-4/spatial PT4')
DefaultAssay(object =pt.img) <- 'Spatial'
img <-pt.img[colnames(x = pt4)]
pt4[['pt4']] <- img
dim(sham1)
#32285  1752
dim(sham2)
#32285  2031
dim(sham3)
#32285  2618
dim(sham4)
#32285  2793
dim(pt1)
#32285  2397
dim(pt2)
#32285  2650
dim(pt3)
#32285  2781
dim(pt4)
#32285  2755
head(sham1@meta.data)
head(sham2@meta.data)
head(sham3@meta.data)
head(sham4@meta.data)
head(pt1@meta.data)
head(pt2@meta.data)
head(pt3@meta.data)
head(pt4@meta.data)

save(sham1,file='./sham1-weichuli.Rdata')
save(sham2,file='./sham2-weichuli.Rdata')
save(sham3,file='./sham3-weichuli.Rdata')
save(sham4,file='./sham4-weichuli.Rdata')
save(pt1,file='./pt1-weichuli.Rdata')
save(pt2,file='./pt2-weichuli.Rdata')
save(pt3,file='./pt3-weichuli.Rdata')
save(pt4,file='./pt4-weichuli.Rdata')

######################################QC,Normalization,Variable#########################################################################
all.merge=merge(sham1,c(sham2,sham3,sham4,pt1,pt2,pt3,pt4))
dim(all.merge)
#32285 19777
all.merge=SCTransform(all.merge,assay = 'Spatial',new.assay.name = "SCT",variable.features.n = 6000,do.scale = T,do.center = T)
dim(all.merge)
#21736 19777
all.merge <- RunPCA(all.merge, verbose = FALSE)  
all.merge@meta.data$stim <- c(rep("sham", ncol(sham1)), rep("sham", ncol(sham2)),rep("sham", ncol(sham3)),
                              rep("sham", ncol(sham4)),rep("pt", ncol(pt1)),rep("pt", ncol(pt2)),
							  rep("pt", ncol(pt3)),rep("pt", ncol(pt4))
							  )
all.merge <- RunHarmony(all.merge,region.by.vars="stim", plot_convergence = TRUE,assay.use="SCT") 

######################################Clustering########################################################################################
all.merge <- FindNeighbors(all.merge,reduction ="harmony" ,dim=1:50,k.param=14)
all.merge <- FindClusters(all.merge, verbose = FALSE,resolution=1.6)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 
col_vector =col_vector [-c(4,5,6,7,8)]

sham1.cluster=as.numeric(all.merge@meta.data$seurat_clusters[which(all.merge@meta.data$orig.ident=='sham1')])
sham1@meta.data=cbind(sham1@meta.data,sham1.cluster)
pdf(file='./sham1_cluster.pdf')
SpatialPlot(sham1,region.by = "sham1.cluster",cols=col_vector[sort(unique(sham1.cluster))])
dev.off()

sham2.cluster=as.numeric(all.merge@meta.data$seurat_clusters[which(all.merge@meta.data$orig.ident=='sham2')])
sham2@meta.data=cbind(sham2@meta.data,sham2.cluster)
pdf(file='./sham2_cluster.pdf')
SpatialPlot(sham2,region.by = "sham2.cluster",cols=col_vector[sort(unique(sham2.cluster))])
dev.off()

sham3.cluster=as.numeric(all.merge@meta.data$seurat_clusters[which(all.merge@meta.data$orig.ident=='sham3')])
sham3@meta.data=cbind(sham3@meta.data,sham3.cluster)
pdf(file='./sham3_cluster.pdf')
SpatialPlot(sham3,region.by = "sham3.cluster",cols=col_vector[sort(unique(sham3.cluster))])
dev.off()

sham4.cluster=as.numeric(all.merge@meta.data$seurat_clusters[which(all.merge@meta.data$orig.ident=='sham4')])
sham4@meta.data=cbind(sham4@meta.data,sham4.cluster)
pdf(file='./sham4_cluster.pdf')
SpatialPlot(sham4,region.by = "sham4.cluster",cols=col_vector[sort(unique(sham4.cluster))])
dev.off()

pt1.cluster=as.numeric(all.merge@meta.data$seurat_clusters[which(all.merge@meta.data$orig.ident=='pt1')])
pt1@meta.data=cbind(pt1@meta.data,pt1.cluster)
pdf(file='./pt1_cluster.pdf')
SpatialPlot(pt1,region.by = "pt1.cluster",cols=col_vector[sort(unique(pt1.cluster))])
dev.off()

pt2.cluster=as.numeric(all.merge@meta.data$seurat_clusters[which(all.merge@meta.data$orig.ident=='pt2')])
pt2@meta.data=cbind(pt2@meta.data,pt2.cluster)
pdf(file='./pt2_cluster.pdf')
SpatialPlot(pt2,region.by = "pt2.cluster",cols=col_vector[sort(unique(pt2.cluster))])
dev.off()

pt3.cluster=as.numeric(all.merge@meta.data$seurat_clusters[which(all.merge@meta.data$orig.ident=='pt3')])
pt3@meta.data=cbind(pt3@meta.data,pt3.cluster)
pdf(file='./pt3_cluster.pdf')
SpatialPlot(pt3,region.by = "pt3.cluster",cols=col_vector[sort(unique(pt3.cluster))])
dev.off()

pt4.cluster=as.numeric(all.merge@meta.data$seurat_clusters[which(all.merge@meta.data$orig.ident=='pt4')])
pt4@meta.data=cbind(pt4@meta.data,pt4.cluster)
pdf(file='./pt4_cluster.pdf')
SpatialPlot(pt4,region.by = "pt4.cluster",cols=col_vector[sort(unique(pt4.cluster))])
dev.off()

######################################Annotating brain regions########################################################################################
region=as.character(all.merge@meta.data$seurat_clusters)
region[which(region%in%c('10'))]='ICA'
region[which(region%in%c('13'))]='PIA_P'
region[which(region%in%c('14'))]='PIA_D'

region[which(region%in%c('17','27'))]='CTX_L1'
region[which(region%in%c('1','31'))]='CTX_L2/3'
region[which(region%in%c('7'))]='CTX_L4'
region[which(region%in%c('5','18'))]='CTX_L5'
region[which(region%in%c('4'))]='CTX_L6'
region[which(region%in%c('29'))]='CTX_L6b'

region[which(region%in%c('26'))]='PirCTX_L1'
region[which(region%in%c('19'))]='PirCTX_L2'
region[which(region%in%c('20'))]='PirCTX_L3'

region[which(region%in%c('15','8'))]='WM'

region[which(region%in%c('23'))]='STR_LSX'
region[which(region%in%c('9'))]='STRv'
region[which(region%in%c('3'))]='STR_CP_1'
region[which(region%in%c('35','0'))]='STR_CP_2'

region[which(region%in%c('24'))]='PAL_1'
region[which(region%in%c('11'))]='PAL_2'
region[which(region%in%c('39'))]='PAL_3'

region[which(region%in%c('12'))]='BST&sAMY'

region[which(region%in%c('21'))]='CTXsp_CLA&EP'
region[which(region%in%c('32'))]='CTXsp_BLA'

region[which(region%in%c('16'))]='HIP_DG'
region[which(region%in%c('37'))]='HIP_CA1'
region[which(region%in%c('33'))]='HIP_CA2/3'

region[which(region%in%c('28'))]='TH_RT'
region[which(region%in%c('36'))]='TH_EPI'
region[which(region%in%c('22'))]='TH_ILM&MTN'
region[which(region%in%c('2'))]='TH_MED&VENT'

region[which(region%in%c('38'))]='HY_MEZ'
region[which(region%in%c('6'))]='HY_LZ'
region[which(region%in%c('25'))]='HY_ARH'

region[which(region%in%c('30','34'))]='VL'

all.merge@meta.data=cbind(all.merge@meta.data,region)
Idents(all.merge)=all.merge@meta.data$region
######################################Saving###############################################
save.image(file = "./all.merge_cluster-result.Rdata")