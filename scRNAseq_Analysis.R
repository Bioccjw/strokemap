################## remotes::install_github('chris-mcginnis-ucsf/DoubletFinder') ### install doubletfinder
#####remove potential doublets##############################################################
rm(list=ls())
library(dplyr)
library(cowplot)
library(Seurat)
library(DoubletFinder)
library(ggplot2) # 
library(Matrix)
library(ggrepel)
doublet_rate<-0.016
########################
PT.data <- Read10X(data.dir = "./PT/")
PT.data <- CreateSeuratObject(counts = PT.data, project = "PT", min.cells = 3)  ##18377  6796
load_number<-dim(PT.data)
PT.data[["percent.mt"]] <- PercentageFeatureSet(PT.data, pattern = "^mt-")
PT.data <- subset(PT.data, subset = nFeature_RNA > 200  & percent.mt < 20)  ##18377  6559
MT_number<-dim(PT.data)
PT.data <- NormalizeData(PT.data)
PT.data <- FindVariableFeatures(PT.data, selection.method = "vst", nfeatures = 2000)
PT.data <- ScaleData(PT.data)
PT.data <- RunPCA(PT.data)
PT.data <- RunUMAP(PT.data, dims = 1:10)
PT.data <- RunTSNE(PT.data, dims = 1:10)
#####pre-process####################
sweep.res.list <- paramSweep_v3(PT.data, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
######## Homotypic Doublet Proportion Estimate    ##ref:  A Spatiotemporal Organ-Wide Gene Expression and Cell Atlas of the Developing Human Heart
nExp_poi <- round(doublet_rate*ncol(PT.data))        ## Assuming 1.6% doublet 
PT.data <- doubletFinder_v3(PT.data, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
pdf("Fig.S4A.PT.pdf")
DimPlot(PT.data,pt.size = 1,label=TRUE, label.size = 5,reduction = "tsne",group.by = "DF.classifications_0.25_0.3_105" )+theme(aspect.ratio = 1)
dev.off()
PT.data<- subset(PT.data, DF.classifications_0.25_0.005_105=="Singlet" )  ##18377  6231
DoubletFinder_number<-dim(PT.data)

write.table(as.matrix(GetAssayData(object = PT.data, slot = "counts")), 
            'PT.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)
Sham.data <- Read10X(data.dir = "./Sham/")
Sham.data <- CreateSeuratObject(counts = Sham.data, project = "Sham", min.cells = 3) 
load_number<-dim(Sham.data)
Sham.data[["percent.mt"]] <- PercentageFeatureSet(Sham.data, pattern = "^mt-")
Sham.data <- subset(Sham.data, subset = nFeature_RNA > 200  & percent.mt < 20)  

MT_number<-dim(Sham.data)
Sham.data <- NormalizeData(Sham.data)
Sham.data <- FindVariableFeatures(Sham.data, selection.method = "vst", nfeatures = 2000)
Sham.data <- ScaleData(Sham.data)
Sham.data <- RunPCA(Sham.data)
Sham.data <- RunUMAP(Sham.data, dims = 1:10)
Sham.data <- RunTSNE(Sham.data, dims = 1:10)
#####pre-process#############################################################################
sweep.res.list <- paramSweep_v3(Sham.data, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
######## Homotypic Doublet Proportion Estimate    ##ref:  A Spatiotemporal Organ-Wide Gene Expression and Cell Atlas of the Developing Human Heart
nExp_poi <- round(doublet_rate*ncol(Sham.data))        ## Assuming 1.6% doublet 
Sham.data <- doubletFinder_v3(Sham.data, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
# table(Sham.data$DF.classifications_0.25_0.005_328)
pdf("Fig.S4A.Sham.pdf")
DimPlot(Sham.data,pt.size = 1,label=TRUE, label.size = 5,reduction = "tsne",group.by = "DF.classifications_0.25_0.07_137" )+theme(aspect.ratio = 1)
dev.off()
Sham.data<- subset(Sham.data, DF.classifications_0.25_0.15_137=="Singlet" )  ##18377  6231
DoubletFinder_number<-dim(Sham.data)

##############
write.table(as.matrix(GetAssayData(object = Sham.data, slot = "counts")), 
            'Sham.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)
##########process data with seurat###################################################################################
rm(list=ls())
library(dplyr)
library(cowplot)
library(Seurat)
resolution<-2  ## set the resolution
####load the expression matrix of the doubletfinder################ 
PT.sparse<-read.csv("./PT/PT.csv")  ## read CSV
Sham.sparse<-read.csv("./Sham/Sham.csv")

common_genes<-intersect(rownames(PT.sparse),rownames(Sham.sparse))
PT_specific_genes<-setdiff(rownames(PT.sparse),rownames(Sham.sparse))
Sham_specific_genes<-setdiff(rownames(Sham.sparse),rownames(PT.sparse))
####merge the two expression matrix########################
PT.sparse_common<-PT.sparse[common_genes,]
PT.sparse_specific<-PT.sparse[PT_specific_genes,]
PT.sparse_zero<-matrix(0,ncol=ncol(PT.sparse),nrow=length(Sham_specific_genes))
colnames(PT.sparse_zero)<-colnames(PT.sparse_common)
rownames(PT.sparse_zero)<-Sham_specific_genes
PT.sparse<-rbind(PT.sparse_common,PT.sparse_specific,PT.sparse_zero)
###############################
Sham.sparse_common<-Sham.sparse[common_genes,]
Sham.sparse_specific<-Sham.sparse[Sham_specific_genes,]
Sham.sparse_zero<-matrix(0,ncol=ncol(Sham.sparse),nrow=length(PT_specific_genes))
colnames(Sham.sparse_zero)<-colnames(Sham.sparse_common)
rownames(Sham.sparse_zero)<-PT_specific_genes
Sham.sparse<-rbind(Sham.sparse_common,Sham.sparse_zero,Sham.sparse_specific)
##############################
colnames(Sham.sparse)<- sub(pattern = "1", replacement = "2", x = colnames(Sham.sparse))##unique colname
####create the seurat object####################################
PT <- CreateSeuratObject(counts = cbind(PT.sparse, Sham.sparse), project = "PT", min.cells = 5) %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(pc.genes = PT@var.genes, npcs = 50, verbose = FALSE)

PT@meta.data$stim <- c(rep("PT", ncol(PT.sparse)), rep("Sham", ncol(Sham.sparse)))#

#######################dimension reduction##############
PT <- PT %>% 
	RunTSNE( dims = 1:20) %>% 
	FindNeighbors( dims = 1:20) %>% 
	FindClusters(resolution = resolution) %>%   identity()
# find all markers of clusters###
PT.markers <- FindAllMarkers(PT,  min.pct = 0.25, logfc.threshold = 0.25)

top20markers <- PT.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)  
#########################remove cluster and reclustering ############################

PT$previous<-PT$RNA_snn_res.2
idents= c(0,2,3,19,22)  ##remove  cells##
cells.use = rownames(PT@meta.data)[!Idents(PT) %in% idents]
PT<- subset(x = PT, cells =cells.use , invert = FALSE)  ##get count mat##

resolution<-2
PT <- PT %>% 
    RunTSNE(dims = 1:20) %>% 
    FindNeighbors(dims = 1:20) %>% 
    FindClusters(resolution = resolution) %>%     identity()
save(PT,file="PT.recluster.rds")
PT.markers <- FindAllMarkers(PT, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top20_markers<-PT.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)  ##top 1 of each cluster 
##rename each cluster based on marker genes#
new.cluster.ids <-c("Microglia","BAMs","Microglia","Microglia","Microglia","Microglia","Astrocytes","Astrocytes","Monocytes","Astrocytes","Oligodendrocytes","Microglia","DCs","Microglia","Granulocytes","Astrocytes","OPCs","NK&T cells","OPCs","Endothelial cells","Astrocytes","Endothelial cells","Fibroblast-like cells","Astrocytes","Oligodendrocytes","Endothelial cells","DCs","Mural cells","Mural cells","OPCs","Mural cells","Mural cells","Oligodendrocytes","Endothelial cells","Microglia","NK&T cells")
names(new.cluster.ids) <- levels(PT)
PT <- RenameIdents(PT, new.cluster.ids) #
save(PT,file="PT.rename.rds") ##
################Figure 3A ######################################################################################
cols<-c("#FED023","#CE5227","#5492CD","#E12F8B","#833C64","#23936A","#FF7F05","#CD9201","#33BEB6","#4B4E4F","#5B5FAA","#9D73B0","#999999")
x<-as.data.frame(PT@reductions$tsne@cell.embeddings)
x$orig.ident<-rownames(x) 
cell.type<-Idents(PT)
cell.type<-as.data.frame(cell.type)
cell.type$orig.ident<-rownames(cell.type)  
c<-merge(x,cell.type,by='orig.ident')
c$cell.type<- factor(c$cell.type,  
                             levels=c("Microglia","BAMs","Astrocytes","OPCs","Oligodendrocytes","Monocytes","DCs","Granulocytes","NK&T cells","Endothelial cells","Mural cells","Fibroblast-like cells"))

xlim<-c(min(c$tSNE_1)-3,max(c$tSNE_1)+2)
ylim<-c(min(c$tSNE_2)-2,max(c$tSNE_2)+2)
bwidth = .05*c(sum(abs(range(c$tSNE_1))),sum(abs(range(c$tSNE_2))))

stat_density <- ggplot(c,aes(tSNE_1,tSNE_2,colour=cell.type)) + xlim(xlim)+ ylim(ylim)+
	stat_density_2d(size = 0.25, colour = "#B5BABE",h=bwidth ,contour_var = "ndensity") 
	
gg =stat_density+geom_point(data = c, aes(x = tSNE_1, y = tSNE_2),size=0.4,stroke=0.8) +  ##node size, stroke size ##
	scale_color_manual(values=cols)+
	       
	guides(fill=guide_legend(title=NULL),legend.position=c(max(c$tSNE_1)+2,mean(ylim)),legend.justification = c(0, 1))+
	guides(colour = guide_legend(override.aes = list(size=2)))+ 
	theme(legend.text=element_text(size=8),legend.title=element_blank())+
	 theme_classic()	

ggsave(gg,file="Fig.3A-tSNE1.pdf",width=7.5,height=6)
#################
PT_PT<- subset(PT, stim == "PT")  ##
x<-as.data.frame(PT_PT@reductions$tsne@cell.embeddings)
x$orig.ident<-rownames(x) 
cell.type<-Idents(PT_PT)
cell.type<-as.data.frame(cell.type)

cell.type$orig.ident<-rownames(cell.type)  

c<-merge(x,cell.type,by='orig.ident')
c$cell.type<- factor(c$cell.type,  
                             levels=c("Microglia","BAMs","Astrocytes","OPCs","Oligodendrocytes","Monocytes","DCs","Granulocytes","NK&T cells","Endothelial cells","Mural cells","Fibroblast-like cells"))

bwidth = .05*c(sum(abs(range(c$tSNE_1))),sum(abs(range(c$tSNE_2))))
##
gg =stat_density+geom_point(data = c, aes(x = tSNE_1, y = tSNE_2),size=0.4,stroke=0.8) +  ##node size, stroke size ##
	scale_color_manual(values=cols)+
	       
	guides(fill=guide_legend(title=NULL),legend.position=c(max(c$tSNE_1)+2,mean(ylim)),legend.justification = c(0, 1))+
	guides(colour = guide_legend(override.aes = list(size=2)))+ 
	theme(legend.text=element_text(size=8),legend.title=element_blank())+
	 theme_classic()	
	
ggsave(gg,file="Fig.3A.PT_tSNE.pdf",width=7.5,height=6)
###############
PT_Sham<- subset(PT, stim == "Sham")  ##
x<-as.data.frame(PT_Sham@reductions$tsne@cell.embeddings)
x$orig.ident<-rownames(x) 
cell.type<-Idents(PT_Sham)
cell.type<-as.data.frame(cell.type)

cell.type$orig.ident<-rownames(cell.type)  

c<-merge(x,cell.type,by='orig.ident')
c$cell.type<- factor(c$cell.type,     #######sham  remove                                                           Granulocytes######
                             levels=c("Microglia","BAMs","Astrocytes","OPCs","Oligodendrocytes","Monocytes","DCs","Granulocytes","NK&T cells","Endothelial cells","Mural cells","Fibroblast-like cells"))
cols<-                              c("#FED023",         "#5492CD", "#E12F8B", "#833C64",       "#23936A","#FF7F05",             "#33BEB6",   "#4B4E4F",          "#5B5FAA",    "#9D73B0")
bwidth = .05*c(sum(abs(range(c$tSNE_1))),sum(abs(range(c$tSNE_2))))
##
gg =stat_density+geom_point(data = c, aes(x = tSNE_1, y = tSNE_2),size=0.4,stroke=0.8) +  ##node size, stroke size ##
	scale_color_manual(values=cols)+
	       
	guides(fill=guide_legend(title=NULL),legend.position=c(max(c$tSNE_1)+2,mean(ylim)),legend.justification = c(0, 1))+
	guides(colour = guide_legend(override.aes = list(size=2)))+ 
	theme(legend.text=element_text(size=8),legend.title=element_blank())+
	 theme_classic()
ggsave(gg,file="Figure3A.Sham_tSNE.pdf",width=7.5,height=6)
############
cell_order<-c("Microglia","BAMs","Astrocytes","OPCs","Oligodendrocytes","Monocytes","DCs","Granulocytes","NK&T cells","Endothelial cells","Mural cells","Fibroblast-like cells")

PT@active.ident <- factor(PT@active.ident, 
                            levels=cell_order)
Num<-table(Idents(PT),PT$stim)
Group <- c(rep("PT" , length(cell_order)) , rep("Sham" ,  length(cell_order))  )
cell.type <- c(cell_order,cell_order)
percent <- c(Num[,2],Num[,1])
data <- data.frame(Group,cell.type ,percent)
data$cell.type <- factor(data$cell.type,levels =cell_order)
# Stacked + percent
cols<-c("#FED023","#CE5227","#5492CD","#E12F8B","#833C64","#23936A","#FF7F05","#CD9201","#33BEB6","#4B4E4F","#5B5FAA","#9D73B0","#999999")
gg<-ggplot(data, aes(fill=cell.type, y=percent, x=Group)) + 
    scale_fill_manual(values=cols)+
	geom_bar( position="fill",stat="identity",width = 0.7)+
	theme_classic()
		
ggsave(gg,file="Figure3A.proportion.pdf")	
#########Figure 3B Dotplot######################################################################################
PT@active.ident <- factor(PT@active.ident,  
                             levels=c("Microglia","BAMs","Astrocytes","OPCs","Oligodendrocytes","Monocytes","DCs","Granulocytes","NK&T cells","Endothelial cells","Mural cells","Fibroblast-like cells"))
dotplot_markergene<-c("Tmem119","Siglech","Ms4a7","Pf4","Mfge8","Aqp4","Pdgfra","Nnat","Ptgds","Plp1","Ccr2","Ly6c2","H2-Aa","Cd74","S100a9","S100a8","Cd3e","Cd3g","Ly6c1","Igfbp7","Tagln","Acta2","Col1a1","Col3a1")
pdf("Fig.3B.Dotplot.pdf",width=11,height=5) ######dotplot of marker gene
DotPlot(object = PT, features = dotplot_markergene,cols =c("#4F7936","#FECA0A"))+theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.4,
 colour=c("#FED023","#FED023","#CE5227","#CE5227","#5492CD","#5492CD","#E12F8B","#E12F8B","#833C64","#833C64","#23936A","#23936A","#FF7F05","#FF7F05","#CD9201","#CD9201","#33BEB6","#33BEB6","#4B4E4F","#4B4E4F","#5B5FAA","#5B5FAA","#9D73B0","#9D73B0","#999999","#999999")),   
  axis.text.y = element_text(  colour=c("#FED023","#CE5227","#5492CD","#E12F8B","#833C64","#23936A","#FF7F05","#CD9201","#33BEB6","#4B4E4F","#5B5FAA","#9D73B0","#999999")))
dev.off()

####### draw the barplot###################################
Count_No_Feature<-function(cell_type){
	new_data <- subset(x = PT, idents = c(cell_type), invert = FALSE)
	return(c(mean(new_data$nFeature_RNA),##gene##
			mean(new_data$nCount_RNA))  )##UMI ##
}
cell_list<-as.character(unique(Idents(PT)))
cell_type_nFeature_result<-c()
for(i in 1:length(cell_list)){
	cell_type_nFeature_result<-rbind(cell_type_nFeature_result,c(Count_No_Feature(cell_list[i])))
}
colnames(cell_type_nFeature_result)<-c("mean.detected.gene","mean.UMI"); rownames(cell_type_nFeature_result)<-cell_list;
barplot_data <- data.frame(cell.type=as.character(rownames(cell_type_nFeature_result)),mean.UMI=cell_type_nFeature_result[,2])
barplot_data$cell.type <- factor(barplot_data$cell.type,levels = rev(c("Microglia","BAMs","Astrocytes","OPCs","Oligodendrocytes","Monocytes","Macrophages","DCs","Granulocytes","NK&T cells","Endothelial cells","Mural cells","Fibroblast-like cells")))
gg<-ggplot(barplot_data, aes( x=cell.type,y=mean.UMI)) + 
    geom_bar(stat = "identity",width = 0.7,fill=c("#999999"))+
	theme_classic()
ggsave(gg,file="Fig.3B.barplot.pdf",width=12,height=5)
markergene_list<-c("Siglech","Ms4a7","Aqp4","Pdgfra","Ptgds","Ly6c2","H2-Aa","S100a9","Cd3e","Ly6c1","Acta2","Col1a1")
pdf("Figure3C.FeaturePlot.pdf")
FeaturePlot(PT, reduction="tsne",features = markergene_list,pt.size = .001, cols = c("#a1a3a6", "#840228") )
dev.off()  
###################################################################
pltdata<-PT$nFeature_RNA
sd<-sd(pltdata)
pdf("Fig.S4B.nFeature.sd.pdf")
hist(pltdata, breaks = 60, main = "Distribution of unique genes per cell", col = "#E5FBB3", border="#74c69d",cex.main=2,
        xlab = expression("Number of unique genes per cell"), ylab = "Number of cells", cex.lab = 1.7, cex.axis = 1.7)
abline(v=mean(pltdata), col="black",lty=5,lwd=3)	
abline(v=c(mean(pltdata)-sd,mean(pltdata)+sd),lty=5, col="red",lwd=3)	
legend("topright", c(paste("mean =",round(mean(pltdata)),sep=" "), paste("SD =",round(sd),sep=" ")), col=c("black", "red"), bty="n",lwd=3)
dev.off()
########################
pltdata<-PT$nCount_RNA
sd<-sd(pltdata)
pdf("Fig.S4C.nCount_RNA.sd.pdf")

hist(pltdata, breaks = 60, main = "Distribution of unique transcripts per cell", col = "#f0ead2", border="#dde5b6",cex.main=2,
        xlab = expression("Number of unique transcripts per cell"), ylab = "Number of cells", cex.lab = 1.7, cex.axis = 1.7)
abline(v=mean(pltdata), col="black",lty=5,lwd=3)
abline(v=c(mean(pltdata)-sd,mean(pltdata)+sd),lty=5, col="red",lwd=3)
legend("topright", c(paste("mean =",round(mean(pltdata)),sep=" "), paste("SD =",round(sd),sep=" ")), col=c("black", "red"), bty="n",lwd=3)
dev.off()
###
PT@active.ident <- factor(PT@active.ident, 
                            levels=c("Microglia","BAMs","Astrocytes","OPCs","Oligodendrocytes","Monocytes","DCs","Granulocytes","NK&T cells","Endothelial cells","Mural cells","Fibroblast-like cells"))
Idents(PT)<-PT$seurat_clusters
PT@active.ident <- factor(PT@active.ident, levels=c("0","2","3","4","5","11","13","34","1","6","7","9","15","20","23","16","18","29","10","24","32","8","12","26","14","17","35","19","21","25","33","27","28","30","31","22"))
marker_list<-c("Fcrls","C1qa","Cx3cr1","Hexb","Siglech","Mt3","Clu","Plpp3","Slc1a2","Atp1a2","Ms4a7","Pf4","Apoc2","Lyz2","Spp1","Cd74","H2-Aa","H2-Ab1","H2-Eb1","Ifitm1","Plp1","Mal","Ptgds","Cldn11","Mag","","Pdgfra","Nnat","Gng3","Rtn1","Maged1","Ly6c1","Igfbp7","Itm2a","Cldn5","Ctla2a","S100a9","S100a8","Il1r2","Igfbp6","G0s2","Cd3e","Cd3d","Gzma","Nkg7","Ms4a4b","Col1a1","Col3a1","Col1a2","Dcn","Bgn","Myh11","Acta2","Tagln","Myl9")
pdf("Fig.S4D.heatmap.pdf",width=12,height=6) 
DoHeatmap(PT, 
   angle = 90,size=2, features = marker_list
   # ,slot="counts"
   )+scale_fill_gradientn(colors = c("#50114F", "#000000", "#FDFD0F"))+ theme(axis.text.y = element_text(size = 6))
 dev.off()   ###

#############################
x<-as.data.frame(PT@reductions$tsne@cell.embeddings)
x$orig.ident<-rownames(x)  ##dimension reduction result#######

Data<-PT@meta.data
Data<-as.data.frame(Data)
Data$orig.ident<-rownames(Data)  ##
c<-merge(x,Data,by='orig.ident')
xlim<-c(min(c$tSNE_1)-5,max(c$tSNE_1)+5)
ylim<-c(min(c$tSNE_2)-5,max(c$tSNE_2)+5)
bwidth = .05*c(sum(abs(range(c$tSNE_1))),sum(abs(range(c$tSNE_2))))
gg = ggplot(c,aes(tSNE_1,tSNE_2,colour=stim)) + xlim(xlim)+ ylim(ylim)+
	stat_density_2d(size = 0.25, colour = "#838C93",h=bwidth ,contour_var = "ndensity") +
	geom_point(size=0.6,stroke=0.6) +             
	guides(colour = guide_legend(override.aes = list(size=3)))  +  #legend size
	scale_color_manual(values=c("#4F7936","#FECA0A"))+
	guides(fill=guide_legend(title=NULL),legend.position=c(max(c$tSNE_1)+5,mean(ylim)),legend.justification = c(0, 1))+
    theme_classic()	
ggsave(gg,file="Fig.S4E.stim.pdf",width=8,height=7)

##################################################################################################################################
##########cell cell communication analysis##############################################################################
##update the cellchatDB###########

library(CellChat)
options(stringsAsFactors = FALSE)
CellChatDB <- CellChatDB.mouse # set CellChatDB <- CellChatDB.human if working on the human dataset
interaction_input <- CellChatDB$interaction
complex_input <- CellChatDB$complex
cofactor_input <- CellChatDB$cofactor
geneInfo <- CellChatDB$geneInfo
write.csv(interaction_input, file = "interaction_input_CellChatDB.csv")
write.csv(complex_input, file = "complex_input_CellChatDB.csv")
write.csv(cofactor_input, file = "cofactor_input_CellChatDB.csv")
write.csv(geneInfo, file = "geneInfo_input_CellChatDB.csv")

###read the updated data###
options(stringsAsFactors = FALSE)
interaction_input <- read.csv(file = 'interaction_input_CellChatDB.csv', row.names = 1)
complex_input <- read.csv(file = 'complex_input_CellChatDB.csv', row.names = 1)
cofactor_input <- read.csv(file = 'cofactor_input_CellChatDB.csv', row.names = 1)
geneInfo <- read.csv(file = 'geneInfo_input_CellChatDB.csv', row.names = 1)
CellChatDB <- list()
CellChatDB$interaction <- interaction_input
CellChatDB$complex <- complex_input
CellChatDB$cofactor <- cofactor_input
CellChatDB$geneInfo <- geneInfo
####
CellChatDB.mouse <- CellChatDB
usethis::use_data(CellChatDB.mouse, overwrite = TRUE)


library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
options(stringsAsFactors = FALSE)



PT@active.ident <- factor(PT@active.ident,  
                             # levels=c("Microglia","BAMs","Astrocytes","OPCs","Oligodendrocytes","Monocytes","Macrophages","DCs","Granulocytes","NK&T cells","Endothelial cells","Mural cells","Fibroblast-like cells"))
                             levels=c("Microglia","BAMs","Astrocytes","OPCs","Oligodendrocytes","Monocytes","DCs","Granulocytes","NK&T cells","Endothelial cells","Mural cells","Fibroblast-like cells"))



PT_PT<-subset(x = PT, subset= stim =="PT" , invert = FALSE)  ##get count mat##
data_input_PT <- GetAssayData(PT_PT, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(PT_PT)
meta_PT <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
##create CellChat object#
cellchat_PT <- createCellChat(object = data_input_PT, meta = meta_PT, group.by = "labels")
######add meta data##
cellchat_PT <- addMeta(cellchat_PT, meta = meta_PT, meta.name = "labels")
cellchat_PT <- setIdent(cellchat_PT, ident.use = "labels") # set "labels" as default cell identity

CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling  "Secreted Signaling" "ECM-Receptor"       "Cell-Cell Contact" 
cellchat_PT@DB <- CellChatDB.use
cellchat_PT <- subsetData(cellchat_PT) # subset the expression data of signaling genes for saving computation cost
######
future::plan("multiprocess", workers = 4) # do parallel

cellchat_PT <- identifyOverExpressedGenes(cellchat_PT)
cellchat_PT <- identifyOverExpressedInteractions(cellchat_PT)
cellchat_PT <- projectData(cellchat_PT, PPI.mouse)
###2 Inference of cell-cell communication network ###########################################################################################################
cellchat_PT <- computeCommunProb(cellchat_PT, raw.use = TRUE)
# Filter
cellchat_PT <- filterCommunication(cellchat_PT, min.cells = 10)
#######
cellchat_PT <- computeCommunProbPathway(cellchat_PT)
cellchat_PT <- netAnalysis_computeCentrality(cellchat_PT, slot.name = "netP")

cellchat_PT <- aggregateNet(cellchat_PT)
###############################################################################################
PT_Sham<-subset(x = PT, subset= stim =="Sham" , invert = FALSE)  ##get count mat##
data_input_Sham <- GetAssayData(PT_Sham, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(PT_Sham)
meta_Sham <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
##create CellChat object#
cellchat_Sham <- createCellChat(object = data_input_Sham, meta = meta_Sham, group.by = "labels")
######add meta data##
cellchat_Sham <- addMeta(cellchat_Sham, meta = meta_Sham, meta.name = "labels")
cellchat_Sham <- setIdent(cellchat_Sham, ident.use = "labels") # set "labels" as default cell identity

###############
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling  "Secreted Signaling" "ECM-Receptor"       "Cell-Cell Contact" 
cellchat_Sham@DB <- CellChatDB.use
cellchat_Sham <- subsetData(cellchat_Sham) # subset the expression data of signaling genes for saving computation cost
######
future::plan("multiprocess", workers = 4) # do parallel

cellchat_Sham <- identifyOverExpressedGenes(cellchat_Sham)
cellchat_Sham <- identifyOverExpressedInteractions(cellchat_Sham)
cellchat_Sham <- projectData(cellchat_Sham, PPI.mouse)
###2 Inference of cell-cell communication network ###########################################################################################################
cellchat_Sham <- computeCommunProb(cellchat_Sham, raw.use = TRUE)
# Filter
cellchat_Sham <- filterCommunication(cellchat_Sham, min.cells = 10)
#######
cellchat_Sham <- computeCommunProbPathway(cellchat_Sham)
cellchat_Sham <- netAnalysis_computeCentrality(cellchat_Sham, slot.name = "netP")
cellchat_Sham <- aggregateNet(cellchat_Sham)

save(cellchat_PT,file="cellchat_PT.rds")
save(cellchat_Sham,file="cellchat_Sham.rds")


######################### Lift up CellChat object and merge together
group.new = levels(cellchat_PT@idents)  ###
cellchat_Sham <- liftCellChat(cellchat_Sham, group.new)
#> The CellChat object will be lifted up using the cell labels FIB-A, FIB-B, FIB-P, DC, Pericyte, MYL, Immune, ENDO, Muscle, MELA, Basal-P, Basal, Spinious
#> Update slots object@net, object@netP, object@idents in a single dataset...
object.list <- list(PT = cellchat_PT, Sham = cellchat_Sham)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "structural")

cols<-c("#FED023","#CE5227","#5492CD","#E12F8B","#833C64","#23936A","#FF7F05","#CD9201","#33BEB6","#4B4E4F","#5B5FAA","#9D73B0","#999999")
weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
par(mfrow = c(1,2), xpd=TRUE)
pdf("Fig.3E.Interaction weights.pdf")
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight,color.use = cols, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interaction weights/strength - ", names(object.list)[i]))
}
dev.off()

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
pdf("Fig.3F.Number of interactions.pdf")
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count,color.use = cols, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()
library(ComplexHeatmap)
i=1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
cols<-c("#FED023","#CE5227","#5492CD","#E12F8B","#833C64","#23936A","#FF7F05","#CD9201","#33BEB6","#4B4E4F","#5B5FAA","#9D73B0")
names(cols)<-c("Microglia","BAMs","Astrocytes","OPCs","Oligodendrocytes","Monocytes","DCs","Granulocytes","NK&T cells","Endothelial cells","Mural cells","Fibroblast-like cells")
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]],color.use=cols, pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i],font.size = 6, width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]],color.use=cols, pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1],font.size = 6, width = 5, height = 6)
pdf("Fig.3G.outgoing.pdf")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()
#####################################
gg1 <- rankNet(cellchat, mode = "comparison", font.size = 14, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", font.size = 14, stacked = F, do.stat = TRUE)
pdf("Fig.S4F.comparison.pdf",width=18,height=12)
	print(gg1 + gg2);
dev.off()

#########Differetial expression genes analsysis######################################################################################################
interaction_input <- read.csv(file = 'interaction_input_CellChatDB.csv', row.names = 1)
complex_input <- read.csv(file = 'complex_input_CellChatDB.csv', row.names = 1)
cofactor_input <- read.csv(file = 'cofactor_input_CellChatDB.csv', row.names = 1)
####  get the secreted signaling ligand gene ###
interaction_input<-interaction_input[which(interaction_input$annotation=="Secreted Signaling"),]
protein_complex<-intersect(rownames(complex_input),interaction_input$ligand)  ##ligand
gene_complex<-complex_input[which(rownames(complex_input)==protein_complex),] ##ligand_complex 
ligand_gene_complex<-c(gene_complex[,1],gene_complex[,2]) 
ligand_gene<-setdiff(interaction_input$ligand,rownames(complex_input))
ligand_gene<-unique(c(ligand_gene,ligand_gene_complex))

#######R3.6 ####https://www.jianshu.com/p/41375176906e
retained_genes<-c(ligand_gene)
features <- rownames(PT)[rownames(PT) %in% retained_genes]
PT <- subset(PT, features = features)

Idents(PT)<-PT$stim
DEG_method<-"wilcox"
PT <- FindVariableFeatures(PT, selection.method = "vst", nfeatures = 2000)
ALL_DEGs <- FindMarkers(PT, ident.1 = "PT", ident.2 = NULL,test.use = DEG_method,logfc.threshold =0,    min.pct = 0)
sig_DEGs <- FindMarkers(PT, ident.1 = "PT", ident.2 = NULL,test.use = DEG_method,logfc.threshold = 0.25,min.pct = 0.1)

ALL_DEGs1<-subset(ALL_DEGs, rownames(ALL_DEGs) %in% ligand_gene)
sig_DEGs1<-subset(sig_DEGs, rownames(sig_DEGs) %in% ligand_gene)

sig_DEGs1<- subset(x = ALL_DEGs1, p_val_adj<0.05 &  (avg_logFC > 0.1 | avg_logFC < 0.1))
ALL_DEGs1$logFC =  ifelse(ALL_DEGs1$avg_logFC> 0 ,ALL_DEGs1$avg_logFC,abs(ALL_DEGs1$avg_logFC))
sig_DEGs1$threshold = as.factor( ifelse(sig_DEGs1$avg_logFC> 0 ,'#335A28','#AC6DAD'))
up_DEGs<- subset(x = sig_DEGs1, p_val_adj<0.05 &  avg_logFC > 0.25)
down_DEGs<- subset(x = sig_DEGs1, p_val_adj<0.05 &  avg_logFC < (-0.25))
selected_genes<-c("Cxcl2","Spp1","Lgals3","Lgals9","Il1b","Lgals1","Mif","Cxcl16","Anxa1","Pf4","Il1a","Il18","Ccl7","Ccl2","C4b","Osm","Nampt","Enho","Tnfsf9","Ccl5","Gdf15")
sig_DEGs1<-subset(sig_DEGs1, rownames(sig_DEGs1) %in% selected_genes)

####################################################
gg2<-ggplot(data = ALL_DEGs1, aes(x = avg_logFC, y = -log10(p_val_adj))) +
  # geom_point(alpha=0.5,  aes(size=logFC,colour = avg_logFC)) +
  geom_point( alpha=0.8,aes(size=logFC,colour = avg_logFC)) +
  # scale_color_manual(values=c("blue", "grey","red"))+ 
  xlim(c(-3.5,3.5))+
  labs(x="log2(fold change)",y="-log10(p.adj)",title="PT vs Sham (secreted ligand genes)")+
  theme( axis.title.x = element_text(size = 36),  axis.text.x = element_text(size = 36),
		axis.text.y = element_text(size = 36),  axis.title.y = element_text(size = 36))+
  geom_text_repel(
    data = sig_DEGs1,
    aes(x = avg_logFC,  y = -log10(p_val_adj),label = rownames(sig_DEGs1)),
    size = 5.2,
	col= sig_DEGs1$threshold,box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"),
	segment.color = "black", show.legend = FALSE )+   ##https://ggplot2.tidyverse.org/reference/scale_gradient.html
	scale_colour_gradientn(colours=c("#141378","#30AF95","#F0EA07"))+
	scale_size_continuous(range = c(1, 6))+ #https://community.rstudio.com/t/ggplot2-is-it-possible-to-combine-color-fill-and-size-legends/17072
	theme_classic()

ggsave(gg2,device="pdf",file="Fig.3H.Secreted ligand gene.pdf")












