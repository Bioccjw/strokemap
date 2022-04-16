#######subtypes re-clustering#############################################
######Microglia####################################################################################################################################
library(dplyr)
library(Seurat)
require(Matrix)
require(magrittr)
library(ggplot2)
load("PT.rename.rds")
subset_data<-subset(x = PT, idents =c("Microglia") , invert = FALSE)  ##
subset_data[["percent.mt"]] <- PercentageFeatureSet(subset_data, pattern = "^mt-")
#######remove mt genes
mt.genes <- rownames(subset_data)[grep("^mt-",rownames(subset_data))]
remove_genes<-c(mt.genes)
features <- rownames(subset_data)[!rownames(subset_data) %in% remove_genes]
subset_data <- subset(subset_data, features = features)

subset_data <- FindVariableFeatures(subset_data, selection.method = "vst", nfeatures = 2000)
subset_data <- RunPCA(subset_data, features = VariableFeatures(object = subset_data))
subset_data <- RunTSNE(subset_data, dims = 1:10)

subset_data <- FindNeighbors(subset_data, dims = 1:10)
subset_data <- FindClusters(subset_data, resolution = 0.25)

new.cluster.ids <-c("MG1","MG2","MG3","MG4","MG5") 
names(new.cluster.ids) <- levels(subset_data)
subset_data <- RenameIdents(subset_data, new.cluster.ids)  #rename cell clusters 
subset_data@active.ident <- factor(subset_data@active.ident, 
                            levels=c("MG1","MG2","MG3","MG4","MG5") )
x<-as.data.frame(subset_data@reductions$tsne@cell.embeddings)
x$orig.ident<-rownames(x)  
head(x)
######
cell.type<-Idents(subset_data)
cell.type<-as.data.frame(cell.type)
head(cell.type)
cell.type$orig.ident<-rownames(cell.type)  ##
head(cell.type)
c<-merge(x,cell.type,by='orig.ident')
xlim<-c(min(c$tSNE_1)-1,max(c$tSNE_1)+1.8)
ylim<-c(min(c$tSNE_2)-5,max(c$tSNE_2)+1.2)
bwidth = .05*c(sum(abs(range(c$tSNE_1))),sum(abs(range(c$tSNE_2))))
##
gg = ggplot(c,aes(tSNE_1,tSNE_2,colour=cell.type)) + xlim(xlim)+ ylim(ylim)+
	stat_density_2d(size = 0.25, colour = "#B5BABE",h=bwidth ,contour_var = "ndensity") +
	geom_point(size=0.6,stroke=0.6) +  ##
	guides(colour = guide_legend(override.aes = list(size=2)))  +  #legend size
	scale_color_manual(values=c("#5492CD","#23936A","#CE5227","#CD9201","#833C64","#E12F8B"))+
	guides(fill=guide_legend(title=NULL),legend.position=c(max(c$tSNE_1)+5,mean(ylim)),legend.justification = c(0, 1))+
    theme_classic()	

ggsave(gg,file="Fig.4B.MG.pdf",width=6,height=5)
table(Idents(subset_data),subset_data$stim)  ##Fig.4C##

x<-as.data.frame(subset_data@reductions$tsne@cell.embeddings)
x$orig.ident<-rownames(x)  ##
head(x)
######
Data<-subset_data@meta.data
Data<-as.data.frame(Data)
Data$orig.ident<-rownames(Data)  ##
c<-merge(x,Data,by='orig.ident')
xlim<-c(min(c$tSNE_1)-5,max(c$tSNE_1)+5)
ylim<-c(min(c$tSNE_2)-5,max(c$tSNE_2)+5)
bwidth = .05*c(sum(abs(range(c$tSNE_1))),sum(abs(range(c$tSNE_2))))
gg = ggplot(c,aes(tSNE_1,tSNE_2,colour=stim)) + xlim(xlim)+ ylim(ylim)+
	stat_density_2d(size = 0.25, colour = "#838C93",h=bwidth ,contour_var = "ndensity") +
	geom_point(size=0.6,stroke=0.6) +             ##
	guides(colour = guide_legend(override.aes = list(size=3)))  +  #legend size
	scale_color_manual(values=c("#4F7936","#FECA0A"))+
	guides(fill=guide_legend(title=NULL),legend.position=c(max(c$tSNE_1)+5,mean(ylim)),legend.justification = c(0, 1))+
    theme_classic()	
ggsave(gg,file="Fig.S5C.stim.pdf",width=6,height=5)
###################################
marker_list_manual<-c("Egr1","Junb","Dusp1","Fos","Flf4","P2ry12","Selplg","Malat1","Cx3cr1","Gpr34",
"Ifi27l2a","Tspo","Ifitm3","Cd52","Apoe","Lgals3","Vim","Spp1","Pgk1","Mif","Ube2c","Stmn1","Hist1h2ap","H2afz","Tubb5")
pdf("Fig.S5B.manual.MG.pdf",width=15,height=9) 
DoHeatmap(subset_data, disp.max = 3, disp.min = -1,slot = "data",
	angle = 0,size=5.52,marker_list_manual
)+scale_fill_gradientn(colors = c("#335A28", "#000000", "#FDFD0F"))+ theme(
axis.text.y = element_text(size = 16,colour=rep("#000000",length(marker_list_manual)),face="italic"))
dev.off() 
################################################################################
normalized_exp<-as.matrix(subset_data@assays$RNA@data)
input_markers<-c("Trem2","P2ry12","Tmem119","Klf4","Sertad1","Tagap","Tspo","Ifi27l2a","Ifitm3","Lgals3","Anxa2","Hilpda","Pclaf","Birc5","Tk1")

marker_profile<-c() ###get the barplot 
for(i in 1:length(input_markers)){
	marker_expr<-data.frame(cell.type=Idents(subset_data),expr=normalized_exp[which(rownames(normalized_exp)==input_markers[i]),])
	marker_expr<-marker_expr[order(marker_expr$cell.type),]
	num_cells=1:nrow(marker_expr)
	gene.symbol=rep(input_markers[i],nrow(marker_expr))
	marker_expr<-cbind(marker_expr,num_cells,gene.symbol)
	marker_profile<-rbind(marker_profile,marker_expr)
}

gg<-ggplot(marker_profile, aes(x = num_cells, y = expr, fill = cell.type)) +
	facet_grid(rows = vars(gene.symbol))+ylab("log(normalized counts)")+xlab("")+
	geom_col(position = "identity")+
	scale_fill_manual(values =  c("#5492CD","#23936A","#CE5227","#CD9201","#833C64"))+
	scale_y_continuous(breaks=seq(0, 6, 4))+
	scale_x_continuous(breaks=seq(0, 10000, 10000))+
	theme(axis.text.y = element_text(color = "grey20", size = 14), 
		 axis.title.y = element_text(color = "grey20", size = 18, face = "plain"))+	 
	theme(axis.text.x = element_blank())+ 
	theme(panel.border = element_blank(),panel.grid=element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black")) 
ggsave(gg,file="Fig.4D.MG_markers.pdf",width=9.5,height=8)
##############################################
pdf("Fig.4F.vlnPlot.pdf",width=15,height=9) 
VlnPlot(subset_data, features = c("Lgals1","Lgals3","Lgals9"),pt.size=.1,
cols=c("#5492CD","#23936A","#CE5227","#CD9201","#833C64"))
dev.off()
pdf("Fig.4G.DotPlot.pdf",width=9,height=15) 
DotPlot(object = subset_data, features = c("Lgals1","Lgals3","Lgals9"),cols =c("#008792", "#ed1941"))+theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.4))
dev.off()
#######BAMs############################################################################################

rm(list=ls())
library(dplyr)
library(Seurat)
subset_data<-subset(x = PT, idents =c("BAMs") , invert = FALSE)  ##

subset_data[["percent.mt"]] <- PercentageFeatureSet(subset_data, pattern = "^mt-")
#######remove mt gene
mt.genes <- rownames(subset_data)[grep("^mt-",rownames(subset_data))]
remove_genes<-c(mt.genes)
features <- rownames(subset_data)[!rownames(subset_data) %in% remove_genes]
subset_data <- subset(subset_data, features = features)

subset_data$cell_type<-Idents(subset_data)
subset_data <- FindVariableFeatures(subset_data, selection.method = "vst", nfeatures = 5000)
subset_data <- RunPCA(subset_data, features = VariableFeatures(object = subset_data))
subset_data <- RunTSNE(subset_data, dims = 1:11)

subset_data <- FindNeighbors(subset_data, dims = 1:11)
subset_data <- FindClusters(subset_data, resolution = 0.5)

new.cluster.ids <-c("BAM1","BAM2","BAM3","BAM4") 
names(new.cluster.ids) <- levels(subset_data)
subset_data <- RenameIdents(subset_data, new.cluster.ids)  #combine and rename cell clusters 
# DimPlot(subset_data,reduction="tsne",label=T)

subset_data@active.ident <- factor(subset_data@active.ident, 
                            levels=c("BAM1","BAM2","BAM3","BAM4")  )

x<-as.data.frame(subset_data@reductions$tsne@cell.embeddings)
x$orig.ident<-rownames(x)  
head(x)
######
cell.type<-Idents(subset_data)
cell.type<-as.data.frame(cell.type)
head(cell.type)
cell.type$orig.ident<-rownames(cell.type)  
head(cell.type)
c<-merge(x,cell.type,by='orig.ident')
xlim<-c(min(c$tSNE_1)-1,max(c$tSNE_1)+1.2)
ylim<-c(min(c$tSNE_2)-1,max(c$tSNE_2)+1.2)
bwidth = .05*c(sum(abs(range(c$tSNE_1))),sum(abs(range(c$tSNE_2))))
##
gg = ggplot(c,aes(tSNE_1,tSNE_2,colour=cell.type)) + xlim(xlim)+ ylim(ylim)+
	stat_density_2d(size = 0.25, colour = "#B5BABE",h=bwidth ,contour_var = "ndensity") +
	geom_point(size=0.8,stroke=0.8) +  ## 
	guides(colour = guide_legend(override.aes = list(size=2)))  +  #legend size
	scale_color_manual(values=c("#5492CD","#FED023","#CE5227","#E12F8B"))+
	guides(fill=guide_legend(title=NULL),legend.position=c(max(c$tSNE_1)+5,mean(ylim)),legend.justification = c(0, 1))+
    theme_classic()	
ggsave(gg,file="Fig.4I.tSNE.pdf",width=5,height=5)
####################################
cluster.averages <- AverageExpression(subset_data, 
                                      assays = "RNA", 
                                      return.seurat = TRUE,
                                      slot = "data")
cluster.averages <- ScaleData(cluster.averages)

subset_data.markers <- FindAllMarkers(subset_data,  min.pct = 0.25, logfc.threshold = 0.25)
topmarkers<-subset_data.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)  ##top 1 of each cluster 
mat <- GetAssayData(cluster.averages, slot = "scale.data")
mat <- as.matrix(mat[topmarkers$gene, ])
library("reshape2")
data1 <- melt(mat)
gene_order<-rev(unique(as.character(data1$Var1)))
data1$Var1 <- factor(data1$Var1, levels=gene_order)
colnames(data1)<-c("Var1","Var2","Expression")

gg<-ggplot(data1,aes(x = Var2, y = Var1, fill = Expression))+
	geom_tile()+
	xlab("")+ylab("")+
	theme_classic()+
	scale_fill_gradientn(colors = c("#50114F", "#000000", "#FDFD0F"),limits=c(-2, 2))+
	theme(text = element_text(size=20,),
        axis.text.x = element_text(angle=0,colour=rep("#000000",4)),
		axis.text.y = element_text(angle=0,face="italic",colour=rep("#000000",20))) 
	
ggsave(gg,file="Fig.4J.AverageExpression.pdf")
##########################################################
pdf("Fig.4L.vlnPlot.pdf",width=15,height=9) 
VlnPlot(subset_data, features = c("Lgals1","Lgals3","Lgals9"),pt.size=.1,
cols=c("#5492CD","#FED023","#CE5227","#E12F8B"))
dev.off()

##################################################################################################################################
#######Astrocytes######################################################################################
pdf("Fig.6A.Cd44.pdf",width=8,height=8.6)
FeaturePlot(PT, reduction="tsne",features = "Cd44",pt.size = 1.2, cols = c("#a1a3a6", "#840228"))
dev.off()
subset_data<-subset(x = PT, idents =c("Astrocytes") , invert = FALSE)  ##
subset_data[["percent.mt"]] <- PercentageFeatureSet(subset_data, pattern = "^mt-")
#######remove mt gene
mt.genes <- rownames(subset_data)[grep("^mt-",rownames(subset_data))]
remove_genes<-c(mt.genes)
features <- rownames(subset_data)[!rownames(subset_data) %in% remove_genes]
subset_data <- subset(subset_data, features = features)

subset_data <- FindVariableFeatures(subset_data, selection.method = "vst", nfeatures = 2000)
subset_data <- RunPCA(subset_data, features = VariableFeatures(object = subset_data))

subset_data <- FindNeighbors(subset_data, dims = 1:9)
subset_data <- FindClusters(subset_data, resolution = 0.4)
subset_data <- RunTSNE(subset_data, dims = 1:9)
###############################
new.cluster.ids <-c("AS1","AS2","AS4","AS5","AS3","AS6") 
names(new.cluster.ids) <- levels(subset_data)
subset_data <- RenameIdents(subset_data, new.cluster.ids)  #combine and rename cell clusters 
subset_data@active.ident <- factor(subset_data@active.ident, 
                            levels=c("AS1","AS2","AS3","AS4","AS5","AS6") )
x<-as.data.frame(subset_data@reductions$tsne@cell.embeddings)
x$orig.ident<-rownames(x)  
head(x)
######
cell.type<-Idents(subset_data)
cell.type<-as.data.frame(cell.type)
head(cell.type)
cell.type$orig.ident<-rownames(cell.type)  ##
head(cell.type)
c<-merge(x,cell.type,by='orig.ident')
####
xlim<-c(min(c$tSNE_1)-2,max(c$tSNE_1)+1)
ylim<-c(min(c$tSNE_2)-1,max(c$tSNE_2)+1.2)
bwidth = .05*c(sum(abs(range(c$tSNE_1))),sum(abs(range(c$tSNE_2))))
##
gg = ggplot(c,aes(tSNE_1,tSNE_2,colour=cell.type)) + xlim(xlim)+ ylim(ylim)+
	stat_density_2d(size = 0.25, colour = "#B5BABE",h=bwidth ,contour_var = "ndensity") +
	geom_point(size=0.6,stroke=0.6) +  
	guides(colour = guide_legend(override.aes = list(size=2)))  +  #legend size
	scale_color_manual(values=c("#5492CD","#833C64","#4B4E4F","#E12F8B","#CD9201","#FF7F05"))+
	guides(fill=guide_legend(title=NULL),legend.position=c(max(c$tSNE_1)+5,mean(ylim)),legend.justification = c(0, 1))+
    theme_classic()	
# gg

ggsave(gg,file="Fig.6B.Astrocytes.pdf",width=6,height=5)

x<-as.data.frame(subset_data@reductions$tsne@cell.embeddings)
x$orig.ident<-rownames(x)  ##降维结果##
head(x)
######
Data<-subset_data@meta.data
Data<-as.data.frame(Data)
Data$orig.ident<-rownames(Data)  ##细胞注释###
c<-merge(x,Data,by='orig.ident')
xlim<-c(min(c$tSNE_1)-5,max(c$tSNE_1)+5)
ylim<-c(min(c$tSNE_2)-5,max(c$tSNE_2)+5)
bwidth = .05*c(sum(abs(range(c$tSNE_1))),sum(abs(range(c$tSNE_2))))
gg = ggplot(c,aes(tSNE_1,tSNE_2,colour=stim)) + xlim(xlim)+ ylim(ylim)+
	stat_density_2d(size = 0.25, colour = "#838C93",h=bwidth ,contour_var = "ndensity") +
	geom_point(size=0.6,stroke=0.6) +             ##填充大小、线条粗细##
	guides(colour = guide_legend(override.aes = list(size=3)))  +  #legend size
	scale_color_manual(values=c("#4F7936","#FECA0A"))+
	guides(fill=guide_legend(title=NULL),legend.position=c(max(c$tSNE_1)+5,mean(ylim)),legend.justification = c(0, 1))+
    theme_classic()	

ggsave(gg,file="Fig.S6A.Astrocytes.stim.pdf",width=6,height=5)
##########################################################################
Num<-table(Idents(subset_data),subset_data$stim)
Group <- c(rep("PT" , length(cell_order)) , rep("Sham" ,  length(cell_order))  )
cell.type <- c(cell_order,cell_order)
percent <- c(Num[,2],Num[,1])
data <- data.frame(Group,cell.type ,percent)
data$cell.type <- factor(data$cell.type,levels =cell_order)
# Stacked + percent
gg<-ggplot(data, aes(fill=cell.type, y=percent, x=Group)) + 
    scale_fill_manual(values=c("#5492CD","#833C64","#4B4E4F","#E12F8B","#CD9201","#FF7F05"))+
	geom_bar(position = "fill",stat="identity",width = 0.5)+
	theme_classic()
ggsave(gg,file="Fig.6C.proportion.pdf",width=6,height=4)	
#################
cluster_order<-rev(c("AS1","AS2","AS3","AS4","AS5","AS6"))
subset_data@active.ident <- factor(subset_data@active.ident,  
                              levels=c(cluster_order))
pdf("Fig.6A.Cd44.dotplot.pdf",width=3.2,height=6) ######dotplot of marker gene
DotPlot(object = subset_data, features = "Cd44",cols =c("#4F7936","#FECA0A"))+theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.4))
dev.off()

pdf("Fig.6D.Cd44.pdf",width=5,height=6)
FeaturePlot(subset_data, reduction="tsne",features = "Cd44",pt.size = 1.2, cols = c("#008792", "#ed1941"))
dev.off()



subset_data@active.ident <- factor(subset_data@active.ident, 
                            levels=c("AS1","AS2","AS3","AS4","AS5","AS6") )

marker_list_manual<-c("Aldoc","Htra1","Cxcl14","S1pr1","Gpr37l1","Gdf10","Agt","Hopx","Hhatl","Id3","Gm26917","Syne1","Gria2","Dclk1","Nwd1","Lcn2","Serpina3n","Gfap","Timp1","A2m","Pclaf","Birc5","Top2a","Rrm2","Cdk1","Clcf1","Flnc","Tnfrsf12a","Fosl1","Tuba1c")
pdf("Doheatmap.manual.pdf",width=15,height=9) 
DoHeatmap(subset_data, disp.max = 2, disp.min = -1,
   # draw.lines=FALSE,
   angle = 0,size=5.52,marker_list_manual
   )+scale_fill_gradientn(colors = c("#335A28", "#000000", "#FDFD0F"))+ theme(
   axis.text.y = element_text(size = 16,colour=rep("#000000",length(marker_list_manual)),face="italic"))
dev.off()   ###
###############WGCNA analysis of Astrocytes subtypes######################################
#######construct the cell type binary matrix
cell.id<-as.data.frame(Idents(subset_data))
cell.type<-sort(unique(new.cluster.ids))
trait<-matrix(0,nrow=nrow(cell.id),ncol=length(cell.type))
for(i in 1:length(cell.type)){
	posi<-which(cell.id[,1]==cell.type[i])
	trait[posi,i]<-1
}
colnames(trait)<-cell.type
rownames(trait)<-rownames(cell.id)
setwd("/mnt/zhoushunheng/stroke/2seurat/result10_Astrocyte/5WGCNA")
subset_data <- subset_data %>% NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 5000) %>%
    ScaleData(verbose = FALSE) 
	
subset_data.matrix<-as.matrix(subset_data@assays$RNA@scale.data)
write.table(subset_data.matrix,file="ExpData.txt",sep="\t",quote=FALSE)
write.table(trait,file="Sam_info.txt",sep="\t",quote=FALSE)
##########WGCNA analysis######################################################

library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

samples=read.csv('Sam_info.txt',sep = '\t',row.names = 1)
expro=read.csv('ExpData.txt',sep = '\t',row.names = 1)

datExpr=as.data.frame(t(expro));
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


##soft threshold selection##
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")


plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

##One-step network construction and module detection##
# we choose 2 for network construction 
power = sft$powerEstimate
power
net = blockwiseModules(datExpr, power = power, maxBlockSize = 6000,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "AS-green-FPKM-TOM",
                       verbose = 3)
table(net$colors)

##result presentation##
# open a graphics window
#sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
moduleLabels = net$colors
# Plot the dendrogram and the module colors underneath

plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


MEs = net$MEs

MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
##########module correlation##
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)


moduleColors = labels2colors(net$colors)
table(moduleColors)

pdf("3.geneTree.pdf")
geneTree = net$dendrograms[[1]];
dev.off()
#########

load(net$TOMFiles[1], verbose=T)


TOM <- as.matrix(TOM)

dissTOM = 1-TOM
# Transform dissTOM with a power to make moderately strong 
# connections more visible in the heatmap
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function

##########
allmodules<-unique(labels2colors(net$colors))
for (i in 1:length(allmodules)){
modules =allmodules[i];
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,                               
                               nodeAttr = moduleColors[inModule]);
}


# table(rownames(trainDt) == rownames(datExpr))
moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic
MEs0 = moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, samples, use = "p")
colnames(MEsWW)
modlues=MEsWW
modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)


pdf("Fig.S6C.module_trait.correlation.pdf")
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(samples), yLabels = names(MEsWW), cex.lab = 0.5,  yColorWidth=0.01, 
               xColorWidth = 0.03,
               ySymbols = colnames(modlues), colorLabels = FALSE, 
			    textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1)
               , main = paste("Module-trait relationships"),
			 colors=c(colorRampPalette(colors = c('#2769A2', 'white'))( 100 ),colorRampPalette(colors = c('white','#7D0224'))( 100 ))
	
			   )
dev.off()


#######R3.6 
rm(list=ls())
library(ggplot2)
setwd("/mnt/zhoushunheng/stroke/2seurat/result10_Astrocyte/5WGCNA/")
enrichment_result=read.table(file="selected_term_updated.txt",header=TRUE,sep="\t",row.names=1,stringsAsFactors = F)
enrichment_result$module
head(enrichment_result)
library(Hmisc)
enrichment_result$Description<-capitalize(enrichment_result$Description)


enrichment_result$module<-factor(enrichment_result$module, levels = rev(unique(enrichment_result$module)))
library(dplyr)
enrichment_result<-enrichment_result %>% group_by(module)  %>% arrange(desc(pvalue),.by_group = TRUE)

enrichment_result<-as.data.frame(enrichment_result)

enrichment_result$pvalue<-as.numeric(enrichment_result$pvalue)
head(enrichment_result)

module_color<- rev(c("#22749D","#6A994E","#495057","#EF476F","#BC4749","#297E85","#000814","#dbb42c","#7f5539") ) ##从上往下设置颜色
module_color<- module_color[1:length(unique(enrichment_result$module))]##
# text_color<-rep(rev(module_color),table(enrichment_result$module))
text_color<-rep(module_color,table(enrichment_result$module))
p<-ggplot(enrichment_result, aes(x = Description,   y = -log10(pvalue),fill=module)) +
        geom_bar(stat = "identity",width=0.8) + 
		geom_text(aes(label=round(-log10(pvalue),2)), position=position_dodge(width=0.9), hjust=-0.1)+  ##text size 
		coord_flip()+  ##翻转############
		scale_y_continuous(limits = c(0, max(-log10(enrichment_result$pvalue))+10))+       ##y axis limit
		theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ 
		scale_fill_manual(values = module_color )+guides(fill = guide_legend(reverse=TRUE))+            ##reverse the legend 
		theme(text = element_text(size=15),axis.text.y = element_text(angle = 0,hjust=1,vjust=0.5, color=text_color),legend.position = c(0.8,0.382))+             ##text color
		ylab('-log10(p value)')+xlab('GO terms')+
		scale_x_discrete(limits=enrichment_result$Description)   ## fix the x axis
ggsave(p,file="Fig.S6D.Selected_terms.pdf",width=15,height=8)


















