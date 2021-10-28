#ref https://satijalab.org/seurat/v3.1/immune_alignment.html
#youtube:https://www.youtube.com/watch?v=PSuiBaUqh0Y&t=32s
library(Seurat)
library(SeuratData)
library(cowplot)
library(patchwork)
setwd("/gpfs/data/Chatv01_Lyny/NGS-class/FinalProject_May12/Seurat_tutorial/filtered_gene_bc_matrices")
                                  


                                    ##------------scRNA-seq FINAL PROJECT--------------
                                        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#1.................-------------reading data and cleaning------------------------------------
#--Exp1-control
exp1_control.data <- Read10X(data.dir ="./exp1_control/")
#-----Initialize the Seurat object with the raw (non-normalized data).
exp1_control.original <- CreateSeuratObject(counts = exp1_control.data, project = "exp1_control", min.cells = 3, min.features = 200)
exp1_control=exp1_control.original
#-----removing unwanted cells
exp1_control[["percent.mt"]] <- PercentageFeatureSet(exp1_control, pattern = "^Mt-")
# Visualize QC metrics as a violin plot
VlnPlot(exp1_control, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#identify outliers
summary(exp1_control$nFeature_RNA)
(upper.exp1_control=2568+1.5*IQR(exp1_control$nFeature_RNA))
lower.exp1_control=330-1.5*IQR(exp1_control$nFeature_RNA) #no lower bound
exp1_control <- subset(exp1_control, subset =  nFeature_RNA < upper.exp1_control & percent.mt < 10)

#--exp4-control
exp4_control.data <- Read10X(data.dir ="./exp4_control/")
#-----Initialize the Seurat object with the raw (non-normalized data).
exp4_control.original <- CreateSeuratObject(counts = exp4_control.data, project = "exp4_control", min.cells = 3, min.features = 200)
exp4_control=exp4_control.original
#removing unwanted cells
exp4_control[["percent.mt"]] <- PercentageFeatureSet(exp4_control, pattern = "^Mt-")
# Visualize QC metrics as a violin plot
VlnPlot(exp4_control, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#identify outliers
summary(exp4_control$nFeature_RNA)
(upper.exp4_control=530+1.5*IQR(exp4_control$nFeature_RNA))
(lower.exp4_control=251-1.5*IQR(exp4_control$nFeature_RNA)) #no lower bound
exp4_control <- subset(exp4_control, subset =  nFeature_RNA < upper.exp4_control & percent.mt < 10)

#--LyDMSO
Ly.data <- Read10X(data.dir ="./Ly/")
#-----Initialize the Seurat object with the raw (non-normalized data).
Ly.original <- CreateSeuratObject(counts = Ly.data, project = "Ly", min.cells = 3, min.features = 200)
Ly=Ly.original
#-----removing unwanted cells
Ly[["percent.mt"]] <- PercentageFeatureSet(Ly, pattern = "^Mt-")
# Visualize QC metrics as a violin plot
VlnPlot(Ly, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#identify outliers
summary(Ly$nFeature_RNA)
(upper.Ly=3087+1.5*IQR(Ly$nFeature_RNA))
lower.Ly=245-1.5*IQR(Ly$nFeature_RNA) #no lower bound
Ly <- subset(Ly, subset =  nFeature_RNA < upper.Ly & percent.mt < 10)

#--mirin
mirin.data <- Read10X(data.dir ="./mirin-exp4/")
#-----Initialize the Seurat object with the raw (non-normalized data).
mirin.original <- CreateSeuratObject(counts = mirin.data, project = "mirin", min.cells = 3, min.features = 200)
mirin=mirin.original
#-----removing unwanted cells
mirin[["percent.mt"]] <- PercentageFeatureSet(mirin, pattern = "^Mt-")
# Visualize QC metrics as a violin plot
VlnPlot(mirin, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#identify outliers: That the interquartile range can be used to identify outliers in data regardless of the distribution.
summary(mirin$nFeature_RNA)
upper.mirin=1159+1.5*IQR(mirin$nFeature_RNA)
lower.mirin=246-1.5*IQR(mirin$nFeature_RNA) #no lower bound
mirin <- subset(mirin, subset =  nFeature_RNA>lower.mirin &nFeature_RNA <upper.mirin & percent.mt < 10)


#2----------------------------------------------------------Integration----------------------------------
exp1_control@meta.data[,"Dataset"]<- "exp1-control"
exp4_control@meta.data[,"Dataset"]<- "exp4-control"
Ly@meta.data[,"Dataset"]<- "Ly"
mirin@meta.data[,"Dataset"]<- "mirin"

merged=merge(exp1_control,y=c(exp4_control,Ly,mirin),add.cell.ids=c("exp1-control","exp4-control","Ly","mirin"),project="Dataset")
merged.list<-SplitObject(merged,split.by = "Dataset")  

#Normalize
merged.list<- lapply(X = merged.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})  

##integration
anchors <- FindIntegrationAnchors(object.list = merged.list, dims = 1:30)  
combined<- IntegrateData(anchorset=anchors,dims=1:30)
DefaultAssay(combined)<- "integrated"

#3---------------------------------------Run the standard workflow for visualization and clustering-------------
set.seed(123456)
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 50, verbose = FALSE)
#Clustering
ElbowPlot(combined)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:15)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:15)

#what resolutions to choose?
#ref: https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html
library(clustree)
combined <- FindClusters(combined, resolution = c(0,0.02,0.04,0.06,0.08,0.1))
clustree(combined, prefix = "integrated_snn_res.",node_colour = "sc3_stability",layout = "sugiyama")
clustree(combined, prefix = "integrated_snn_res.") +guides(edge_colour = FALSE, edge_alpha = FALSE) +
  theme(legend.position = "bottom")

##---resolution FINAL= 0.04
combined <- FindClusters(combined, resolution = 0.04) #re-run the code and pick resolution =0.04 for stable clustering
# Visualization
DimPlot(combined, reduction = "umap", label = TRUE)
DimPlot(combined, reduction = "umap", split.by = "Dataset",label = TRUE)

#------------------Identify conserved cell type markers to ANNOTATE CLUSTERS
DefaultAssay(combined) <- "RNA"
for (i in 0:4){
  conserved_markers=FindConservedMarkers(combined, ident.1 = i, grouping.var = "Dataset", verbose = FALSE)
  conserved_markers$gene=rownames(conserved_markers)
  write.table(conserved_markers,file = paste("C",i,".markers_FINAL",sep=""),row.names = F,col.names = T,quote=F)
}

#-------------------labeling each cluster
#ref: http://mousebrain.org/celltypes/
# https://www.abcam.com/neuroscience/neural-markers-guide

combined <- RenameIdents(combined, `0` = "Satellite_Glial", `1` = "Endothelial_cell", `2` = "Macrophage", 
                                `3` = "Mature_Neuron", `4` = "Schwann_Cells")
DimPlot(combined, label = TRUE)
DimPlot(combined, reduction = "umap", split.by = "Dataset",label = TRUE,
        label.size = 3,repel = TRUE,ncol = 2)

#----------------Identify differential expressed genes across conditions
combined$celltype.Dataset <- paste(Idents(combined), combined$Dataset, sep = "_")
combined$celltype <- Idents(combined)
Idents(combined) <- "celltype.Dataset"

#will combine both exp1 and exp-4 control only for all downstream analysis
#-----------------------Control vs Ly-treated
ident.Ly= c("Endothelial_cell_Ly","Satellite_Glial_Ly","Macrophage_Ly","Mature_Neuron_Ly","Schwann_Cells_Ly")
exp4_control=c("Endothelial_cell_exp4-control","Satellite_Glial_exp4-control","Macrophage_exp4-control","Mature_Neuron_exp4-control","Schwann_Cells_exp4-control")
exp1_control=c("Endothelial_cell_exp1-control","Satellite_Glial_exp1-control","Macrophage_exp1-control","Mature_Neuron_exp1-control","Schwann_Cells_exp1-control")

file_names_Ly=paste(ident.Ly,"_control_Final",sep = "")
for (i in 1:5){
    results <- FindMarkers(combined,  ident.1 = ident.Ly[i],ident.2 = c(exp4_control[i],exp1_control[i]), verbose = FALSE)
    write.table(results,file=file_names_Ly[i],quote = F,row.names = T,col.names = T)
}
#-----------------------Control vs Mirin-treated
ident.mirin= c("Endothelial_cell_mirin","Satellite_Glial_mirin","Macrophage_mirin","Mature_Neuron_mirin","Schwann_Cells_mirin")
file_names_mirin=paste(ident.mirin,"_control_Final",sep = "")
for (i in 1:5){
  results <- FindMarkers(combined,  ident.1 = ident.mirin[i],ident.2 = c(exp4_control[i],exp1_control[i]), verbose = FALSE)
  write.table(results,file=file_names_mirin[i],quote = F,row.names = T,col.names = T)
}

#---Mirin vs Ly: endothelial and satelite cells
file_names_mirin_Ly=paste(ident.mirin,"_Ly_Final",sep = "")
for (i in 1:5){
  results <- FindMarkers(combined,  ident.1 = ident.mirin[i],ident.2 = ident.Ly[i], verbose = FALSE)
  write.table(results,file=file_names_mirin_Ly[i],quote = F,row.names = T,col.names = T)
}


##Volcano plot in terminal coz Rstudio did not work
library(EnhancedVolcano)
file.lists=c("Endothelial_cell_Ly_control_Final","Satellite_Glial_Ly_control_Final","Macrophage_Ly_control_Final",
             "Mature_Neuron_Ly_control_Final","Schwann_Cells_Ly_control_Final",
             "Endothelial_cell_mirin_control_Final","Satellite_Glial_mirin_control_Final","Macrophage_mirin_control_Final",
             "Mature_Neuron_mirin_control_Final","Schwann_Cells_mirin_control_Final")
for (i in 3:5){
  data=read.table(file.lists[i],header=T)
  jpeg(paste(file.lists[i],".jpg",sep = ""))
  EnhancedVolcano(data,
                  lab = rownames(data),
                  x = 'avg_logFC',
                  y = 'p_val_adj',
                  xlim = c(-2.5,2.5),
                  xlab = 'logFC',
                  title = file.lists[i],
                  pCutoff = 0.01,
                  FCcutoff = 0.5,
                  pointSize = 3.0,
                  labSize = 3.0,
                  labCol='black',
                  labFace='bold',
                  colAlpha = 3/5,
                  legend=c('NS','Log (base 2) fold-change','P value',
                           'P value & Log (base 2) fold-change'),
                  legendPosition = 'top',
                  legendLabSize = 14,
                  legendIconSize = 5.0,
                  drawConnectors = TRUE,
                  widthConnectors = 0.2,
                  colConnectors = 'grey30')
  dev.off() 
}
#-----------------------Signature of different treatments from most afftected genes? (by cell type? together?)-----------------


##--finding signature
genes=read.table("macrophage_mirin",header=T)
sub=subset(combined, subset = celltype == "Macrophage")
test=subset(sub, idents = c("Macrophage_exp4-control", "Macrophage_Ly","Macrophage_mirin"))
sub_genes=c("CTSB","CTSL","CTSD","CTSK","CD68","CD63","CTSA","GM2A","PSAP",
            "NPC2","CTSZ","FCER1G","ATP6VOC","IFI30","B2M","LGLMN") # these are genes in the top three enrichned pathways for Mirin Macrophage: Lysosome, antigen processing and apoptosis
DoHeatmap(test,features = sub_genes,slot = "data")

sat=read.table("satelite_mirin",header=F)
sub_sat=subset(combined, subset = celltype == "Satellite_Glial")
DoHeatmap(sub_sat,features = sat$V1,slot = "data")


##----------------ploting top 3 most significant pathways
pathway=read.table("pathways",header=T)
mirin_pathway=pathway[pathway$Treatment=="Mirin",]
Ly_pathway=pathway[pathway$Treatment=="Ly",]
mirin_pathway$pathway <- factor(mirin_pathway$pathway,levels = c("Cardiac_Muscle_contraction","Oxidative_phosphorylation",
                                                                 "Ribosome","Apoptosis",
                                                                 "Antigen_Processing_and_presentation","Lysosome","FC-gamma_R-mediated_phagocytosis",
                                                                 "Protein_processing_in_ER","Protein_export"))
label=mirin_pathway$pathway
p1=ggplot(mirin_pathway,aes(x=Cell,y=-log10(adj_pvalue),fill=pathway))+
  geom_bar(stat="identity", position=position_dodge(),width = 0.8)+
  geom_text(aes(label=label),position = position_dodge(0.6),size=3, hjust = "inward")+
  scale_y_continuous(breaks=seq(0,38,5))+geom_hline(yintercept=1.30103, linetype="dashed", 
                                                    color = "red", size=0.3)+
  coord_flip()

p2=ggplot(Ly_pathway,aes(x=Cell,y=log,fill=pathway))+
  geom_bar(stat="identity", position=position_dodge(),width = 0.2)+ylab("-log10(adj_pvalue)")+
  geom_text(aes(label=pathway),position = position_dodge(0.2),size=3, hjust = "inward")+
  scale_y_continuous(breaks=seq(0,38,5))+geom_hline(yintercept=1.30103, linetype="dashed", 
                                                   color = "red", size=0.3)+ylim(0,38)+coord_flip()

ggarrange(p2,p1,legend = "none",align = "h",labels = c("A.Ly-treated","B.Mirin-treated"))

#-----Test of differential proportion across conditions
##--------freq cells
(freq_table <- prop.table(x = table(combined@active.ident, combined@meta.data[, "Dataset"]),
                         margin = 2))
barplot(height = freq_table)
freq=as.data.frame(freq_table)
freq <- freq %>%
  spread(key = Var2, value = Freq)
freq$control_avg=(freq$`exp1-control`+freq$`exp4-control`)/2
write.table(freq,file="freq_cells",quote = F,row.names = F,col.names = T)

(table <- table(combined$celltype, combined$Dataset))
df_ori=as.data.frame(table)
df <- as.data.frame(table) %>%
  spread(key = Var2, value = Freq)
df$control_avg=(df$`exp1-control`+df$`exp4-control`)/2
colsum=apply(df[-1], 2, sum)
df=rbind(df[-1],colsum)
write.table(df,file="raw_count_cells",quote = F,row.names = F,col.names = T)

#--test of proportion
(pro.control.Ly=prop.test(x=c(18,305.5),n=c(2283,4612.5),alternative = "less",correct=TRUE))
pro.control.Ly$p.value
(Mirin.endo=prop.test(x=c(336,825.5),n=c(2625,4612.5),alternative = "less",correct=TRUE))

#----------------plotting proportion of cells in to barplots
#load packages
library(ggplot2)
install.packages("tidyr",dep=T)
library(tidyr)
library(plyr)

names(df_ori)= c("Cell_type","Data_set","Number of Cells")
# Get the levels for type in the required order
df_ori$Cell_type= factor(df_ori$Cell_type, levels = c("Satellite_Glial","Endothelial_cell","Macrophage","Mature_Neuron","Schwann_Cells"))
tbl = arrange(df_ori, df_ori$Data_set, desc(df_ori$Cell_type))

# Calculate the percentages
tbl_n = ddply(tbl, .(Data_set), transform, percent = `Number of Cells`/sum(`Number of Cells`) * 100)
# Format the labels 
tbl_n$label = paste0(sprintf("%.2f", tbl_n$percent), "%")

##plotting
ggplot(tbl_n, aes(x = Data_set, y =  Number.of.Cells, fill = Cell_type)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  geom_text_repel(aes(label=label),position = position_stack(vjust = 0.5), size = 3.5)
  

