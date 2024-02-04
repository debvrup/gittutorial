# Set the directory

setwd("C:/Users/aruna/Documents/SC-RNA-SEQ")
#*****************************************************
#loading the libraries
library(SingleR)
library(scran)
library(ggplot2)
library(BiocParallel)
library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(hdf5r)
library(ggsci)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(patchwork)
library(glmGamPoi)
library(reticulate)
library(DoubletFinder)

#********************************************************
# 1) Input data and meta data
sparse_matrix <- Seurat::Read10X_h5(filename = "UVM_GSE169609_expression.h5", use.names = T)
metadata      <- readr::read_tsv("UVM_GSE169609_CellMetainfo_table.tsv")

# Create a Seurat object
# include only genes that are are expressed in 3 or more cells and cells with complexity of 200 genes or more
seurat_obj <- CreateSeuratObject(counts = sparse_matrix, meta.data = metadata ,project = "TISCH2",min.cells = 3, min.features = 200 )
seurat_obj
str(seurat_obj)

#remove data to save memory
#sparse_matrix <- NULL

#separating meta data
meta <- seurat_obj@meta.data
head(meta)
summary(meta$nCount_RNA)
summary(meta$nFeature_RNA)

#Adding Mitochondrial genes and ribosomal proteins
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")


#violin plots-VlnPlot
VlnPlot(seurat_obj, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))

#_Plot Metadata Features_
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.rb")


#*************************************************

#2)__QC__improving data
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
str(seurat_obj)

#*************************************************

#3)Normalization and filtering
seurat_obj <- NormalizeData(seurat_obj,normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj <- FindVariableFeatures(seurat_obj,selection.method = "vst", nfeatures = 2000)

#Identify the 10 most highly variable genes:
top10 <- head(VariableFeatures(seurat_obj), 10)
top10
#Plot variable features with and without labels:
plot1 <- VariableFeaturePlot(seurat_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot1 + plot2

#______Doublet Detection________
set.seed(123)  # Set a seed for reproducibility
seurat_obj <- doubletFinder(seurat_obj, PCs = 1:20, pK = 30, nExp = 500)

# Identify doublet scores and plot
doublet_scores <- seurat_obj$doublet_scores
hist(doublet_scores, breaks = 50, main = "Doublet Scores")

# Remove identified doublets
seurat_obj <- subset(seurat_obj, cells = !seurat_obj$doublet)

# Update metadata after removing doublets
metadata <- seurat_obj@meta.data


#*************************************************

#4)Scaling the data and dimensional reduction
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)
seurat_obj <- RunPCA(seurat_obj,features = VariableFeatures(object = seurat_obj))
print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)

#plot
VizDimLoadings(seurat_obj, dims = 1:9, reduction = "pca") & 
    theme(axis.text=element_text(size=5), axis.title=element_text(size=8,face="bold"))

#plot
# allows for easy exploration of the primary sources of heterogeneity in a dataset
# and can be useful when trying to decide which PCs to include for further downstream analyses
DimHeatmap(seurat_obj, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seurat_obj, dims = 1:6, nfeatures = 20, cells = 500, balanced = T)

#plot
DimPlot(seurat_obj, reduction = "pca")
DimPlot(seurat_obj, reduction = "pca") + NoLegend()
ElbowPlot(seurat_obj)

# ___to dertermine "dimensionality" of the dataset -------
# essentially determine how many PCs to consider - we would ideally want to consider PCs that show maximum variations

# JackStraw Procedure!
# identify ‘significant’ PCs as those who have a strong enrichment of low p-value features.
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time

seurat_obj <- JackStraw(seurat_obj, num.replicate = 100)
seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:20)

JackStrawPlot(seurat_obj, dims = 1:15)
# The JackStrawPlot() function provides a visualization tool for comparing the distribution of p-values for each PC with a uniform distribution (dashed line). 
# ‘Significant’ PCs will show a strong enrichment of features with low p-values (seurat_objlid curve above the dashed line).
# An alternative heuristic method generates an ‘Elbow plot’: a ranking of principle components based on the percentage of variance explained by each one (ElbowPlot() function).
ElbowPlot(seurat_obj)

#********************************************************

#5)FindNeighbours For Clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10, reduction = "pca")

#clustering
#seurat_obj <- SCTransform(seurat_obj)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5,graph.name = "RNA_snn")

#generate UMAP reduced dimensionality representation
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10, verbose = F)
table(seurat_obj@meta.data$seurat_clusters)

DimPlot(seurat_obj,label.size = 4,repel = T,label = T) #Dimplot UMAP By default

#control for clustering resolution and other possible artifacts, we will take a close look at two minor
#cell populations: 1) dendritic cells (DCs), 2) platelets, aka thrombocytes. Let’s visualise two markers
#for each of this cell type: LILRA4 and TPM2 for DCs, and PPBP and GP1BB for platelets.

#FeaturePlot(seurat_obj, features = c("LILRA4", "TPM2", "PPBP", "GP1BB"))

#visualize other confounders:
FeaturePlot(seurat_obj, features = "percent.mt") & theme(plot.title = element_text(size=10))
FeaturePlot(seurat_obj, features = "nFeature_RNA") & theme(plot.title = element_text(size=10))
DimPlot(seurat_obj,label.size = 4,repel = T,label = T)

#************************************************************

#6)cell cycle scores
cc.genes.updated.2019
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes)
table(seurat_obj[[]]$Phase)

# Visualize the distribution of cell cycle markers across
RidgePlot(seurat_obj, features = c("SLBP", "RRM1", " CKS2", " POLR1B"), ncol = 2)
# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by
# phase
seurat_obj <- RunPCA(seurat_obj, features = c(s.genes, g2m.genes))
DimPlot(seurat_obj)

#clusters defined by any of the technical differences
FeaturePlot(seurat_obj,features = "percent.mt",label.size = 4,repel = T,label = T) & 
  theme(plot.title = element_text(size=10))
VlnPlot(seurat_obj,features = "percent.mt") & theme(plot.title = element_text(size=10))
FeaturePlot(seurat_obj,features = "percent.rb",label.size = 4,repel = T,label = T) & theme(plot.title = element_text(size=10))
VlnPlot(seurat_obj,features = "percent.rb") & theme(plot.title = element_text(size=10))
VlnPlot(seurat_obj,features = c("nCount_RNA","nFeature_RNA")) & 
  theme(plot.title = element_text(size=10))

FeaturePlot(seurat_obj,features = c("S.Score","G2M.Score"),label.size = 4,repel = T,label = T) & 
  theme(plot.title = element_text(size=10))
VlnPlot(seurat_obj,features = c("S.Score","G2M.Score")) & 
  theme(plot.title = element_text(size=10))

#************************************************************

#7)SCTransform normalization and clustering
seurat_obj <- SCTransform(seurat_obj, method = "glmGamPoi", ncells = 8824, 
                    vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = F)
seurat_obj

#PCA, UMAP, and clustering
seurat_obj <- RunPCA(seurat_obj, verbose = F)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = F)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = F)
seurat_obj <- FindClusters(seurat_obj, verbose = F)
table(seurat_obj[[]]$seurat_clusters)


DimPlot(seurat_obj, label = T)

#identification of rare cell populations,fewer neighbors in the KNN graph, Leiden algorithm 
#seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, k.param = 15, verbose = F)
#seurat_obj <- FindClusters(seurat_obj, verbose = F, algorithm = 4, resolution = 0.9)

table(seurat_obj[[]]$seurat_clusters)

DimPlot(seurat_obj, label = T)

#****************************************************************

#8)Differential expression and marker selection || CLUSTER BIOMARKERS
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)
#Find markers for every cluster compared to all remaining cells, report only the positive ones
all.markers <- FindAllMarkers(seurat_obj, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
seurat_obj <- FindMarkers(seurat_obj, group.by = "condition")

# ___find all markers of cluster 1 --------
cluster1.markers <- FindMarkers(seurat_obj, ident.1 = 2, min.pct = 0.25)
head(cluster1.markers, n = 5)
# ___find all markers distinguishing cluster 5 from clusters 0 and 3 --------
cluster5.markers <- FindMarkers(seurat_obj, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

#### experiment labels on a UMAP visualization
#DimPlot(seurat_obj, reduction = "umap", group.by = "Celltype..major.lineage.", cells = sample(colnames(seurat_obj))) + scale_color_igv()

### cell type labels on a UMAP visualization
#DimPlot(seurat_obj, reduction = "umap", group.by = "experiment", cells = sample(colnames(seurat_obj))) +
 # scale_color_nejm()

# Save the aggregated markers as a CSV file in the current working directory
write.csv(all.markers, file = "seurat_markers.csv", row.names = FALSE)

#****************************************************************

#9) Cell type annotation using SingleR
monaco.ref <- celldex::MonacoImmuneData()
#hpca.ref <- celldex::HumanPrimaryCellAtlasData()
#dice.ref <- celldex::DatabaseImmuneCellExpressionData()

#converting Seurat object to single cell experiment (seurat_obj)
seurat_obj <- as.SingleCellExperiment(DietSeurat(seurat_obj))
seurat_obj

monaco.main <- SingleR(test = seurat_obj,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.main)
monaco.fine <- SingleR(test = seurat_obj,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.fine)
# hpca.main <- SingleR(test = seurat_obj,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.main)
# hpca.fine <- SingleR(test = seurat_obj,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.fine)
# dice.main <- SingleR(test = seurat_obj,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.main)
# dice.fine <- SingleR(test = seurat_obj,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.fine)

table(monaco.main$pruned.labels)
#table(hpca.main$pruned.labels)
#table(dice.main$pruned.labels)

#finer cell types annotations are you after, the harder they are to get reliably.
#This is where comparing many databases, as well as using individual markers from literature, would be very valuable.
table(monaco.fine$pruned.labels)
# table(hpca.fine$pruned.labels)
# table(dice.fine$pruned.labels)

#add the annotations to the Seurat object metadata
seurat_obj@meta.data$monaco.main <- monaco.main$pruned.labels
seurat_obj@meta.data$monaco.fine <- monaco.fine$pruned.labels

# seurat_obj@meta.data$hpca.main   <- hpca.main$pruned.labels
# seurat_obj@meta.data$dice.main   <- dice.main$pruned.labels
# seurat_obj@meta.data$hpca.fine   <- hpca.fine$pruned.labels
# seurat_obj@meta.data$dice.fine   <- dice.fine$pruned.labels

seurat_obj <- SetIdent(seurat_obj, value = "monaco.fine")

#****************************************************************

#10)___Visualization___
# VlnPlot() (shows expression probability distributions across clusters)
# FeaturePlot() (visualizes feature expression on a tSNE or PCA plot) are our most commonly used visualizations. 
# RidgePlot(), CellScatter(), and DotPlot() as additional methods to view your dataset.
#__feature genes can be found in saved biomarkers files
DimPlot(seurat_obj, label = T , repel = T, label.size = 3) + NoLegend()

str(seurat_obj)
features <- seurat_obj@commands$RunPCA.RNA$features
VlnPlot(seurat_obj, features = c("CA1", "AHSP"))

# ___VlnPlot() - you can plot raw counts as well ---------
VlnPlot(seurat_obj, features = c("CA1", "AHSP"), slot = "counts", log = TRUE)

# ___FeaturePlot()- visualize feature expression in low-dimensional space ---------
FeaturePlot(seurat_obj, features = features[1:5])

# Visualize co-expression of two features simultaneously
FeaturePlot(seurat_obj, features = features[1:2], blend = TRUE)

# ___interactive plots --------
# Include additional data to display alongside cell names by passing in a data frame of
# information Works well when using FetchData
# works only with one feature
plot <- FeaturePlot(seurat_obj, features = c("MCM5"))
HoverLocator(plot = plot, information = FetchData(seurat_obj, vars = c("ident", "PC_1", "nFeature_RNA")))

# ___doHeatmap() --------
top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(seurat_obj, features = top10$gene) + NoLegend()

DoHeatmap(seurat_obj, features = top10$gene, group.by = "cluster") + NoLegend()

# ___RidgePlot() - Visualize single cell expression distribution in each cluster -------
RidgePlot(seurat_obj, features = features[1:5], ncol=2)
#___________________NEED TO BE FIXEDD______________________________________________
# ___Dot plots - the size of the dot corresponds to the percentage of cells expressing the feature --------
# in each cluster. The color represents the average expression level
DotPlot(seurat_obj, features = features[1:5]) + RotatedAxis()
DotPlot(object = seurat_obj, features = c ("SLBP","RRM1"," CKS2"," POLR1B","MCM5") + RotatedAxis())

# ___Single cell heatmap of feature expression -------
DoHeatmap(subset(seurat_obj, downsample = 100), features = features[1:5], size = 3)

#*********************************************************************

#11) saving processed data ------
saveRDS(seurat_obj, file = "seurat_obj.rds")

seurat_obj<-readRDS("seurat_obj.rds")


