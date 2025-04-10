# Load packages
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

# Install data
# Output of spaceranger pipeline
#InstallData("stxBrain")

# Load data
brain <- LoadData("stxBrain", type = "anterior1")

# STEP: Check content of Seurat object
brain@meta.data %>% head()
brain@assays
DefaultAssay(brain)
brain@meta.data$orig.ident %>% table()
brain@meta.data$slice %>% table()
brain@meta.data$region %>% table()

# STEP: Investigate data
VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
VlnPlot(brain, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()

# Spatial plots 
SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right") 

# STEP: Pre-processing 
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)

# Gene expression visualization
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))

# STEP: Dimensionality reduction, clustering, and visualization
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
ElbowPlot(brain)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

# Check number of cells in each cluster 
brain@meta.data$seurat_clusters %>% table()

# STEP: Visualize all clusters 
p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2

# Visualize individual clusters 
SpatialDimPlot(brain, 
               cells.highlight = CellsByIdentities(object = brain, idents = c(2, 1, 4, 3, 5, 8)), 
               facet.highlight = TRUE, 
               ncol = 3)

# STEP: Identification of Spatially Variable Features
de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)

SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)
