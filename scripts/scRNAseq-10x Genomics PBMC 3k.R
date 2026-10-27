library(Seurat)

download.file("https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz", 
              destfile = "pbmc3k.tar.gz")

untar("pbmc3k.tar.gz")

pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")

# percentage of reads that map to the mitochondrial genome
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Subset the data to remove dead cells and doublets
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Log-normalize the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Selection and Scaling
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

pbmc <- ScaleData(pbmc)

# PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

# UMAP
pbmc <- FindNeighbors(pbmc, dims = 1:10)

pbmc <- FindClusters(pbmc, resolution = 0.5)

pbmc <- RunUMAP(pbmc, dims = 1:10)

DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("UMAP of 3k PBMCs")

# Annotation
# Cluster 1
cluster1_markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1_markers, n = 5)

# for EVERY cluster automatically
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Extract the top 2 genes per cluster
library(dplyr)
pbmc.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 2, wt = avg_log2FC)