# How to use the package
The following text provides a sample.
```
# Here is a sample
library(Seurat)
library(scAnnotationBot)

# Seurat tutorial Data get from: https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
pbmc.data = Read10X(data.dir = "./filtered_gene_bc_matrices/hg19/")
pbmc = CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc = subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc = NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc = FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc = ScaleData(pbmc)
pbmc = RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pcNums = 1:10
pbmc = FindNeighbors(pbmc, dims = pcNums)
pbmc = FindClusters(pbmc, graph.name = "RNA_snn", resolution = 0.5)
pbmc = RunUMAP(pbmc, dims = pcNums)
DimPlot(pbmc, reduction = "umap", label = TRUE)
# FindAllMarkers
pbmc.markers = FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# ERNIEBot
account = list(client_id = 'MTtw1WwBdn0GCZ5xHTdeoEPk', client_secret = '18mxgswFWmmxnUIIQLL0FGL9CfDj6aWW') # This account for test only
cellType_ERNIEBot = ERNIEBotCellType(account = account, 
                                     model = "completions_pro", 
                                     inputTable = pbmc.markers, tissueName = 'Human peripheral blood mononuclear cell', 
                                     useGeneNumber = 10)
pbmc@meta.data$celltype_ERNIEBot = as.factor(cellType_ERNIEBot[as.character(Idents(pbmc))])
DimPlot(pbmc, 
        label = T, 
        group.by='celltype_ERNIEBot', 
        reduction = "umap")
# DashScope
cellType_DashScope = DashScopeCellType(api_key = 'sk-ce382596cab642bfa80b0907af0e7daf', # This key for test only
                                       model = "qwen-max-longcontext", 
                                       inputTable = pbmc.markers, tissueName = 'Human peripheral blood mononuclear cell', 
                                       useGeneNumber = 10)
pbmc@meta.data$celltype_DashScope = as.factor(cellType_DashScope[as.character(Idents(pbmc))])
DimPlot(pbmc, 
        label = T, 
        group.by='celltype_DashScope', 
        reduction = "umap")
```
