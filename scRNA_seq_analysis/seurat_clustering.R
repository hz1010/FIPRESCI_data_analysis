library(Seurat)
data=readRDS("merged_gene_mat.rds")
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
#all.genes <- rownames(data)
data <- ScaleData(data)
print("1")
data <- RunPCA(data, features = VariableFeatures(object = data))
data <- RunTSNE(data, dims = 1:20, tsne.method = "FIt-SNE", nthreads = 4, max_iter = 2000)
saveRDS(data,file="after_cluster_20PC.rds")
library(ggplot2)
p1 <- DimPlot(data, reduction = "tsne", pt.size = 0.1) + ggtitle(label = "FIt-SNE")
p1


