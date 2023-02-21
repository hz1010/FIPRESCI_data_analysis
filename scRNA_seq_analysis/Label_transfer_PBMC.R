### following Seurat label transfer guidelines
library(Seurat)
library(SeuratDisk)
library(uwot)
reference=LoadH5Seurat("pbmc_multimodal.h5seurat")
query=readRDS("after_cluster_20PC.rds")
query <- SCTransform(query, verbose = FALSE)
anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  normalization.method = "SCT",
  reduction="pcaproject",
  reference.reduction="spca",
  dims = 1:30
)
saveRDS(anchors,file="anchors.rds")
query <- MapQuery(
  anchorset = anchors,
  query = query,
  reference = reference,
  refdata = "celltype.l2",
  reference.reduction="spca",
  reduction.model = "wnn.umap"
)
saveRDS(query,file="integration.rds")
