library(monocle3)
library(SeuratWrappers)
library(Seurat)
data=readRDS("~/LT_20dim_SCT_res04.rds")
DefaultAssay(data)<-"RNA"
idx=which(data@meta.data$prediction.score.max>0.4)
data=data[,idx]


brain_related=c("Neural Tube","Radial glia","Neural progenitor cells","Granule neurons","Postmitotic premature neurons","Inhibitory interneurons","Cholinergic neurons",
"Excitatory neurons","Inhibitory neuron progenitors","Inhibitory neurons","Sensory neurons","Premature oligodendrocyte","Melanocytes","Ependymal cell","Notochord cells",
"Isthmic organizer cells","Oligodendrocyte Progenitors")

idx=which((data@meta.data$predicted.id %in% brain_related))
data=data[,idx]

sc=read.table("umap.txt") ##umap.txt is PAGA umap from scanpy script
data@reductions$umap@cell.embeddings[,1]=sc$V1
data@reductions$umap@cell.embeddings[,2]=sc$V2
cds=as.cell_data_set(data)

cds=cluster_cells(cds,resolution=0.00025)
saveRDS(cds,file="seurat_cds1_bor.rds")


cds <- learn_graph(cds)

idx=which(cds@colData$predicted.id=="Radial glia")
cell=colnames(cds)[idx]

cds <- order_cells(cds, root_cells=cell)
p2=plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_groups_by_cluster=TRUE,
           graph_label_size=1.5)
p2
saveRDS(cds,file="seurat_cds2.rds")

