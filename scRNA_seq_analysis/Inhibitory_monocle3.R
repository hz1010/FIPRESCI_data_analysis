library(monocle3)
library(ggplot2)

data=readRDS("monocle_seurat_cds.rds") ##monocle3 data in seurat_object
Inhibi=c("Neural Tube","Neural progenitor cells","Radial glia","Inhibitory interneurons","Inhibitory neuron progenitors","Inhibitory neurons") ##Inhibitory ralated neurons
idx=which(data@colData$predicted.id %in% Inhibi)
cell_name=colnames(data)[idx]
time=data@principal_graph_aux$UMAP$pseudotime[cell_name]
data=data[,idx]
df=data.frame(celltype=data@colData$predicted.id,pseudo_time=time)
levels=c("Neural Tube","Radial glia","Neural progenitor cells","Inhibitory neuron progenitors","Inhibitory interneurons","Inhibitory neurons")
df$celltype=factor(df$celltype,levels=levels)

p=ggplot(df, aes(x = celltype, y = pseudo_time))+geom_boxplot(aes(fill = celltype))+theme_bw()+
theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.grid=element_blank())
p

p1=plot_cells(data, color_cells_by = "pseudotime",label_cell_groups = FALSE,label_groups_by_cluster = FALSE,
label_branch_points = FALSE,
       label_roots = FALSE,
       label_leaves = FALSE)
p1

p2=plot_cells(data,
           color_cells_by = "predicted.id",
           label_groups_by_cluster=FALSE,
           group_label_size=5,
           label_branch_points=FALSE,
           label_leaves=FALSE,
           label_roots = FALSE,
           label_principal_points = FALSE,
           rasterize=TRUE
        )
p2

saveRDS(data,file="Inhibi.rds")
