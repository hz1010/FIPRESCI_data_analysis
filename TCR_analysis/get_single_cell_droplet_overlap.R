library(Seurat)
library(Matrix)

##prediction score > 0.4
data=readRDS("integration2.rds")
idx=which(data@meta.data$predicted.id.score>=0.4)
data=data[,idx]

##VDJ-seq barcode match
VDJ=read.csv("~/outs/all_contig_annotations.csv")
VDJ_cell=unique(VDJ$barcode)
RNA_cell=colnames(data)
RNA_cell=substr(RNA_cell,8,25)
tmp=table(RNA_cell)
barcode1=names(tmp[tmp==1])
VDJ_1cell=VDJ_cell[VDJ_cell %in% barcode1]
saveRDS(VDJ_1cell,file="single_cell_droplet_overlap.rds")
