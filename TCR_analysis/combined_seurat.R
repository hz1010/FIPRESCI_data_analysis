##following scRepertoire guidelines (https://ncborcherding.github.io/vignettes/vignette.html)
pdf("combined2.pdf")
library(Seurat)
library(scRepertoire)
library(circlize)
library(scales)
colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF",
              "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
              "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))
seurat=readRDS("donor_seurat.rds")
ind=seurat@meta.data$donor
split=strsplit(seurat@meta.data$donor,"-")
for (i in 1:length(ind))
{
    
    if(split[[i]][2]=="healthy")
    {
        ind[i]=split[[i]][2]
    } 
    else
    {
        ind[i]=paste0(split[[i]][2],"-",split[[i]][3])
    }
}
#print(table(ind))
seurat@meta.data$ill=ind
Idents(object=seurat)=seurat@meta.data$predicted.id
clonalOverlay(seurat,freq.cutpoint=1,bins=10,reduction="ref.umap")#+guides(color=FALSE)+scale_fill_gradient()
clonalOverlay(seurat,freq.cutpoint=1,bins=10,reduction="ref.umap",facet="donor")#+guides(color=FALSE)+scale_fill_gradient()
clonalOverlay(seurat,freq.cutpoint=1,bins=10,reduction="ref.umap",facet="ill")#+guides(color=FALSE)+scale_fill_gradient()

circles=getCirclize(seurat,groupBy="donor")
saveRDS(circles,file="circles.rds")
grid.cols=scales::hue_pal()(length(unique(seurat@meta.data$donor)))
print(unique(seurat@meta.data$donor))
circlize::chordDiagram(circles,self.link=1,grid.col=grid.cols)
title(main="all-donors")

circles=getCirclize(seurat,groupBy="predicted.id")
grid.cols=scales::hue_pal()(length(unique(seurat@meta.data$predicted.id)))
circlize::chordDiagram(circles,self.link=1,grid.col=grid.cols)
title(main="all-celltype")

sample=unique(seurat@meta.data$donor)
for (i in 1:length(sample))
{
    idx=which(seurat@meta.data$donor==sample[i])
    tmp=seurat[,idx]
    circles=getCirclize(tmp,groupBy="predicted.id")
    grid.cols=scales::hue_pal()(length(unique(tmp@meta.data$predicted.id)))
    circlize::chordDiagram(circles,self.link=1,grid.col=grid.cols)
    title(main=sample[i])
    print(unique(tmp@meta.data$predicted.id))
    print(grid.cols)

}
