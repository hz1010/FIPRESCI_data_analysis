library(RColorBrewer)
temp <- RColorBrewer::brewer.pal(11, "Spectral")
temp[6] <- "gold"
rbPal <- colorRampPalette(temp)
colors=rev(rbPal(50))
library(Seurat)
library(Matrix)
library(ggplot2)
pdf("multiplet.pdf")

##classify mouse or human gene
data=readRDS("mixture_filter.rds")
mat=data@assays$RNA@counts
name=rownames(data)
name=substr(name,1,4)
idx1=which(name!="mm10")
name[idx1]="hg38"
idx2=which(name=="mm10")
mat=as.matrix(mat)
mat_hg38=mat[idx1,]
mat_mm10=mat[idx2,]

hg38_UMI_count=colSums(mat_hg38)
mm10_UMI_count=colSums(mat_mm10)

##define human and mouse droplet
species=mm10_UMI_count/(hg38_UMI_count+mm10_UMI_count)
idx_mm10=which(species>=0.9)
idx_multiplet=which((species>0.1) & (species<0.9))
idx_hg38=which(species<=0.1)
species[idx_mm10]="mm10"
species[idx_hg38]="hg38"
species[idx_multiplet]="multiplet"


##plot droplet based on human genome (x axis) and mouse genome (y axis)
df=data.frame(hg38=hg38_UMI_count,mm10=mm10_UMI_count,type=species)
ggplot(df,aes(x=hg38,y=mm10,colour=type))+geom_point(size=0.5)+geom_hline(yintercept=0)+geom_vline(xintercept=0)+theme(panel.grid=element_blank())+theme_bw()+theme(panel.border = element_blank())


