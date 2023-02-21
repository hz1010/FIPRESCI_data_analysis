##combined with scRepertoire guidelines (https://ncborcherding.github.io/vignettes/vignette.html)

pdf("clonetype.pdf")#,width=20,height=30)
library(ggplot2)
library(Seurat)
library(scRepertoire)
colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", "#7301A8FF", "#9C179EFF",
              "#BD3786FF", "#D8576BFF","#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))

##define main T cell
T=c("CD8 Naive","CD8 TEM","CD4 Naive","CD4 TCM","Treg","CD4 TEM")

##read data 
overlap_cell=readRDS("single_cell_droplet_overlap.rds")
data=readRDS("integration.rds")
VDJ=read.csv("~/outs/all_contig_annotations.csv")
idx=which(data@meta.data$predicted.id.score>=0.4)
data=data[,idx]
###########################################change sample names
name=data@meta.data$donor
print(table(name))
split=strsplit(name,split="-")
for (i in 1:length(name))
{
    if (split[[i]][2]=="healthy")
    {
        name[i]=paste0(split[[i]][2],"-",split[[i]][1])        
    }
    else
    {
        name[i]=paste0(split[[i]][2],"-",split[[i]][3],"-",split[[i]][1])
    }
}

data@meta.data$donor=name
RNA_cell=colnames(data)
RNA_cell=substr(RNA_cell,8,25)
##########################################get ID. cancer types or health are IDs here 
sample=unique(data@meta.data$donor)
ID=sample
for (i in 1:length(sample))
{   
    if (strsplit(sample,split="-")[[i]][1]=="healthy")
    {
        ID[i]=strsplit(sample,split="-")[[i]][1]
    }
    else
    {
        ID[i]=paste0(strsplit(sample,split="-")[[i]][1],"-",strsplit(sample,split="-")[[i]][2])
    }
}



contig_list=list()
a=0
for (i in 1:length(sample))
{
    
    idx1=which(data@meta.data$donor==sample[i])
    RNA_cell_selected=RNA_cell[idx1]
    overlap=overlap_cell[overlap_cell %in% RNA_cell_selected]
    contig_list[i]=list(VDJ[which(VDJ$barcode %in% overlap),])
    a=a+length(which(unique(VDJ$barcode) %in% overlap))
}
saveRDS(contig_list,file="contig_list.rds")
#saveRDS(contig_list,file="contig_list.rds")

######## get object for scRepertoire input
combined <- combineTCR(contig_list,
                samples = sample,ID=ID,
                 cells ="T-AB")
saveRDS(combined,file="combined.rds")
#order=sample[c(1,2,6,7,8,9,10,11,12,13,14,3,4,5)]

######quantContig for unique clonetype
#quantContig(combined, cloneCall="gene+nt", scale = T)+theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="none")+ggtitle("chain=both")##+scale_x_discrete(limits=order)
#quantContig(combined, cloneCall="gene+nt", scale = T,chain="TRA")+theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="none")+ggtitle("chain=TRA")
#quantContig(combined, cloneCall="gene+nt", scale = T,chain="TRB")+theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="none")+ggtitle("chain=TRB")
#p=quantContig(combined, cloneCall="gene+nt", scale = F, exportTable = T)+theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="none")+ggtitle("chain=both")
#saveRDS(p,file="unique_clotype.rds")
#####abundanceContig plot to show clonetype abundance
#abundanceContig(combined, cloneCall = "gene+nt", scale = F)

#####length distribution of the CDR3 sequences 
#lengthContig(combined, cloneCall="nt", chain = "both")+ggtitle("Both")+theme(legend.text=element_text(size=rel(0.5)))
a1=lengthContig(combined, cloneCall="nt", chain = "TRA")+ggtitle("TRA")+theme(legend.text=element_text(size=rel(0.5)))
lengthContig(combined, cloneCall="nt", chain = "TRB")+ggtitle("TRB")+theme(legend.text=element_text(size=rel(0.5)))
saveRDS(a1,file="lengthcontig.rds")
#####show dynamic changes between different samples
#compareClonotypes(combined, numbers = 2, samples = c("healthy-H1_healthy","healthy-H3_healthy"), 
#                    cloneCall="gene+nt", graph = "alluvial")+theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.text=element_text(size=rel(0.2)))

#####visulize gene usage in different samples
##subset=subsetContig(combined,name="sample",variables=c("H1-healthy","H3-healthy"))
#vizGenes(combined, gene = "V", chain = "TRB", plot = "bar", order = "variance", scale = TRUE)+ theme(axis.text=element_text(size=26))
#vizGenes(combined, gene = "V", chain = "TRB", plot = "heatmap", order = "variance", scale = TRUE)+ theme(axis.title.x = element_text(size = 16))
#vizGenes(combined, gene = "J", chain = "TRB", plot = "bar", order = "variance", scale = TRUE)+ theme(axis.text=element_text(size=26))
#vizGenes(combined, gene = "V", chain = "TRA", plot = "bar", order = "variance", scale = TRUE)+ theme(axis.text=element_text(size=26))
#vizGenes(combined, gene = "J", chain = "TRA", plot = "bar", order = "variance", scale = TRUE)+ theme(axis.text=element_text(size=26))

####show clonetype Abundance 
#clonalHomeostasis(combined,cloneCall="gene+nt")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
#clonalProportion(combined,cloneCall="gene+nt")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
a=clonalProportion(combined,cloneCall="gene+nt")+theme(axis.text.x = element_text(angle = 90, hjust = 1))+theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="none")
mat_melt=a$data
    col <- length(unique(mat_melt$Var2))
    levels=c("healthy-H1_healthy","healthy-H3_healthy","liver-cancer-S1_liver-cancer","liver-cancer-S2_liver-cancer",
"breast-cancer-S3_breast-cancer","breast-cancer-S4_breast-cancer","breast-cancer-S5_breast-cancer","stomach-cancer-S11_stomach-cancer",
"stomach-cancer-S12_stomach-cancer","pancreatic-cancer-S6_pancreatic-cancer",
"pancreatic-cancer-S7_pancreatic-cancer","pancreatic-cancer-S8_pancreatic-cancer","endometrial-cancer-S9_endometrial-cancer","endometrial-cancer-S10_endometrial-cancer")

levels=c("healthy-H1_healthy","healthy-H3_healthy","liver-cancer-S1_liver-cancer","liver-cancer-S2_liver-cancer",
"breast-cancer-S3_breast-cancer","breast-cancer-S4_breast-cancer","breast-cancer-S5_breast-cancer","pancreatic-cancer-S6_pancreatic-cancer",
"pancreatic-cancer-S7_pancreatic-cancer","pancreatic-cancer-S8_pancreatic-cancer",
"endometrial-cancer-S9_endometrial-cancer","endometrial-cancer-S10_endometrial-cancer","stomach-cancer-S11_stomach-cancer","stomach-cancer-S12_stomach-cancer")

mat_melt$Var1=factor(mat_melt$Var1,levels=levels)
    plot <- ggplot(mat_melt, aes(x=Var1, y=value, fill=Var2)) +
        geom_bar(stat = "identity", position="fill",
                    color = "black", lwd= 0.25) +
        scale_fill_manual(name = "Clonal Indices",
                        values = colorblind_vector(col)) +
        xlab("Samples") +
        ylab("Occupied Repertoire Space") +
        theme_classic()+theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="none")+ggtitle(paste0("chain=both,all cell"))
print(plot)

#clonalOverlap(combined,cloneCall="gene",method="overlap")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
#saveRDS(p,file="clonalpropo.rds")
####show clonetype diversity 
#p=clonalDiversity(combined,cloneCall="gene+nt",group.by="sample",x.axis="ID",n.boots=1000)+theme(axis.text.x = element_text(angle = 90, hjust = 1))
#p





