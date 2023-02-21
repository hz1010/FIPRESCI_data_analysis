## following guidelines from https://github.com/s175573/GIANA.

pdf("TCR_TREE.pdf",height=15)
library(ape)
##colors
COLs= c("dodgerblue2","gold", # red
        "green4",
        "#6A3D9A", # purple
        "#FF7F00", # orange
        "cyan",
        "skyblue2","#FB9A99", # lt pink
        "palegreen2",
        "#CAB2D6", # lt purple
        "#FDBF6F", # lt orange
        "gray70", "khaki2",
        "maroon","orchid1","deeppink1","blue1","steelblue4",
        "darkturquoise","green1","yellow4","yellow3",
        "darkorange4","brown")

ffc='../aa_pre_both_all--RotationEncodingBL62.txt'  ## GIANA standard output
ffm='../aa_pre_both_all--RotationEncodingBL62.txt_EncodingMatrix.txt'  ## Output matrix with -M option
cluster=c()
d1=read.table(ffc, header=F,sep='\t',stringsAsFactors = F)
Mat=read.table(ffm, header=F,sep='\t',stringsAsFactors = F)
cluster=c(5,6,27)
idx=which(d1[,2] %in% cluster)
print(unique(d1[,2]))
#d1=d1[idx,]
print(unique(d1[,2]))
Mat=Mat[,c(1,4:99)]
rownames(Mat)=Mat[,1]
Mat=Mat[,2:97]
rownames(Mat[d1[,1],])

d1=cbind(d1, Mat[d1[,1],])
#rownames(d1)=rownames(Mat[d1[,1],])
Mat.cor=cor(t(d1[, 5:100]))
rownames(Mat.cor)
Mat.dist=sqrt(1-Mat.cor)
str(Mat.dist)
tr = nj(Mat.dist)
str(tr)
tr$tip.label = gsub('\\.[1-9]{1}$','', tr$tip.label)
#par(mar=c(0,0,0,0))
table(rownames(Mat[d1[,1],])==tr$tip.label)
col.gr=COLs[c(1,10,12,14)]
par(mar=c(0,0,0,0))
plot.phylo(tr, font=2, cex=0.9,
               ## color the CDR3s according to antigen-specificity
           align.tip.label=TRUE, x.lim=c(0,1.1))
