library(Matrix)
library(ggplot2)
library(reshape)

##One index represent One Condition
round1=c("GGACGA","TCAGTG","TTGCTC","GCTGAG","GGTGCT","TTAAGA")
con=c("polyT","random","mix","polyT","random","mix")
cell=c("cell","cell","cell","nuclei","nuclei","nuclei")
mat=matrix(data=NA, nrow = 10, ncol = 6, byrow = FALSE, dimnames = NULL)
for (i in 1:6)
{
    print(i)
    for(j in 0:9)
    {
     	path=paste0("/xtdisk/jiangl_group/huangzh/SuperLoading/T22-add/Saturation/Reper/",round1[i],"/size_",j,"/",round1[i],"_",j,"/outs/filtered_feature_bc_matrix")
        print(path)
        data=Read10X(data.dir = path)
        data=CreateSeuratObject(counts = data, min.cells = 0, min.features = 0)
        mat[j+1,i]=mean(data@meta.data$nFeature_RNA)
    }
}
colnames(mat)=round1
rownames(mat)=seq(50,5,by=-5)

df=data.frame(con=mat)
df=melt(df)
df$x=rep(seq(50,5,by=-5),times=6)
df$RT_Condition=rep(con,each=10)
df$Preparation=rep(cell,each=10)
df$CT=paste0(df$cell,"_",df$RT)


p=ggplot(data = df, mapping = aes(x = x, y = value,color=RT_Condition,linetype=Preparation)) +
geom_line()+geom_point()+theme_bw()+labs(x="Sequencing depth per cell\n(thound read pairs)",y="Average number of genes detected")+
ggtitle("Sensitivity per cell")+theme(plot.title = element_text(hjust = 0.5))

p1=ggplot(data = df, mapping = aes(x = x, y = value,color=RT_Condition)) +
geom_line()+geom_point()+theme_bw()+labs(x="Sequencing depth per cell\n(thound read pairs)",y="Average number of genes detected")+
ggtitle("Sensitivity per cell")+theme(plot.title = element_text(hjust = 0.5))

