##this script is used to select ATAC peaks whose distance to any TSS is greater than 2000bp (distal peaks)
data=read.table("Hela_dis2k.bed")
idx=which(data$V18>2000)
print(length(idx))
data1=data[idx,]
data1=data1[,c(1,2,3,4)]
write.table(data1,file="enhancer_tss2k.bed",row.names=FALSE,col.names=FALSE,quote=F,sep="\t")

