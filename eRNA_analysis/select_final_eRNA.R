con=c("Random","PolyT","Mix")
for (i in 1:length(con)){
print(con[i])
eRNA=read.table(paste0("~/",con[i],"/Overlap_",con[i],".bed"))
RPKM=readRDS(paste0("~/",con[i],"/RPKM_",con[i],".rds"))
idx1=which(RPKM>1)
test=readRDS(paste0(con[i],"test_BH.rds"))
idx2=which(test<0.01)
eRNA=eRNA[intersect(idx1,idx2),]
print(length(idx1))
print(length(idx2))
print(nrow(eRNA))

write.table(eRNA,file=paste0("./",con[i],"eRNA.bed"),sep="\t",row.name=F,col.name=F,quote=F)
}
