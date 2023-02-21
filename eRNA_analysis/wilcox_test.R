##BH ajusted wilcox test for core peak region and backgroud region

library(Matrix)
con=c("Random","PolyT","Mix")
for (j in 1:length(con)){
data=read.table(paste0("mat_",con[j],".txt"))
data=as.matrix(data)
dim(data)
test=vector(mode="numeric",length=nrow(data))
backgroud=vector(mode="numeric",length=ncol(data))
for (i in 1:length(backgroud)){backgroud[i]=median(data[,i])}
for (i in 1:nrow(data)){
test[i]=wilcox.test(data[i,150:250],backgroud[c(1:149,251:400)],alternative="greater")$p.value
}
test=p.adjust(test,"BH")
saveRDS(test,file=paste0(con[j],"test_BH.rds"))

