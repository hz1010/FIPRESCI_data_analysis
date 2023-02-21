library(Matrix)
data=readRDS("TSS_Proportion_mat.rds")
tss_gene=readRDS("gene_by_tss_order.rds")
unique_gene=readRDS("unique_gene.rds")
gene_mat=readRDS("gene_mat.rds")
data=ceiling((data)*100+1)
gene_p=vector(mode="numeric",length=length(unique_gene))
gene_p=rep(1,length(unique_gene))
for (i in 1:length(unique_gene))
{
    idx=which(tss_gene==unique_gene[i])
    if (length(idx)==1)
    {
        gene_p[i]=1
        next
    }
    sub_data=data[idx,]
    ll=c()
    for (l in 1:ncol(gene_mat))
    {
        if ((quantile(gene_mat[,l],0.5)<gene_mat[unique_gene[i],l]) & (quantile(gene_mat[,l],0.99)>gene_mat[unique_gene[i],l]))
        {
            ll=c(ll,l)
        }
    }
    if(length(ll)<=1){next}
    sub_data=sub_data[,ll]
    sub_data=t(sub_data)
    print(sub_data)
    p_list=vector(mode="numeric",length=nrow(sub_data)*(nrow(sub_data)-1)/2)
    s=0
    for (j in 1:(nrow(sub_data)-1))
    {   
        for (k in (j+1):nrow(sub_data))
        {
            s=s+1
            p_list[s]=fisher.test(sub_data[c(j,k),],simulate.p.value=TRUE)$p.value
           
        }
    }
    print(s)
    print(length(p_list))
    p_list_a=p.adjust(p_list,method="BH")
    gene_p[i]=min(p_list_a)
}
names(gene_p)=unique_gene
saveRDS(gene_p,file="gene_p_value.rds")
