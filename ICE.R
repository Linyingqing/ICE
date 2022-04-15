# preprocess the MAF file
library(data.table)
#获取数据
load('data/new_InfluenceGraph.Rdata')
load('data/Lung.Rdata')###prostate/breast/lung
compartment=read.csv('data/Gen_compartment_together.csv')
#功能交互网络影响图


pro = c(0.2,0.4,0.6,0.8)
res =c()
for (q in 1:4){
for(ii in 1:30){
  res_tem = c()

#随机取出基因连接边
inF=influenceGraph[intersect(colnames(patMutMatrix),row.names(influenceGraph)),intersect(colnames(patMutMatrix),row.names(influenceGraph))]
prob=runif(dim(inF)[1], 0, 1)
inF[which(prob <= pro[q]),which(prob <= pro[q])] = 0
inf <- inF
rm(inF)
inf[lower.tri(inf)] <- t(inf)[lower.tri(inf)]


total_graph = inf 
patMatrix = patMutMatrix
gen_fr=colSums(patMatrix)
gene_name = colnames(patMutMatrix)


#计算突变频率
gene_frequence=gen_fr/nrow(patMatrix)
gene_frequence=cbind(names(gene_frequence),gene_frequence)
total_graph=(total_graph-min(total_graph))/(max(total_graph)-min(total_graph))
tt=list(total_graph=total_graph,frequency=gene_frequence)

seperate = compartment
frequency = tt$frequency

#计算不在细胞隔间和网络中的孤立点
#notcom_net = gene_name[-which(gene_name %in% row.names(influenceGraph))]
#notcom_sep = gene_name[-(which(gene_name %in% seperate[,1])) ]
#notcom = c(notcom_sep,notcom_sep)
#notcom = notcom[!duplicated(notcom)]
#notcom_fre = gene_frequence[which(gene_frequence[,1]%in%notcom),]
#notcom_ei = -log(as.numeric(as.numeric(notcom_fre[,2])))*as.numeric(notcom_fre[,2])
#notcom_fre = cbind(notcom_fre,notcom_ei)

com=unique(seperate[,2])
group=list()

for(i in 1:length(com))#####divide matrix into 11 compartment subgroups将矩阵划分为11个隔室子群
{
  com_i=seperate[which(seperate[,2]==com[i]),1]
  group_i=intersect(colnames(total_graph),com_i)
  group[[i]]=cbind(Gene=group_i,compartment=as.character(com[i]))
}

H_G=c()

CER = c()
for(j in 1:length(group))###calculate the entrophy for each sub-group计算每个子组的熵
{
  sub_group=total_graph[group[[j]][,1],group[[j]][,1]]
  row.names(frequency)=frequency[,1]
  i_entropy=c()
  
  #k=2
  N = length(rownames(sub_group[,]))
  ##遍历基因数
  for(k in 1:nrow(sub_group))
  {
    
    
    n_1 =length(which(sub_group[k,]!=0))####number of the neighbor nodes邻居基因数
    
    #########################below is the cross_entropy of variation frequency
    
    #wij=sub_group[k,which(sub_group[k,]!=0)]
    #n_2 =length(which(sub_group[which(sub_group[k,]!=0),]!=0))+n_1
    
    yj=frequency[colnames(sub_group)[which(sub_group[k,]!=0 )],]
    
    if(is.null(dim(yj))==TRUE)
    {
      #Wi=sum(wij*log(as.numeric(frequency[row.names(sub_group)[k],2])/as.numeric(yj[2]))*as.numeric(frequency[row.names(sub_group)[k],2]))
      CEi=sum(abs(log(as.numeric(frequency[row.names(sub_group)[k],2])/as.numeric(yj[2]))*as.numeric(frequency[row.names(sub_group)[k],2])))
      
    }else
    {
      #=sum(wij*log(as.numeric(frequency[row.names(sub_group)[k],2])/as.numeric(yj[,2]))*as.numeric(frequency[row.names(sub_group)[k],2]))
      CEi=sum(abs(log(as.numeric(frequency[row.names(sub_group)[k],2])/as.numeric(yj[,2]))*as.numeric(frequency[row.names(sub_group)[k],2])))
      
    }
    
    CC = c()
    CC = cbind(CC , row.names(sub_group)[k])
    CC = cbind(CC , CEi)
    CC = cbind(CC , (- log(as.numeric(as.numeric(frequency[row.names(sub_group)[k],2])))*as.numeric(frequency[row.names(sub_group)[k],2])))
      
    CEi = CEi*0.01 + (- log(as.numeric(as.numeric(frequency[row.names(sub_group)[k],2])))*as.numeric(frequency[row.names(sub_group)[k],2]))
    
    Ei = (1+n_1/(N-1))*CEi
    CC = cbind(CC , Ei)
    
    i_entropy=rbind(i_entropy,entropy=Ei)
    CER = rbind(CER , CC)
  }
  
  group[[j]]=cbind(group[[j]],entrophy=i_entropy)
  H_k=sum(as.numeric(group[[j]][,3]))
  H_G=rbind(H_G,cbind(as.character(group[[j]][1,2]),H_k))
}

ent = list(subgroup=group,H_G=H_G)

finalresult = ent$subgroup

finalresult1=c()
for(i in 1:length(finalresult))
{
  finalresult1=rbind(finalresult1,finalresult[[i]])
}
#######weight of compartment size

######below is weight of the compartment edges
comp=compartment[which(duplicated(compartment[,2])==FALSE),2]
com=c()
ECC_tissue_matrix = inf

for(i in 1:length(comp))
{
  gen_nam=finalresult1[which(comp[i]==finalresult1[,2]),1]
  com_graph=ECC_tissue_matrix[gen_nam,gen_nam]
  com_value=length(which(com_graph!=0))/length(which(ECC_tissue_matrix!=0))
  com=rbind(com,cbind(as.character(comp[i]),com_value))
}

gene=unique(finalresult1[,1])
total_gene=c()
for(i in 1:length(gene))
{
  G_v=max(as.numeric(finalresult1[which(finalresult1[,1]==gene[i]),3]))###the maximize value as the final result
  total_gene=rbind(total_gene,cbind(gene[i],G_v))
}
#total_gene = rbind(total_gene,cbind(notcom_fre[,1],notcom_fre[,3]))
total_gene1=total_gene[order(as.numeric(total_gene[,2]),decreasing=T),]




driver_gene = read.csv('NCG6_cancergenes.csv')

a1 = length(total_gene1[,1][which( total_gene1[,1][1:50] %in%driver_gene$symbol)]) 

a2 = length(total_gene1[,1][which( total_gene1[,1][1:100] %in%driver_gene$symbol)]) 

a3 = length(total_gene1[,1][which( total_gene1[,1][1:200] %in%driver_gene$symbol)]) 

res_tem = cbind(a1,a2,a3,pro[q])
res = rbind(res,res_tem)


#write.table (total_gene1, file ="lung.csv",row.names = FALSE, col.names =TRUE, quote =FALSE)
#write.table (pre, file ="Breast_pre.csv",row.names = FALSE, col.names =TRUE, quote =FALSE)
#write.table (CER, file ="cer.csv",row.names = FALSE, col.names =TRUE, quote =FALSE)

}
}
res 
