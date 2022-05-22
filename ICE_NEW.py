# preprocess the MAF file
library(data.table)
#获取数据
load('data/new_InfluenceGraph.Rdata')
#load('data/Lung.Rdata')###prostate/breast/lung
patMutMatrix = read.csv('data/TCGA_OV.csv')
compartment=read.csv('data/Gen_compartment_together.csv')
#功能交互网络影响图

TT = c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1)

pro = c(0,0.2,0.4,0.6,0.8)
res =c()
res0 =c()
res1 =c()
res2 =c()
res3 =c()

for (t in 3:3){
for (q in 1:1){
  for(ii in 1:1){
    res_tem = c()
    
    #随机取出基因连接边
    inF=influenceGraph[intersect(colnames(patMutMatrix),row.names(influenceGraph)),intersect(colnames(patMutMatrix),row.names(influenceGraph))]
    prob=runif(dim(inF)[1], 0, 1)
    inF[which(prob <= pro[q]),which(prob < pro[q])] = 0
    inf <- inF
    rm(inF)
    inf[lower.tri(inf)] <- t(inf)[lower.tri(inf)]
    
    
    total_graph = inf
    patMatrix = patMutMatrix
    gen_fr=colSums(patMatrix[,2:length(colnames(patMutMatrix))])
    gene_name = colnames(patMutMatrix)
    
    
    #计算突变频率
    gene_frequence=gen_fr/nrow(patMatrix)
    gene_frequence=cbind(names(gene_frequence),gene_frequence)
    total_graph=(total_graph-min(total_graph))/(max(total_graph)-min(total_graph))
    tt=list(total_graph=total_graph,frequency=gene_frequence)
    
    seperate = compartment
    frequency = tt$frequency
    
    #计算不在细胞隔间和网络中的孤立点

    notcom_sep = gene_name[-(which(gene_name %in% seperate[,1])) ]
    notcom_sep = notcom_sep[2:length(notcom_sep)]
    notcom_fre = gene_frequence[which(gene_frequence[,1]%in%notcom_sep),]
    notcom_group=influenceGraph[intersect(notcom_sep,row.names(influenceGraph)),intersect(notcom_sep,row.names(influenceGraph))]
    
    N_com = length(notcom_group)
    
    
    #notcom_ei = -log(as.numeric(as.numeric(notcom_fre[,2])))*as.numeric(notcom_fre[,2])
    #notcom_fre = cbind(notcom_fre,notcom_ei)
    
    #notcom_fre[which(notcom_fre[,3]==NaN),3] = 0
    

    com=unique(seperate[,2])
    group=list()

    for(i in 1:length(com))#####divide matrix into 11 compartment subgroups将矩阵划分为11个隔室子群
    {
      com_i=seperate[which(seperate[,2]==com[i]),1]
      group_i=intersect(colnames(total_graph),com_i)
      group[[i]]=cbind(Gene=group_i,compartment=as.character(com[i]))
      
    }
    
    com_i=gene_name[-(which(gene_name %in% seperate[,1])) ]
    group_i=intersect(colnames(total_graph),com_i)
    group[[length(com)+1]]=cbind(Gene=group_i,compartment="notcompartment")
    #不在compartment的基因

    
    H_G=c()
    CER = c()
    
    for(j in 1:length(group))###calculate the entrophy for each sub-group计算每个子组的熵
    {
    
      sub_group = total_graph[group[[j]][,1],group[[j]][,1]]
      
      row.names(frequency) = frequency[,1]
      i_entropy=c()
      
      if (j == length(group)){
        N = length(rownames(total_graph[,]))
      }else{
        N = length(rownames(sub_group[,]))
      }

      ##遍历基因数
      for(k in 1:nrow(sub_group))
      {
        
        if (j == length(group)){
          n_1 =length(which(total_graph[row.names(sub_group)[k],]!=0))
          yj = frequency[colnames(total_graph)[which(total_graph[row.names(sub_group)[k],]!=0)],]
        }else{
          n_1 =length(which(sub_group[k,]!=0))
          yj = frequency[colnames(sub_group)[which(sub_group[k,]!=0 )],]
        }
      
        #n_1 =length(which(sub_group[k,]!=0))####number of the neighbor nodes邻居基因数
        
        #########################below is the cross_entropy of variation frequency
        
        #n_2 =length(which(sub_group[which(sub_group[k,]!=0),]!=0))+n_1
        
        if(is.null(dim(yj))==TRUE)
        {
          #Wi=sum(wij*log(as.numeric(frequency[row.names(sub_group)[k],2])/as.numeric(yj[2]))*as.numeric(frequency[row.names(sub_group)[k],2]))
          RH = sum(abs(log(as.numeric(frequency[row.names(sub_group)[k],2])/as.numeric(yj[2]))*as.numeric(frequency[row.names(sub_group)[k],2])))
          
        }else
        {
          #=sum(wij*log(as.numeric(frequency[row.names(sub_group)[k],2])/as.numeric(yj[,2]))*as.numeric(frequency[row.names(sub_group)[k],2]))
          RH = sum(abs(log(as.numeric(frequency[row.names(sub_group)[k],2])/as.numeric(yj[,2]))*as.numeric(frequency[row.names(sub_group)[k],2])))
          
        }
        
        if (frequency[row.names(sub_group)[k],2] == 0)
        {
          SH = 0
        }else
        {
          SH = - log(as.numeric(as.numeric(frequency[row.names(sub_group)[k],2])))*as.numeric(frequency[row.names(sub_group)[k],2])
        }
        
        CC = c()
        
        if (j!=length(group)){
          CEi = RH*TT[t] + SH}
        else
        {
          CEi = RH*TT[t] + SH
        }
        
        Ei = (1+n_1/(N-1))*CEi
        
        CC = cbind(CC , row.names(sub_group)[k],RH,SH,Ei,group[[j]][1,2])
        
        i_entropy=rbind(i_entropy,entropy=Ei)
        CER = rbind(CER , CC)
      }
      
      group[[j]]=cbind(group[[j]],entrophy=i_entropy)

    }
    
    ent = list(subgroup=group)
    
    finalresult = ent$subgroup
    
    finalresult1=c()
    for(i in 1:length(finalresult))
    {
      finalresult1=rbind(finalresult1,finalresult[[i]])
    }
    #######weight of compartment size
    
    ######below is weight of the compartment edges

    gene=unique(finalresult1[,1])
    total_gene=c()
    for(i in 1:length(gene))
    {
      G_v=max(as.numeric(finalresult1[which(finalresult1[,1]==gene[i]),3]))###the maximize value as the final result
      total_gene=rbind(total_gene,cbind(gene[i],G_v))
    }
    #total_gene = rbind(total_gene,cbind(notcom_fre[,1],notcom_fre[,3]))
    total_gene1=total_gene[order(as.numeric(total_gene[,2]),decreasing=T),]
    
    drivergene = read.csv('NCG6_cancergenes.csv')
    
    #bem = c("X2020Rule",benchmark,"CGCpointMut","CGC","CTAT", "HCD","MouseMut","oncoGene")
    
    NG0 = length(total_gene1[,1][which( total_gene1[,1] %in%drivergene$symbol)])
    
    a01 = length(total_gene1[,1][which( total_gene1[,1][1:50] %in%drivergene$symbol )])/50
    
    a02 = length(total_gene1[,1][which( total_gene1[,1][1:100] %in%drivergene$symbol)])/100
    
    a03 = length(total_gene1[,1][which( total_gene1[,1][1:200] %in%drivergene$symbol)])/200
    
    b01 = length(total_gene1[,1][which( total_gene1[,1][1:50] %in%drivergene$symbol)])/NG0
    b02 = length(total_gene1[,1][which( total_gene1[,1][1:100] %in%drivergene$symbol )])/NG0
    b03 = length(total_gene1[,1][which( total_gene1[,1][1:200] %in%drivergene$symbol )])/NG0
  
    res_tem0 = cbind(a01,b01,2*a01*b01/(a01+b01),a02,b02,2*a02*b02/(a02+b02),  a03,b03,2*a03*b03/(a03+b03),"NGC")
    
    res0 = rbind(res0,res_tem0)
  
    
     
    drivergene1 = read.csv('Census.csv')
    
    #bem = c("X2020Rule",benchmark,"CGCpointMut","CGC","CTAT", "HCD","MouseMut","oncoGene")
    
    NG1 = length(total_gene1[,1][which( total_gene1[,1] %in%drivergene1$Gene.Symbol)])
    
    a11 = length(total_gene1[,1][which( total_gene1[,1][1:50] %in%drivergene1$Gene.Symbol )])/50
    
    a12 = length(total_gene1[,1][which( total_gene1[,1][1:100] %in%drivergene1$Gene.Symbol)])/100
    
    a13 = length(total_gene1[,1][which( total_gene1[,1][1:200] %in%drivergene1$Gene.Symbol)])/200
    
    b11 = length(total_gene1[,1][which( total_gene1[,1][1:50] %in%drivergene1$Gene.Symbol)])/NG1
    b12 = length(total_gene1[,1][which( total_gene1[,1][1:100] %in%drivergene1$Gene.Symbol )])/NG1
    b13 = length(total_gene1[,1][which( total_gene1[,1][1:200] %in%drivergene1$Gene.Symbol )])/NG1
    
    res_tem1 = cbind(a11,b11,2*a11*b11/(a11+b11),a12,b12,2*a12*b12/(a12+b12),  a13,b13,2*a13*b13/(a13+b13),"CGC")
    
    res1 = rbind(res1,res_tem1)
    
    
    drivergene2 = read.csv('data.csv')
    
    #bem = c("X2020Rule",benchmark,"CGCpointMut","CGC","CTAT", "HCD","MouseMut","oncoGene")
    
    NG2 = length(total_gene1[,1][which( total_gene1[,1] %in%drivergene2$MouseMut)])
    
    a21 = length(total_gene1[,1][which( total_gene1[,1][1:50] %in%drivergene2$MouseMut )])/50
    
    a22 = length(total_gene1[,1][which( total_gene1[,1][1:100] %in%drivergene2$MouseMut)])/100
    
    a23 = length(total_gene1[,1][which( total_gene1[,1][1:200] %in%drivergene2$MouseMut)])/200
    
    b21 = length(total_gene1[,1][which( total_gene1[,1][1:50] %in%drivergene2$MouseMut)])/NG2
    b22 = length(total_gene1[,1][which( total_gene1[,1][1:100] %in%drivergene2$MouseMut )])/NG2
    b23 = length(total_gene1[,1][which( total_gene1[,1][1:200] %in%drivergene2$MouseMut )])/NG2
    
    res_tem2 = cbind(a21,b21,2*a21*b21/(a21+b21),a22,b22,2*a22*b22/(a22+b22),  a23,b23,2*a23*b23/(a23+b23),"MouseMut")
    
    res2 = rbind(res2,res_tem2)
    
    #bem = c("X2020Rule",benchmark,"CGCpointMut","CGC","CTAT", "HCD","MouseMut","oncoGene")
    
    #write.table (total_gene1, file ="ov1.csv",row.names = FALSE, col.names =TRUE, quote =FALSE)
    #write.table (pre, file ="Breast_pre.csv",row.names = FALSE, col.names =TRUE, quote =FALSE)
    #write.table (CER, file ="cer.csv",row.names = FALSE, col.names =TRUE, quote =FALSE)
    
  }
}
}
# preprocess the MAF file
library(data.table)
#获取数据
load('data/new_InfluenceGraph.Rdata')
#load('data/Lung.Rdata')###prostate/breast/lung
patMutMatrix = read.csv('data/TCGA_OV.csv')
compartment=read.csv('data/Gen_compartment_together.csv')
#功能交互网络影响图

TT = c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1)

pro = c(0,0.2,0.4,0.6,0.8)
res =c()
res0 =c()
res1 =c()
res2 =c()
res3 =c()

for (t in 3:3){
for (q in 1:1){
  for(ii in 1:1){
    res_tem = c()
    
    #随机取出基因连接边
    inF=influenceGraph[intersect(colnames(patMutMatrix),row.names(influenceGraph)),intersect(colnames(patMutMatrix),row.names(influenceGraph))]
    prob=runif(dim(inF)[1], 0, 1)
    inF[which(prob <= pro[q]),which(prob < pro[q])] = 0
    inf <- inF
    rm(inF)
    inf[lower.tri(inf)] <- t(inf)[lower.tri(inf)]
    
    
    total_graph = inf
    patMatrix = patMutMatrix
    gen_fr=colSums(patMatrix[,2:length(colnames(patMutMatrix))])
    gene_name = colnames(patMutMatrix)
    
    
    #计算突变频率
    gene_frequence=gen_fr/nrow(patMatrix)
    gene_frequence=cbind(names(gene_frequence),gene_frequence)
    total_graph=(total_graph-min(total_graph))/(max(total_graph)-min(total_graph))
    tt=list(total_graph=total_graph,frequency=gene_frequence)
    
    seperate = compartment
    frequency = tt$frequency
    
    #计算不在细胞隔间和网络中的孤立点

    notcom_sep = gene_name[-(which(gene_name %in% seperate[,1])) ]
    notcom_sep = notcom_sep[2:length(notcom_sep)]
    notcom_fre = gene_frequence[which(gene_frequence[,1]%in%notcom_sep),]
    notcom_group=influenceGraph[intersect(notcom_sep,row.names(influenceGraph)),intersect(notcom_sep,row.names(influenceGraph))]
    
    N_com = length(notcom_group)
    
    
    #notcom_ei = -log(as.numeric(as.numeric(notcom_fre[,2])))*as.numeric(notcom_fre[,2])
    #notcom_fre = cbind(notcom_fre,notcom_ei)
    
    #notcom_fre[which(notcom_fre[,3]==NaN),3] = 0
    

    com=unique(seperate[,2])
    group=list()

    for(i in 1:length(com))#####divide matrix into 11 compartment subgroups将矩阵划分为11个隔室子群
    {
      com_i=seperate[which(seperate[,2]==com[i]),1]
      group_i=intersect(colnames(total_graph),com_i)
      group[[i]]=cbind(Gene=group_i,compartment=as.character(com[i]))
      
    }
    
    com_i=gene_name[-(which(gene_name %in% seperate[,1])) ]
    group_i=intersect(colnames(total_graph),com_i)
    group[[length(com)+1]]=cbind(Gene=group_i,compartment="notcompartment")
    #不在compartment的基因

    
    H_G=c()
    CER = c()
    
    for(j in 1:length(group))###calculate the entrophy for each sub-group计算每个子组的熵
    {
    
      sub_group = total_graph[group[[j]][,1],group[[j]][,1]]
      
      row.names(frequency) = frequency[,1]
      i_entropy=c()
      
      if (j == length(group)){
        N = length(rownames(total_graph[,]))
      }else{
        N = length(rownames(sub_group[,]))
      }

      ##遍历基因数
      for(k in 1:nrow(sub_group))
      {
        
        if (j == length(group)){
          n_1 =length(which(total_graph[row.names(sub_group)[k],]!=0))
          yj = frequency[colnames(total_graph)[which(total_graph[row.names(sub_group)[k],]!=0)],]
        }else{
          n_1 =length(which(sub_group[k,]!=0))
          yj = frequency[colnames(sub_group)[which(sub_group[k,]!=0 )],]
        }
      
        #n_1 =length(which(sub_group[k,]!=0))####number of the neighbor nodes邻居基因数
        
        #########################below is the cross_entropy of variation frequency
        
        #n_2 =length(which(sub_group[which(sub_group[k,]!=0),]!=0))+n_1
        
        if(is.null(dim(yj))==TRUE)
        {
          #Wi=sum(wij*log(as.numeric(frequency[row.names(sub_group)[k],2])/as.numeric(yj[2]))*as.numeric(frequency[row.names(sub_group)[k],2]))
          RH = sum(abs(log(as.numeric(frequency[row.names(sub_group)[k],2])/as.numeric(yj[2]))*as.numeric(frequency[row.names(sub_group)[k],2])))
          
        }else
        {
          #=sum(wij*log(as.numeric(frequency[row.names(sub_group)[k],2])/as.numeric(yj[,2]))*as.numeric(frequency[row.names(sub_group)[k],2]))
          RH = sum(abs(log(as.numeric(frequency[row.names(sub_group)[k],2])/as.numeric(yj[,2]))*as.numeric(frequency[row.names(sub_group)[k],2])))
          
        }
        
        if (frequency[row.names(sub_group)[k],2] == 0)
        {
          SH = 0
        }else
        {
          SH = - log(as.numeric(as.numeric(frequency[row.names(sub_group)[k],2])))*as.numeric(frequency[row.names(sub_group)[k],2])
        }
        
        CC = c()
        
        if (j!=length(group)){
          CEi = RH*TT[t] + SH}
        else
        {
          CEi = RH*TT[t] + SH
        }
        
        Ei = (1+n_1/(N-1))*CEi
        
        CC = cbind(CC , row.names(sub_group)[k],RH,SH,Ei,group[[j]][1,2])
        
        i_entropy=rbind(i_entropy,entropy=Ei)
        CER = rbind(CER , CC)
      }
      
      group[[j]]=cbind(group[[j]],entrophy=i_entropy)

    }
    
    ent = list(subgroup=group)
    
    finalresult = ent$subgroup
    
    finalresult1=c()
    for(i in 1:length(finalresult))
    {
      finalresult1=rbind(finalresult1,finalresult[[i]])
    }
    #######weight of compartment size
    
    ######below is weight of the compartment edges

    gene=unique(finalresult1[,1])
    total_gene=c()
    for(i in 1:length(gene))
    {
      G_v=max(as.numeric(finalresult1[which(finalresult1[,1]==gene[i]),3]))###the maximize value as the final result
      total_gene=rbind(total_gene,cbind(gene[i],G_v))
    }
    #total_gene = rbind(total_gene,cbind(notcom_fre[,1],notcom_fre[,3]))
    total_gene1=total_gene[order(as.numeric(total_gene[,2]),decreasing=T),]
    
    drivergene = read.csv('NCG6_cancergenes.csv')
    
    #bem = c("X2020Rule",benchmark,"CGCpointMut","CGC","CTAT", "HCD","MouseMut","oncoGene")
    
    NG0 = length(total_gene1[,1][which( total_gene1[,1] %in%drivergene$symbol)])
    
    a01 = length(total_gene1[,1][which( total_gene1[,1][1:50] %in%drivergene$symbol )])/50
    
    a02 = length(total_gene1[,1][which( total_gene1[,1][1:100] %in%drivergene$symbol)])/100
    
    a03 = length(total_gene1[,1][which( total_gene1[,1][1:200] %in%drivergene$symbol)])/200
    
    b01 = length(total_gene1[,1][which( total_gene1[,1][1:50] %in%drivergene$symbol)])/NG0
    b02 = length(total_gene1[,1][which( total_gene1[,1][1:100] %in%drivergene$symbol )])/NG0
    b03 = length(total_gene1[,1][which( total_gene1[,1][1:200] %in%drivergene$symbol )])/NG0
  
    res_tem0 = cbind(a01,b01,2*a01*b01/(a01+b01),a02,b02,2*a02*b02/(a02+b02),  a03,b03,2*a03*b03/(a03+b03),"NGC")
    
    res0 = rbind(res0,res_tem0)
  
    
     
    drivergene1 = read.csv('Census.csv')
    
    #bem = c("X2020Rule",benchmark,"CGCpointMut","CGC","CTAT", "HCD","MouseMut","oncoGene")
    
    NG1 = length(total_gene1[,1][which( total_gene1[,1] %in%drivergene1$Gene.Symbol)])
    
    a11 = length(total_gene1[,1][which( total_gene1[,1][1:50] %in%drivergene1$Gene.Symbol )])/50
    
    a12 = length(total_gene1[,1][which( total_gene1[,1][1:100] %in%drivergene1$Gene.Symbol)])/100
    
    a13 = length(total_gene1[,1][which( total_gene1[,1][1:200] %in%drivergene1$Gene.Symbol)])/200
    
    b11 = length(total_gene1[,1][which( total_gene1[,1][1:50] %in%drivergene1$Gene.Symbol)])/NG1
    b12 = length(total_gene1[,1][which( total_gene1[,1][1:100] %in%drivergene1$Gene.Symbol )])/NG1
    b13 = length(total_gene1[,1][which( total_gene1[,1][1:200] %in%drivergene1$Gene.Symbol )])/NG1
    
    res_tem1 = cbind(a11,b11,2*a11*b11/(a11+b11),a12,b12,2*a12*b12/(a12+b12),  a13,b13,2*a13*b13/(a13+b13),"CGC")
    
    res1 = rbind(res1,res_tem1)
    
    
    drivergene2 = read.csv('data.csv')
    
    #bem = c("X2020Rule",benchmark,"CGCpointMut","CGC","CTAT", "HCD","MouseMut","oncoGene")
    
    NG2 = length(total_gene1[,1][which( total_gene1[,1] %in%drivergene2$MouseMut)])
    
    a21 = length(total_gene1[,1][which( total_gene1[,1][1:50] %in%drivergene2$MouseMut )])/50
    
    a22 = length(total_gene1[,1][which( total_gene1[,1][1:100] %in%drivergene2$MouseMut)])/100
    
    a23 = length(total_gene1[,1][which( total_gene1[,1][1:200] %in%drivergene2$MouseMut)])/200
    
    b21 = length(total_gene1[,1][which( total_gene1[,1][1:50] %in%drivergene2$MouseMut)])/NG2
    b22 = length(total_gene1[,1][which( total_gene1[,1][1:100] %in%drivergene2$MouseMut )])/NG2
    b23 = length(total_gene1[,1][which( total_gene1[,1][1:200] %in%drivergene2$MouseMut )])/NG2
    
    res_tem2 = cbind(a21,b21,2*a21*b21/(a21+b21),a22,b22,2*a22*b22/(a22+b22),  a23,b23,2*a23*b23/(a23+b23),"MouseMut")
    
    res2 = rbind(res2,res_tem2)
    
    #bem = c("X2020Rule",benchmark,"CGCpointMut","CGC","CTAT", "HCD","MouseMut","oncoGene")
    
    #write.table (total_gene1, file ="ov1.csv",row.names = FALSE, col.names =TRUE, quote =FALSE)
    #write.table (pre, file ="Breast_pre.csv",row.names = FALSE, col.names =TRUE, quote =FALSE)
    #write.table (CER, file ="cer.csv",row.names = FALSE, col.names =TRUE, quote =FALSE)
    
  }
}
}
