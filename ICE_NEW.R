library(data.table)
load('data/new_InfluenceGraph.Rdata')
#load('data/Lung.Rdata')###prostate/breast/lung
patMutMatrix = read.csv('data/CancerPan.csv',sep=' ')
compartment=read.csv('data/Gen_compartment_together.csv')

W = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)#Cross-entropy weights
pro = c(0,0.2,0.4,0.6,0.8)#Probability of preserving gene-linked edges

for (t in 1:9){
  for (q in 1:5){

    #Disconnecting gene-to-gene linkage edges with probability p
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
    
    #Calculating mutation frequency
    gene_frequence=gen_fr/nrow(patMatrix)
    gene_frequence=cbind(names(gene_frequence),gene_frequence)
    total_graph=(total_graph-min(total_graph))/(max(total_graph)-min(total_graph))
    tt=list(total_graph=total_graph,frequency=gene_frequence)
    
    seperate = compartment
    frequency = tt$frequency

    com=unique(seperate[,2])
    group=list()

    for(i in 1:length(com))#divide matrix into 11 compartment subgroups
    {
      com_i=seperate[which(seperate[,2]==com[i]),1]
      group_i=intersect(colnames(total_graph),com_i)
      group[[i]]=cbind(Gene=group_i,compartment=as.character(com[i]))
    }
    
    #Genes not in the cell compartment
    com_i=gene_name[-(which(gene_name %in% seperate[,1])) ]
    group_i=intersect(colnames(total_graph),com_i)
    group[[length(com)+1]]=cbind(Gene=group_i,compartment="notcompartment")
  
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

      for(k in 1:nrow(sub_group))
      {
        
        if (j == length(group)){
          n_1 =length(which(total_graph[row.names(sub_group)[k],]!=0))#number of the neighbor nodes
          yj = frequency[colnames(total_graph)[which(total_graph[row.names(sub_group)[k],]!=0)],]
        }else{
          n_1 =length(which(sub_group[k,]!=0))
          yj = frequency[colnames(sub_group)[which(sub_group[k,]!=0 )],]
        }

        #########################below is the cross_entropy of variation frequency
        
        if(is.null(dim(yj))==TRUE)
        {
          RH = sum(abs(log(as.numeric(frequency[row.names(sub_group)[k],2])/as.numeric(yj[2]))*as.numeric(frequency[row.names(sub_group)[k],2])))
          
        }else
        {
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
          CEi = RH*W[t] +  SH}
        else
        {
          CEi = RH*W[t] + SH
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
    gene=unique(finalresult1[,1])
    total_gene=c()
    for(i in 1:length(gene))
    {
      G_v=max(as.numeric(finalresult1[which(finalresult1[,1]==gene[i]),3]))###the maximize value as the final result
      total_gene=rbind(total_gene,cbind(gene[i],G_v))
    }
    total_gene1=total_gene[order(as.numeric(total_gene[,2]),decreasing=T),]
  }
}





