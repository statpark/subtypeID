#### codes for mRMR and MCFS using the golden standard dataset (BREAST, COLON)
#library
library(dplyr)  
library(rmcfs) #monte carlo feature selection method
library(praznik) # minimum Redundancy Maximum Relevance method
library(SNFtool)
library(PINSPlus)
library(survival)
library(ConsensusClusterPlus)
library(parallel)
library(cluster)
library(CancerSubtypes)
library(iClusterPlus)
library(diptest)
library(clusterCrit)
library(pdfCluster)
library(aricode)


#working directory= desktop
setwd("C:/Users/suzuk/Desktop")
#min-max scaling
nor_minmax = function(x){
  result = (x - min(x)) / (max(x) - min(x))
  return(result)
}


BRCA_mRNA <- read.csv("BRCA_mRNA.csv", row.names=1)
BRCA_label <- read.csv("BRCA_label.csv", sep="")
BRCA_mRNA<-as.data.frame(t(BRCA_mRNA))
BRCA_labeled<-cbind(BRCA_label,BRCA_mRNA)
BRCA_labeled_transpose <- as.data.frame(t(BRCA_labeled))
write.csv(BRCA_labeled_transpose,file="BRCA_labeled.csv",fileEncoding = 'cp949')

COAD_mRNA <- read.csv("COAD_mRNA.csv", row.names=1)
COAD_label <- read.csv("COAD_label.csv", sep="")
COAD_mRNA<-as.data.frame(t(COAD_mRNA))
COAD_labeled<-cbind(COAD_label,COAD_mRNA)
COAD_labeled_transpose <- as.data.frame(t(COAD_labeled))
write.csv(COAD_labeled_transpose,file="COAD_labeled.csv",fileEncoding = 'cp949')


##2. mRMR(minimum Redundancy Maximum Relevance Method)
#mRMR using praznik
BRCA_labeled$Label<-as.factor(BRCA_labeled$Label)
COAD_labeled$Label<-as.factor(COAD_labeled$Label)

BRCA_500<-MRMR(BRCA_mRNA, BRCA_labeled[,1], k = 500, positive = FALSE, threads = 0)
BRCA_MRMR_500<- BRCA_labeled %>% select(Label, all_of(BRCA_500$selection))

BRCA_2000<-MRMR(BRCA_mRNA, BRCA_labeled[,1], k = 2000, positive = FALSE, threads = 0)
BRCA_MRMR_2000<- BRCA_labeled %>% select(Label, all_of(BRCA_2000$selection))

COAD_500<-MRMR(COAD_mRNA, COAD_labeled[,1], k = 500, positive = FALSE, threads = 0)
COAD_MRMR_500<- COAD_labeled %>% select(Label,all_of(COAD_500$selection))

COAD_2000<-MRMR(COAD_mRNA, COAD_labeled[,1], k = 2000, positive = FALSE, threads = 0)
COAD_MRMR_2000<- COAD_labeled %>% select(Label, all_of(COAD_2000$selection))

BRCA_MRMR_500_label <- as.matrix(BRCA_MRMR_500[,1])
BRCA_MRMR_2000_label <- as.matrix(BRCA_MRMR_2000[,1])
COAD_MRMR_500_label <- as.matrix(COAD_MRMR_500[,1])
COAD_MRMR_2000_label <- as.matrix(COAD_MRMR_2000[,1])

write.csv(as.data.frame(t(COAD_MRMR_500[,-1])),file="COAD_MRMR_500.csv",fileEncoding = 'cp949')
write.csv(as.data.frame(t(COAD_MRMR_2000[,-1])),file="COAD_MRMR_2000.csv",fileEncoding = 'cp949')
write.csv(as.data.frame(t(BRCA_MRMR_500[,-1])),file="BRCA_MRMR_500.csv",fileEncoding = 'cp949')
write.csv(as.data.frame(t(BRCA_MRMR_2000[,-1])),file="BRCA_MRMR_2000.csv",fileEncoding = 'cp949')
write.csv(as.data.frame(BRCA_MRMR_500_label), file="BRCA_MRMR_500_label.csv",fileEncoding = 'cp949')
write.csv(as.data.frame(BRCA_MRMR_2000_label), file="BRCA_MRMR_2000_label.csv",fileEncoding = 'cp949')
write.csv(as.data.frame(COAD_MRMR_500_label), file="COAD_MRMR_500_label.csv",fileEncoding = 'cp949')
write.csv(as.data.frame(COAD_MRMR_2000_label), file="COAD_MRMR_2000_label.csv",fileEncoding = 'cp949')


BRCA_MRMR_500_mm<-apply(X = BRCA_MRMR_500, MARGIN = 2, FUN = "nor_minmax")
BRCA_MRMR_500_mm<-as.data.frame(BRCA_MRMR_500_mm)
BRCA_MRMR_2000_mm<-apply(X = BRCA_MRMR_2000, MARGIN = 2, FUN = "nor_minmax")
BRCA_MRMR_2000_mm<-as.data.frame(BRCA_MRMR_2000_mm)
COAD_MRMR_500_mm<-apply(X = COAD_MRMR_500, MARGIN = 2, FUN = "nor_minmax")
COAD_MRMR_500_mm<-as.data.frame(COAD_MRMR_500_mm)
COAD_MRMR_2000_mm<-apply(X = COAD_MRMR_2000, MARGIN = 2, FUN = "nor_minmax")
COAD_MRMR_2000_mm<-as.data.frame(COAD_MRMR_2000_mm)

write.csv(as.data.frame(t(COAD_MRMR_500_mm)),file="COAD_MRMR_500_mm.csv",fileEncoding = 'cp949')
write.csv(as.data.frame(t(COAD_MRMR_2000_mm)),file="COAD_MRMR_2000_mm.csv",fileEncoding = 'cp949')
write.csv(as.data.frame(t(BRCA_MRMR_500_mm)),file="BRCA_MRMR_500_mm.csv",fileEncoding = 'cp949')
write.csv(as.data.frame(t(BRCA_MRMR_2000_mm)),file="BRCA_MRMR_2000_mm.csv",fileEncoding = 'cp949')


COAD_MRMR_500 <- read.csv("COAD_MRMR_500_mm.csv", row.names=1)
COAD_MRMR_2000 <- read.csv("COAD_MRMR_2000_mm.csv", row.names=1)
BRCA_MRMR_500 <- read.csv("BRCA_MRMR_500_mm.csv", row.names=1)
BRCA_MRMR_2000 <- read.csv("BRCA_MRMR_2000_mm.csv", row.names=1)


###################################################
# algorithm runs #
run.snf <- function(omics.list, subtype.data) {
  start = Sys.time()
  omics.list = omics.list
  subtype = subtype.data$name
  alpha=0.5
  T.val=30
  num.neighbors = round(ncol(omics.list[[1]]) / 10)
  similarity.data = lapply(omics.list, function(x) {affinityMatrix(dist2(as.matrix(t(x)),as.matrix(t(x))), 
                                                                   num.neighbors, alpha)})
  if (length(similarity.data) == 1) {
    W = similarity.data[[1]]
  } else {
    W = SNF(similarity.data, num.neighbors, T.val)  
  }
  
  num.clusters = estimateNumberOfClustersGivenGraph(W, 2:MAX.NUM.CLUSTERS)[[3]]  
  clustering = spectralClustering(W, num.clusters)
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=clustering, timing=time.taken))
}

run.pins <- function(omics.list, subtype.data) {
  start = Sys.time()
  omics.list = omics.list
  subtype = subtype.data$name
  omics.transposed = lapply(omics.list, t)
  if (length(omics.list) == 1) {
    pins.ret = PINSPlus::PerturbationClustering(data=omics.transposed[[1]],
                                                kMax = MAX.NUM.CLUSTERS)
    clustering = pins.ret$cluster
    
  } else {
    pins.ret = PINSPlus::SubtypingOmicsData(dataList=omics.transposed,
                                            kMax = MAX.NUM.CLUSTERS)
    clustering = pins.ret$cluster2
  }
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=clustering, timing=time.taken))
}

ExecuteCC <- function(clusterNum,d,plot,maxK=8,clusterAlg="hc",distance="pearson",title="ConsensusClusterResult",
                      reps=500, pItem=0.8, pFeature=1,innerLinkage="average", finalLinkage="average",
                      writeTable=FALSE,weightsItem=NULL,weightsFeature=NULL, verbose=FALSE,corUse="everything") {
  start = Sys.time()
  
  if(is.list(d))
  {
    temp=NULL
    for(i in 1: length(d))
    {
      temp=rbind(temp,d[[i]])
    }
    temp=t(scale(t(temp)))
  }
  else
    temp=d
  
  originalResult=ConsensusClusterPlus(temp, maxK=maxK,clusterAlg=clusterAlg,distance=distance,title=title,
                                      reps=reps, pItem=pItem, pFeature=pFeature,plot=plot,innerLinkage=innerLinkage, finalLinkage=finalLinkage,
                                      writeTable=FALSE,weightsItem=weightsItem,weightsFeature=weightsFeature,verbose=verbose,corUse=corUse)
  clustering=originalResult[[clusterNum]][["consensusClass"]]
  distanceMatrix=originalResult[[clusterNum]][["consensusMatrix"]]
  attr(distanceMatrix,'class')="Similarity"
  time.taken = as.numeric(Sys.time() - start, units='secs')
  result=list(clustering=clustering, timing=time.taken)
  result
}


run.cc.gold <- function(plot) {
  result_cc_gold = list()
  j=1
  for (i in 1:length(SUBTYPES.DATA)) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    subtype.raw.data = read.csv(paste0(subtype,"_MCFS_500.csv"),row.names=1)
    #이름 바꿔주기csv파일 
    subtype.raw.data = as.matrix(subtype.raw.data)
    
    cur.iteration.data = subtype.raw.data
    
    if (subtype=="BRCA") {
      k=3
    } else if (subtype=="COAD") { 
      k=4
    }
    
    algorithm.ret = ExecuteCC(clusterNum=k, d=cur.iteration.data, plot=plot, maxK=8, clusterAlg="hc", distance="pearson")
    algorithm.ret$clustering.name = paste("clusteirng_", subtype, "_cc")
    
    result_cc_gold[[j]]=list(clustering.name=algorithm.ret$clustering.name, clustering=algorithm.ret$clustering, 
                             timing=algorithm.ret$timing)
    j=j+1
  } 
  save(result_cc_gold, file="result_cc_gold.RData") #이름바꾸기
}

run.nemo <- function(omics.list, subtype.data, num.clusters=NULL, is.missing.data=F, num.neighbors=NA) {
  omics.list = omics.list
  clustering = NEMO::nemo.clustering(omics.list, num.clusters, num.neighbors)
  return(list(clustering=clustering, timing=time.taken))
}

MAX.NUM.CLUSTERS=10
run.nmf <- function(omics.list, subtype.data) {
  total.time.taken = 0
  start = Sys.time()
  omics.list = omics.list
  time.taken = as.numeric(Sys.time() - start, units='secs')
  total.time.taken = total.time.taken + time.taken
  subtype = subtype.data$name
  
  for (k in 1:MAX.NUM.CLUSTERS) {
    start = Sys.time()
    file.name = paste0(subtype, '_', k)
    nmf.ret = nmf(omics.list[[1]], k, method='lee')
    coef.mat = t(coef(nmf.ret))
    time.taken = as.numeric(Sys.time() - start, units='secs')
    total.time.taken = total.time.taken + time.taken
    write.table(coef.mat, file=file.name, quote=F, row.names=F, col.names=F, sep=',')
  }
  
  explained.vars = c()
  clustering.per.num.clusters = list()
  for (k in 1:MAX.NUM.CLUSTERS) {
    file.name = paste0(subtype, '_', k)
    consensus.mat = read.csv(file.name, header=F)
    
    start = Sys.time()
    cur.clustering = apply(consensus.mat, 1, which.max)
    explained.var = sum(unlist(apply(consensus.mat, 2, var)))
    explained.vars = c(explained.vars, explained.var)
    clustering.per.num.clusters[[k]] = cur.clustering
  }
  
  dimension = get.elbow(explained.vars, is.max=F)
  nmf.clustering = clustering.per.num.clusters[[dimension]]
  return(list(clustering=nmf.clustering))  
}


ALGORITHM.NAMES = c("snf","pins","nmf","nemo")


run.benchmark <- function() {  result = list()
j=1
for (i in 1:length(SUBTYPES.DATA)) {
  current.subtype.data = SUBTYPES.DATA[[i]]
  subtype = current.subtype.data$name
  subtype.raw.data = read.csv(paste0(subtype,"_MRMR_2000_mm.csv"),row.names=1)
  #이거 csv 파일 바꿔오며 사용 
  subtype.raw.data = list(subtype.raw.data)  #리스트로 바꾸기
  
  for (algorithm.name in ALGORITHM.NAMES) {
    set.seed(42)
    print(paste('subtype', subtype, 'running algorithm', algorithm.name, 'exp'))
    
    algorithm.func.name = paste0('run.', algorithm.name)
    algorithm.func = get(algorithm.func.name)
    clustering.name = paste0('clustering_', subtype,'_', algorithm.name)
    
    cur.iteration.data = subtype.raw.data
    algorithm.ret = algorithm.func(cur.iteration.data, current.subtype.data)
    
    result[[j]]=list(clustering.name=clustering.name, clustering=algorithm.ret$clustering, timing=algorithm.ret$timing)
    j=j+1
  }
}
save(result, file="MRMR_2000.RData")  #여기 바꿔주며 이름 바꿔 저장
} 
run.benchmark()


## testing with true labels  # part1=BRCA(BREAST), part2=COLON(COAD)
truelabel<- read.table("BRCA_MRMR_2000_label.csv", quote="\"", comment.char="")
unique(truelabel$V1)
part1 = truelabel$V1
for (i in 1:length(part1)) {
  if (part1[i]=="LumA") {
    part1[i] = 1
  } else if (part1[i]=="Her2") {
    part1[i] = 2 
  } else if (part1[i]=="LumB") {
    part1[i] = 3
  } else if (part1[i]=="Normal") {
    part1[i] = 4
  } else if (part1[i]=="Basal") {
    part1[i] = 5
  }
} 
part1= as.integer(part1)

truelabel<- read.table("COAD_MRMR_2000_label.csv", quote="\"", comment.char="")
part2 = truelabel$V1
for (i in 1:length(part2)) {
  if (part2[i]=="CIN") {
    part2[i] = 1
  } else if (part2[i]=="MSI") {
    part2[i] = 2
  } else if (part2[i]=="GS") {
    part2[i] = 3
  } else if (part2[i]=="POLE") {
    part2[i] = 4
  }
}
part2 = as.integer(part2)
###################################
#### ARI and NMI results #### 
FS= c("MRMR","MCFS")

ALGORITHM.NAMES = c("snf","pins","nmf","nemo")

for (k in 1:length(FS)) {
  fs=FS[k]
  BRCA_MRMR_500 <- read.csv("COAD_MRMR_2000_mm.csv", row.names=1) #change BRCA,COAD,fs,500,2000
  data= BRCA_MRMR_500
  load("MRMR_2000.RData") #load change
  
  all.result.fs = analyze.benchmark(result) 
  all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  rownames(all.surv.pvalues) = ALGORITHM.NAMES
  colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
  all.clusterings = all.result.fs$all.clusterings
  subtype = colnames(all.surv.pvalues)[2]                         ## 1=BRCA, 2=COAD
  subtype.clusterings = all.clusterings[[subtype]]
  
  extindex = list()
  for (j in 1:length(subtype.clusterings)) {
    extindex[[j]] = list()
    
    clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
    extindex[[j]]$ARI = adj.rand.index(part2, clustering)    ## part1=BRCA, part2=COAD
    extindex[[j]]$NMI = NMI(part2, clustering)               ## part1=BRCA, part2=COAD
  }
  save(extindex, file=paste0("result/extindex2000")) 
  
}
View(extindex)
