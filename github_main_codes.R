library(SNFtool)
library(PINSPlus)
library(survival)
library(ConsensusClusterPlus)
library(parallel)
library(cluster)
library(CancerSubtypes)
library(iClusterPlus)
library(diptest)

ALGORITHM.NAMES = c('cc', 'nmf', 'pins', 'icluster', 'snf', 'nemo')
ALGORITHM.DISPLAY.NAMES = as.list(c('CC', 'CNMF', 'PINS', 'iClusterBayes', 'SNF', 'NEMO'))

OMIC.SUBSETS = list('exp')
names(OMIC.SUBSETS) = c('exp')

######## FUNCTIONS RUN ############
get.raw.data <- function(subtype.name,datasets.path = get.dataset.dir.path(), only.primary=NA, intersect.patients=T) {
  omics.dir = file.path(datasets.path, subtype.name)
  omics.files = list.files(omics.dir)
  omics.files = setdiff(omics.files, c('methy', 'mirna', 'survival'))  # list ([1] exp, [2] methy, [3] mirna)로 만든다
  raw.data = lapply(file.path(omics.dir, omics.files), read.table) # 위 리스트 순서대로 데이터 불러옴 
  
  if (!is.na(only.primary)) {
    raw.data = lapply(raw.data, function(x) filter.non.tumor.samples(x, only.primary = only.primary))
  }
  if (intersect.patients) {
    name.corrected.data = fix.patient.names(raw.data)
    patients.intersection = Reduce(intersect, lapply(name.corrected.data, colnames))
    ret.data = lapply(name.corrected.data, function(datum) datum[,patients.intersection])  
  } else {
    ret.data = raw.data
  }
  return(ret.data)
}
get.dataset.dir.path <- function() {
  return('/Users/kyle/Downloads/data')
}
set.omics.list.attr <- function(subtype.raw.data, subtype.data) {
  attr(subtype.raw.data[[1]], 'is.seq') = subtype.data$is.rna.seq
  return(subtype.raw.data)
}
filter.non.tumor.samples <- function(raw.datum, only.primary=only.primary) {
  # 01 is primary, 06 is metastatic, 03 is blood derived cancer
  if (!only.primary)
    return(raw.datum[,substring(colnames(raw.datum), 14, 15) %in% c('01', '03', '06')])
  else
    return(raw.datum[,substring(colnames(raw.datum), 14, 15) %in% c('01')])
}
get.fixed.names <- function(patient.names, include.type=F) {
  # fix the TCGA names to only include the patient ids
  if (include.type) {
    return(gsub('-', '\\.', toupper(substring(patient.names, 1, 15))))
  } else {
    return(gsub('-', '\\.', toupper(substring(patient.names, 1, 12))))  
  }
}
fix.patient.names <- function(subtype.raw.data, include.type=F) {
  for (i in 1:length(subtype.raw.data)) {
    colnames(subtype.raw.data[[i]]) = get.fixed.names(colnames(subtype.raw.data[[i]]),
                                                      include.type)
  }
  return(subtype.raw.data)
}
get.elbow <- function(values, is.max) {
  second.derivatives = c()
  for (i in 2:(length(values) - 1)) {
    second.derivative = values[i + 1] + values[i - 1] - 2 * values[i]
    second.derivatives = c(second.derivatives, second.derivative)
  }
  print(second.derivatives)
  if (is.max) {
    return(which.max(second.derivatives) + 1)
  } else {
    return(which.min(second.derivatives) + 1)
  }
}
get.clustering.silhouette <- function(raw.data, clustering) {
  sils = c()
  for (i in 1:length(raw.data)) {
    x = raw.data[[i]]
    distmatrix = dist2(as.matrix(t(x)),as.matrix(t(x)))
    sil = silhouette(clustering, dmatrix = distmatrix)[,3]
    sils = c(sils, mean(sil))
  }
  return(mean(sils))
}
get.subtype.survival.path <- function(subtype) {
  datasets.path = get.dataset.dir.path()
  survival.file.path = file.path(datasets.path, subtype, 'survival')
  return(survival.file.path)
}
get.clinical.params.dir <- function() {
  return('/Users/kyle/Downloads/data/clinical/')
}
get.clinical.params <- function(subtype.name) {
  clinical.data.path = paste(get.clinical.params.dir(), subtype.name, sep = '')
  clinical.params = read.table(clinical.data.path,
                               sep='\t', header=T, row.names = 1, stringsAsFactors = F)
  rownames.with.duplicates = get.fixed.names(rownames(clinical.params))  
  clinical.params = clinical.params[!duplicated(rownames.with.duplicates),]
  rownames(clinical.params) = rownames.with.duplicates[!duplicated(rownames.with.duplicates)]
  return(clinical.params)
}
get.dataset.dir.path <- function() {
  return('/Users/kyle/Downloads/data')
}
get.subtype.survival.path <- function(subtype) {
  datasets.path = get.dataset.dir.path()
  survival.file.path = file.path(datasets.path, subtype, 'survival')
  return(survival.file.path)
}
get.clustering.results.dir.path <- function() {
  return('/Users/kyle/Desktop/result_final')
}
get.plots.dir.path <- function() {
  results.dir.path = get.clustering.results.dir.path()
  return(file.path(results.dir.path, 'plots'))
}
get.tables.dir.path <- function() {
  results.dir.path = get.clustering.results.dir.path()
  return(file.path(results.dir.path, 'tables'))
}
subtype.to.display.name <- function(subtype) {
  for (i in 1:length(SUBTYPES.DATA)) {
    if (SUBTYPES.DATA[[i]]$name == subtype) {
      return(SUBTYPES.DATA[[i]]$display.name)
    }
  }
}
print.matrix.latex.format <- function(mat) {
  print(do.call(paste, as.list(c(colnames(mat), sep=' & '))))
  for (i in 1:nrow(mat)) {
    print(do.call(paste, as.list(c(rownames(mat)[i], round(mat[i,], digits=2), sep=' & '))))
  }
}
log.and.normalize <- function(omics.data, subtype.data, normalize=T, filter.var=F) {
  # filter features with no variance at all
  for (i in 1:length(omics.data)) {
    omics.data[[i]] = omics.data[[i]][apply(omics.data[[i]], 1, var) > 0,]
  }
  
  for (i in 1:length(omics.data)) {
    if (attr(omics.data[[i]], 'is.seq')) {
      omics.data[[i]] = log(1+omics.data[[i]])
    }
  }
  
  if (filter.var) {
    omics.data = lapply(omics.data, keep.high.var.features)
  }
  
  if (normalize) {
    omics.data = lapply(omics.data, normalize.matrix)    
  }
  
  return(omics.data)
}
normalize.matrix <- function(data.matrix) {
  temp = data.matrix - rowMeans(data.matrix)
  should.keep = (apply(temp, 1, sd) != 0)
  return ((temp / apply(temp, 1, sd))[should.keep, ])
}
MAX.NUM.CLUSTERS = 10
# evaluation
analyze.benchmark <- function(result_data) {
  all.clusterings = list()
  all.timings = list()
  j=1
  for (i in 1:length(SUBTYPES.DATA)) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    #subtype.raw.data = get.raw.data(subtype,only.primary=current.subtype.data$only.primary)
    
    all.clusterings[[subtype]] = list()
    all.timings[[subtype]] = list()
    
    for (algorithm.name in ALGORITHM.NAMES) {
      all.clusterings[[subtype]][[algorithm.name]] = list()
      all.timings[[subtype]][[algorithm.name]] = list()
      
      clustering = result_data[[j]]$clustering
      timing = result_data[[j]]$timing
      
      #if (!any(is.na(clustering))) {
      #  names(clustering) = colnames(subtype.raw.data[[1]])
      #}
      
      all.clusterings[[subtype]][[algorithm.name]] = clustering
      all.timings[[subtype]][[algorithm.name]] = timing
      j=j+1
    }
  }
  return(list(all.clusterings=all.clusterings, all.timings=all.timings))
}
check.survival <- function(groups, subtype, survival.file.path) {
  if (missing(survival.file.path)) {
    survival.file.path = get.subtype.survival.path(subtype)
  }
  survival.data = read.table(survival.file.path, header = TRUE)
  patient.names = names(groups)
  patient.names.in.file = as.character(survival.data[, 1])
  patient.names.in.file = toupper(gsub('-', '\\.', substring(patient.names.in.file, 1, 12)))
  
  stopifnot(all(patient.names %in% patient.names.in.file))
  
  indices = match(patient.names, patient.names.in.file)
  ordered.survival.data = survival.data[indices,]
  ordered.survival.data["cluster"] <- groups
  ordered.survival.data$Survival[is.na(ordered.survival.data$Survival)] = 0
  ordered.survival.data$Death[is.na(ordered.survival.data$Death)] = 0
  return(survdiff(Surv(Survival, Death) ~ cluster, data=ordered.survival.data))
}
get.empirical.surv <- function(clustering, subtype) {
  set.seed(42)
  surv.ret = check.survival(clustering, subtype)
  orig.chisq = surv.ret$chisq
  orig.pvalue = get.logrank.pvalue(surv.ret)
  # The initial number of permutations to run
  num.perms = round(min(max(10 / orig.pvalue, 1000), 1e6))
  should.continue = T
  
  total.num.perms = 0
  total.num.extreme.chisq = 0
  
  while (should.continue) {
    print(num.perms)
    perm.chisq = as.numeric(mclapply(1:num.perms, function(i) {
      cur.clustering = sample(clustering)
      names(cur.clustering) = names(clustering)
      cur.chisq = check.survival(cur.clustering, subtype)$chisq
      return(cur.chisq)
    }, mc.cores=60))
    
    total.num.perms = total.num.perms + num.perms
    total.num.extreme.chisq = total.num.extreme.chisq + sum(perm.chisq >= orig.chisq)
    
    binom.ret = binom.test(total.num.extreme.chisq, total.num.perms)
    cur.pvalue = binom.ret$estimate
    cur.conf.int = binom.ret$conf.int
    
    print(c(total.num.extreme.chisq, total.num.perms))
    print(cur.pvalue)
    print(cur.conf.int)
    
    sig.threshold = 0.05
    is.conf.small = ((cur.conf.int[2] - cur.pvalue) < min(cur.pvalue / 10, 0.01)) & ((cur.pvalue - cur.conf.int[1]) < min(cur.pvalue / 10, 0.01))
    is.threshold.in.conf = cur.conf.int[1] < sig.threshold & cur.conf.int[2] > sig.threshold
    if ((is.conf.small & !is.threshold.in.conf) | (total.num.perms > 2e7)) {
      should.continue = F
    } else {
      num.perms = 1e5
    }
  }
  
  return(list(pvalue = cur.pvalue, conf.int = cur.conf.int, total.num.perms=total.num.perms, 
              total.num.extreme.chisq=total.num.extreme.chisq))
}
get.logrank.pvalue <- function(survdiff.res) {
  1 - pchisq(survdiff.res$chisq, length(survdiff.res$n) - 1)  
}
benchmark.omics.time <- function(benchmark.results) {
  all.alg.times = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  rownames(all.alg.times) = ALGORITHM.NAMES
  colnames(all.alg.times) = sapply(SUBTYPES.DATA, function(x) x$name)
  all.timings = benchmark.results$all.timings
  for (i in 1:length(all.timings)) {
    subtype = colnames(all.alg.times)[i]
    subtype.timings = all.timings[[subtype]]
    for (j in 1:length(subtype.timings)) {
      timing = subtype.timings[[ALGORITHM.NAMES[j]]]
      all.alg.times[j, i] = timing
    }
  }
  return(all.alg.times)
}
benchmark.omics.num.clusters <- function(benchmark.results) {
  num.clusters = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  rownames(num.clusters) = ALGORITHM.NAMES
  colnames(num.clusters) = sapply(SUBTYPES.DATA, function(x) x$name)
  all.clusterings = benchmark.results$all.clusterings
  for (i in 1:length(all.clusterings)) {
    subtype = colnames(num.clusters)[i]
    subtype.clusterings = all.clusterings[[subtype]]
    for (j in 1:length(subtype.clusterings)) {
      clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
      num.clusters[j, i] = max(clustering)
    }
  }
  return(num.clusters)
}
benchmark.omics.surv <- function(benchmark.results) {
  all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  rownames(all.surv.pvalues) = ALGORITHM.NAMES
  colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
  all.clusterings = benchmark.results$all.clusterings
  for (i in 1:length(all.clusterings)) {
    subtype = colnames(all.surv.pvalues)[i]
    subtype.clusterings = all.clusterings[[subtype]]
    for (j in 1:length(subtype.clusterings)) {
      surv.path = file.path("/Users/kyle/Desktop/result_final/surv", paste(subtype, ALGORITHM.NAMES[j], 'surv', sep='_'))
      
      if (file.exists(surv.path)) {
        load(surv.path)
        pvalue = empirical.surv.ret$pvalue
      } else {
        clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
        if (length(table(clustering)) > 1) {
          pvalue = -log10(get.logrank.pvalue(check.survival(clustering, subtype)))
          empirical.surv.ret = get.empirical.surv(clustering, subtype)
          save(empirical.surv.ret, file=surv.path)
          pvalue = empirical.surv.ret$pvalue
        } else {
          pvalue = NA
        }  
      }
      all.surv.pvalues[j, i] = pvalue
    }
  }
  return(all.surv.pvalues)
}
perform.all.analyses <- function(benchmark.ret) {
  
  for (i in 1:3) {
    cur.func = list(benchmark.omics.time, benchmark.omics.num.clusters,benchmark.omics.surv)[[i]]
    
    benchmark.data = cur.func(benchmark.ret)
    
    displayed.benchmark.data = benchmark.data      
    colnames(displayed.benchmark.data)[1:ncol(displayed.benchmark.data)] = 
      sapply(as.list(colnames(displayed.benchmark.data)[1:ncol(displayed.benchmark.data)]), 
             subtype.to.display.name)
    rownames(displayed.benchmark.data) = ALGORITHM.DISPLAY.NAMES
    print.matrix.latex.format(displayed.benchmark.data)
    
    table.name = c('runtime', 'num_cluster', 'survival')[i]
    
    file.name = deparse(substitute(benchmark.ret))
    fs.name = substr(file.name, 12, 14)
    if (file.name=="all.result") {
      file=file.path("/Users/kyle/Desktop/result/tables", paste0(table.name, '.csv'))
    } else file = file.path("/Users/kyle/Desktop/result_fs", paste0(fs.name,'/tables/',table.name, '.csv'))
    
    write.csv(displayed.benchmark.data, file=file)
    
    print('------------------------')
  }
}

# FS methods
FSbyDip<-function(data, dataname, cut.type="topk", value) {
  dip.pval=matrix(0,nrow(data))
  dip = apply(data,1,dip.test)
  for (i in 1:length(dip)) dip.pval[[i]] = dip[[i]]$p.value
  rownames(dip.pval) = names(dip)
  if(cut.type=="topk")
  {
    index= sort(dip.pval, decreasing = FALSE, index.return=TRUE)
    if(value>nrow(data))
    {
      value=nrow(data)
      cat("Warning: the feature selection number is beyond the original feature numnber")
    }
    cutoff=index$x[value]
    index=index$ix[1:value]
    selectData=data[index,]
  }
  selectData
}

# ALGORITHMS 
##CC
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
run.cc <- function(plot) {
  result_cc = list()
  j=1
  for (i in 1:length(SUBTYPES.DATA)) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    subtype.raw.data = get.raw.data(subtype, only.primary=current.subtype.data$only.primary)
    subtype.raw.data = set.omics.list.attr(subtype.raw.data, current.subtype.data)
    subtype.raw.data = log.and.normalize(subtype.raw.data, subtype.data)
    subtype.raw.data = as.matrix(subtype.raw.data[[1]])
    
    cur.iteration.data = subtype.raw.data
    wd = paste0("/Users/kyle/Desktop/result_final/",subtype)
    setwd(wd)
    
    if (subtype=="aml") {
      k=5
    } else if (subtype=="breast") { 
      k=3
    } else if (subtype=="colon") {
      k=4
    } else if (subtype=="gbm") {
      k=3
    }
    
    algorithm.ret = ExecuteCC(clusterNum=k, d=cur.iteration.data, plot=plot, maxK=8, clusterAlg="hc", distance="pearson")
    algorithm.ret$clustering.name = paste0("clusteirng_", subtype, "_cc")
    
    result_cc[[j]]=list(clustering.name=algorithm.ret$clustering.name, clustering=algorithm.ret$clustering, 
                        timing=algorithm.ret$timing)
    j=j+1
  } 
  setwd("/Users/kyle/Desktop/result_final")
  save(result_cc, file="result_cc.RData")
}
tune.icluster <- function() {
  for (i in 1:length(SUBTYPES.DATA)) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    subtype.raw.data = get.raw.data(subtype, only.primary=current.subtype.data$only.primary)
    subtype.raw.data = set.omics.list.attr(subtype.raw.data, current.subtype.data)
    subtype.raw.data = log.and.normalize(subtype.raw.data, subtype.data)
    subtype.raw.data = t(as.matrix(subtype.raw.data[[1]]))
    cur.iteration.data = subtype.raw.data
    
    tune = iClusterPlus::tune.iClusterBayes(cpus=8, dt1=cur.iteration.data, K=1:8)
    tune_wd = paste0("/Users/kyle/Desktop/result_final/",subtype,"/icluster/tune")
    save(tune, file=tune_wd)
    
    mat = matrix(0, nrow=8, ncol=2)
    for (i in 1:length(tune$fit)) {
      mat[i,1] = tune$fit[[i]]$BIC
      mat[i,2] = tune$fit[[i]]$dev.ratio
    }
    mat_wd = paste0("/Users/kyle/Desktop/result_final/",subtype,"/icluster/mat")
    save(mat, file=mat_wd)
  }
}
run.icluster <- function() {
  result_icluster = list()
  j=1
  for (i in 1:length(SUBTYPES.DATA)) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    subtype.raw.data = get.raw.data(subtype, only.primary=current.subtype.data$only.primary)
    subtype.raw.data = set.omics.list.attr(subtype.raw.data, current.subtype.data)
    subtype.raw.data = log.and.normalize(subtype.raw.data, subtype.data)
    subtype.raw.data = t(as.matrix(subtype.raw.data[[1]]))
    cur.iteration.data = subtype.raw.data
    
    if (subtype=="aml") {
      K=4
    } else if (subtype=="liver") { 
      K=4
    } else if (subtype=="gbm") {
      K=3
    }
    
    start = Sys.time()
    result = iClusterPlus::iClusterBayes(dt1=cur.iteration.data, K=K)
    result$clustering.name = paste("clusteirng_", subtype, "_icluster")
    time.taken = as.numeric(Sys.time() - start, units='secs')
    result_icluster[[j]]=list(clustering.name=result$clustering.name, clustering=result$clusters, 
                              timing=time.taken)
    j=j+1
  } 
  setwd("/Users/kyle/Desktop/result_final")
  save(result_icluster, file="result_icluster.RData")
}
run.icluster.unnorm <- function() {
  result_icluster_unnorm = list()
  j=1
  for (i in 1:length(SUBTYPES.DATA)) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    subtype.raw.data = get.raw.data(subtype, only.primary=current.subtype.data$only.primary)
    subtype.raw.data = set.omics.list.attr(subtype.raw.data, current.subtype.data)
    subtype.raw.data = log.and.normalize(subtype.raw.data, subtype.data, normalize=FALSE)
    subtype.raw.data = t(as.matrix(subtype.raw.data[[1]]))
    cur.iteration.data = subtype.raw.data
    
    if (subtype=="aml") {
      K=1
    } else if (subtype=="liver") { 
      K=2
    } else if (subtype=="gbm") {
      K=1
    }
    
    start = Sys.time()
    result = iClusterPlus::iClusterBayes(dt1=cur.iteration.data, K=K)
    result$clustering.name = paste("clusteirng_", subtype, "_icluster")
    time.taken = as.numeric(Sys.time() - start, units='secs')
    result_icluster_unnorm[[j]]=list(clustering.name=result$clustering.name, clustering=result$clusters, 
                                     timing=time.taken)
    j=j+1
  } 
  setwd("/Users/kyle/Desktop/result_final")
  save(result_icluster_unnorm, file="result_icluster_unnorm.RData")
}
run.nmf <- function(omics.list, subtype.data) {
  total.time.taken = 0
  start = Sys.time()
  omics.list = log.and.normalize(omics.list, subtype.data, normalize=F)
  time.taken = as.numeric(Sys.time() - start, units='secs')
  total.time.taken = total.time.taken + time.taken
  subtype = subtype.data$name
  
  for (k in 1:MAX.NUM.CLUSTERS) {
    start = Sys.time()
    file.name = paste0('/Users/kyle/Desktop/result_final/nmf/', subtype, '_', k)
    nmf.ret = nmf(omics.list[[1]], k, method='lee')
    coef.mat = t(coef(nmf.ret))
    time.taken = as.numeric(Sys.time() - start, units='secs')
    total.time.taken = total.time.taken + time.taken
    write.table(coef.mat, file=file.name, quote=F, row.names=F, col.names=F, sep=',')
  }
  
  explained.vars = c()
  clustering.per.num.clusters = list()
  for (k in 1:MAX.NUM.CLUSTERS) {
    file.name = paste0('/Users/kyle/Desktop/result_final/nmf/', subtype, '_', k)
    consensus.mat = read.csv(file.name, header=F)
    
    start = Sys.time()
    cur.clustering = apply(consensus.mat, 1, which.max)
    explained.var = sum(unlist(apply(consensus.mat, 2, var)))
    explained.vars = c(explained.vars, explained.var)
    clustering.per.num.clusters[[k]] = cur.clustering
    time.taken = as.numeric(Sys.time() - start, units='secs')
    total.time.taken = total.time.taken + time.taken
  }
  
  dimension = get.elbow(explained.vars, is.max=F)
  nmf.clustering = clustering.per.num.clusters[[dimension]]
  return(list(clustering=nmf.clustering, timing=total.time.taken))  
}
run.pins <- function(omics.list, subtype.data) {
  start = Sys.time()
  omics.list = log.and.normalize(omics.list, subtype.data)
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
run.snf <- function(omics.list, subtype.data) {
  start = Sys.time()
  omics.list = log.and.normalize(omics.list, subtype.data)
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
run.nemo <- function(omics.list, subtype.data, num.clusters=NULL, is.missing.data=F, num.neighbors=NA) {
  start = Sys.time()
  omics.list = log.and.normalize(omics.list, subtype.data)
  clustering = NEMO::nemo.clustering(omics.list, num.clusters, num.neighbors)
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=clustering, timing=time.taken))
}
run.benchmark.sn <- function() {
  result_sn = list()
  j=1
  for (i in 1:length(SUBTYPES.DATA)) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    subtype.raw.data = read.csv(paste0("/Users/kyle/Desktop/",subtype,"_MRMR_500.csv"),row.names=1)
    subtype.raw.data = set.omics.list.attr(subtype.raw.data, current.subtype.data)
    
    for (algorithm.name in ALGORITHM.NAMES) {
      set.seed(42)
      print(paste('subtype', subtype, 'running algorithm', algorithm.name, 'exp'))
      
      algorithm.func.name = paste0('run.', algorithm.name)
      algorithm.func = get(algorithm.func.name)
      clustering.name = paste0('clustering_', subtype,'_', algorithm.name)
      
      cur.iteration.data = subtype.raw.data
      algorithm.ret = algorithm.func(cur.iteration.data, current.subtype.data)
      
      result_sn[[j]]=list(clustering.name=clustering.name, clustering=algorithm.ret$clustering, timing=algorithm.ret$timing)
      j=j+1
    }
  }
  save(result_sn, file="/Users/kyle/Desktop/result_final/result_sn.RData")
} # for SNF, NEMO
run.benchmark.np <- function() {
  result_np = list()
  j=1
  for (i in 1:length(SUBTYPES.DATA)) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    subtype.raw.data = get.raw.data(subtype, only.primary=current.subtype.data$only.primary)
    subtype.raw.data = set.omics.list.attr(subtype.raw.data, current.subtype.data)
    
    for (algorithm.name in ALGORITHM.NAMES) {
      set.seed(42)
      print(paste('subtype', subtype, 'running algorithm', algorithm.name, 'exp'))
      
      algorithm.func.name = paste0('run.', algorithm.name)
      algorithm.func = get(algorithm.func.name)
      clustering.name = paste0('clustering_', subtype,'_', algorithm.name)
      
      cur.iteration.data = subtype.raw.data
      algorithm.ret = algorithm.func(cur.iteration.data, current.subtype.data)
      
      result_np[[j]]=list(clustering.name=clustering.name, clustering=algorithm.ret$clustering, timing=algorithm.ret$timing)
      j=j+1
    }
  }
  save(result_np, file="/Users/kyle/Desktop/result_final/result_np.RData")
} # for NMF, PINS

# ALGORITHMS W/O log.and.normalize (for FS)
run.nmf.fs <- function(omics.list, subtype.data, fs) {
  total.time.taken = 0
  start = Sys.time()
  time.taken = as.numeric(Sys.time() - start, units='secs')
  total.time.taken = total.time.taken + time.taken
  
  #save.subtype.matlab.format(omics.list)
  subtype = subtype.data$name
  
  for (k in 1:MAX.NUM.CLUSTERS) {
    start = Sys.time()
    file.name = paste0('/Users/kyle/Desktop/result_fs/', fs, '/nmf/NMF', k, '_', subtype)
    nmf.ret = nmf(omics.list[[1]], k, method='lee')
    coef.mat = t(coef(nmf.ret))
    time.taken = as.numeric(Sys.time() - start, units='secs')
    total.time.taken = total.time.taken + time.taken
    write.table(coef.mat, file=file.name, quote=F, row.names=F, col.names=F, sep=',')
  }
  
  explained.vars = c()
  clustering.per.num.clusters = list()
  for (k in 1:MAX.NUM.CLUSTERS) {
    file.name = paste0('/Users/kyle/Desktop/result_fs/', fs, '/nmf/NMF', k, '_', subtype)
    consensus.mat = read.csv(file.name, header=F)
    
    start = Sys.time()
    cur.clustering = apply(consensus.mat, 1, which.max)
    explained.var = sum(unlist(apply(consensus.mat, 2, var)))
    explained.vars = c(explained.vars, explained.var)
    clustering.per.num.clusters[[k]] = cur.clustering
    time.taken = as.numeric(Sys.time() - start, units='secs')
    total.time.taken = total.time.taken + time.taken
  }
  
  dimension = get.elbow(explained.vars, is.max=F)
  nmf.clustering = clustering.per.num.clusters[[dimension]]
  return(list(clustering=nmf.clustering, timing=total.time.taken))  
} # normalize = F
run.pins.fs <- function(omics.list, subtype.data) {
  start = Sys.time()
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
} # normalize = F
run.snf.fs <- function(omics.list, subtype.data) {
  start = Sys.time()
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
run.nemo.fs <- function(omics.list, subtype.data, num.clusters=NULL, is.missing.data=F, num.neighbors=NA) {
  start = Sys.time()
  clustering = NEMO::nemo.clustering(omics.list, num.clusters, num.neighbors)
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=clustering, timing=time.taken))
}
# CC, SNF, NEMO 는 log and normalization 했음

##################################### CC & iClusterBayes ##########################################################
############## CC w/o FS
#run.cc(plot="png") => aml: 2 or 5, gbm: 3, lihc : 3
run.cc(plot="NULL")
load("/Users/kyle/Desktop/result_final/result_cc.RData")
View(result_cc)

############## iCluster w/o FS (기준 : dev.ratio 로 고르고 bic 도 종합적으로 봄)
tune.icluster() 
#<UNNORMORMALIZED> aml: 2(K=1), lihc : 3(K=2), gbm: 2(K=1) <NORMALIZED> aml: 2(K=1), lihc : 3(K=2), gbm: 2(K=1)

run.icluster.unnorm()
run.icluster()

# AML 
load("/Users/kyle/Desktop/result_final/aml/icluster/mat")
par(mfrow=c(1,2))
plot(mat[1:8,1], type="b", ylab="BIC", xlab="K for AML")
plot(mat[1:8,2], type="b", ylab="deviance ratio", xlab="K for AML")
# <UNNORM> K = 1 (2 subtypes) <NORM> K = 4 (5 subtypes)

# LIHC
load("/Users/kyle/Desktop/result_final/liver/icluster/mat")
par(mfrow=c(1,2))
plot(mat[1:6,1], type="b", ylab="BIC", xlab="K for LIHC")
plot(mat[1:6,2], type="b", ylab="deviance ratio", xlab="K for LIHC")
# <UNNORM> k = 2 (3 subtypes) <NORM> K = 4 (5 subtypes) (bic기준 K=2)

# GBM
load("/Users/kyle/Desktop/result_final/gbm/icluster/mat")
par(mfrow=c(1,2))
plot(mat[1:6,1], type="b", ylab="BIC", xlab="K for GBM")
plot(mat[1:6,2], type="b", ylab="deviance ratio", xlab="K for GBM")
# <UNNORM> k = 1 (2 subtypes) <NORM> K = 3 (4 subtypes)

#iCluster patient attribute
load("/Users/kyle/Desktop/result_final/result_icluster.RData")
load("/Users/kyle/Desktop/result_final/result_icluster_unnorm.RData")
for (i in 1:3) {
  attr(result_icluster[[i]]$clustering, "names") = names(result_cc[[i]]$clustering)
}

##################################### ALGORITHMS WITH FS #########################################
SUBTYPES.DATA = list(
  list(name='aml', only.primary=F, is.rna.seq=T, is.mirna.seq=T, display.name='AML'),
  list(name='gbm', only.primary=T, is.rna.seq=F, is.mirna.seq=F, display.name='GBM'),
  list(name='breast', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='BIC'),
  list(name='colon', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='COAD'))

FS=c("var","mad","med","dip")
########### CC with FS
run.cc.var <- function(plot) {
  result_cc = list()
  j=1
  for (i in 1:length(SUBTYPES.DATA)) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    subtype.raw.data = get.raw.data(subtype, only.primary=current.subtype.data$only.primary)
    subtype.raw.data = set.omics.list.attr(subtype.raw.data, current.subtype.data)
    subtype.raw.data = log.and.normalize(subtype.raw.data, subtype.data)
    subtype.raw.data = as.matrix(subtype.raw.data[[1]])
    subtype.raw.data = FSbyVar(subtype.raw.data, cut.type="topk", 2000)
    
    cur.iteration.data = subtype.raw.data
    wd = paste0("/Users/kyle/Desktop/result_fs/var/",subtype)
    setwd(wd)
    
    if (subtype=="aml") {
      k=2
    } else if (subtype=="liver") { 
      k=4
    } else if (subtype=="gbm") {
      k=3
    }
    
    algorithm.ret = ExecuteCC(clusterNum=k, d=cur.iteration.data, plot=plot, maxK=8, clusterAlg="hc", distance="pearson")
    algorithm.ret$clustering.name = paste("clustering_", subtype, "_cc")
    
    result_cc[[j]]=list(clustering.name=algorithm.ret$clustering.name, clustering=algorithm.ret$clustering, 
                        timing=algorithm.ret$timing)
    j=j+1
  } 
  save(result_cc, file="/Users/kyle/Desktop/result_fs/var/result_cc.RData")
}
run.cc.mad <- function(plot) {
  result_cc = list()
  j=1
  for (i in 1:length(SUBTYPES.DATA)) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    subtype.raw.data = get.raw.data(subtype, only.primary=current.subtype.data$only.primary)
    subtype.raw.data = set.omics.list.attr(subtype.raw.data, current.subtype.data)
    subtype.raw.data = log.and.normalize(subtype.raw.data, subtype.data)
    subtype.raw.data = as.matrix(subtype.raw.data[[1]])
    subtype.raw.data = FSbyMAD(subtype.raw.data, cut.type="topk", 2000)
    
    cur.iteration.data = subtype.raw.data
    wd = paste0("/Users/kyle/Desktop/result_fs/mad/",subtype)
    setwd(wd)
    
    if (subtype=="aml") {
      k=2
    } else if (subtype=="liver") { 
      k=3
    } else if (subtype=="gbm") {
      k=5
    }
    
    algorithm.ret = ExecuteCC(clusterNum=k, d=cur.iteration.data, plot=plot, maxK=8, clusterAlg="hc", distance="pearson")
    algorithm.ret$clustering.name = paste("clustering_", subtype, "_cc")
    
    result_cc[[j]]=list(clustering.name=algorithm.ret$clustering.name, clustering=algorithm.ret$clustering, 
                        timing=algorithm.ret$timing)
    j=j+1
  } 
  save(result_cc, file="/Users/kyle/Desktop/result_fs/mad/result_cc.RData")
}

log.and.norm <- function(omics.data, subtype, normalize=T, filter.var=F) {
  # filter features with no variance at all
  if (subtype!="gbm") {
    omics.data[[1]] = log(1+omics.data[[1]])
  }
  if (normalize) {
    omics.data = lapply(omics.data, normalize.matrix)    
  }
  return(omics.data)
}

run.cc.fs <- function(plot) {
  result_cc = list()
  for (k in 1:length(FS)) {
    fs=FS[k]
    for (i in 1:length(SUBTYPES.DATA)) {
      current.subtype.data = SUBTYPES.DATA[[i]]
      subtype = current.subtype.data$name
      wd = paste0("/Users/kyle/Downloads/data/fs/", subtype)
      setwd(wd)
      if (fs=="var") {
        load("var2000")
        cur.iteration.data=var2000
        if (subtype=="aml") {
          K=6
        } else if (subtype=="gbm") { 
          K=3
        } else if (subtype=="breast") {
          K=3
        } else if (subtype=="colon") {
          K=3
        }
      } else if (fs=="mad") {
        load("mad2000")
        cur.iteration.data=mad2000
        if (subtype=="aml") {
          K=3
        } else if (subtype=="gbm") { 
          K=4
        } else if (subtype=="breast") {
          K=3
        } else if (subtype=="colon") {
          K=4
        }
      } else if (fs=="med") {
        load("med2000")
        cur.iteration.data=med2000
        if (subtype=="aml") {
          K=2
        } else if (subtype=="gbm") { 
          K=3
        } else if (subtype=="breast") {
          K=3
        } else if (subtype=="colon") {
          K=4
        }
      } else if (fs=="dip") {
        load("dip2000")
        cur.iteration.data=dip2000
        if (subtype=="aml") {
          K=5
        } else if (subtype=="gbm") { 
          K=3
        } else if (subtype=="breast") {
          K=3
        } else if (subtype=="colon") {
          K=4
        }
      }
      
      cur.iteration.data = log.and.norm(cur.iteration.data, subtype=subtype)
      cur.iteration.data = cur.iteration.data[[1]]
      
      wd = paste0("/Users/kyle/Desktop/result_fs/",fs,"/",subtype)
      setwd(wd)
      
      algorithm.ret = ExecuteCC(clusterNum=K, d=cur.iteration.data, plot=plot, maxK=8, clusterAlg="hc", distance="pearson")
      algorithm.ret$clustering.name = paste0("clustering_", subtype, "_cc")
      
      result_cc[[i]]=list(clustering.name=algorithm.ret$clustering.name, clustering=algorithm.ret$clustering, 
                          timing=algorithm.ret$timing)
    }
    save(result_cc, file=paste0("/Users/kyle/Desktop/result_fs/",fs,"/result_cc.RData"))
  }
}
run.cc.fs(plot="png")
load("/Users/kyle/Desktop/result_fs/var/result_cc.RData")
View(result_cc)

############## iCluster with FS
createdata.fs.icluster <- function() {
  for (i in 1:length(SUBTYPES.DATA)) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    subtype.raw.data = get.raw.data(subtype, only.primary=current.subtype.data$only.primary)
    subtype.raw.data = set.omics.list.attr(subtype.raw.data, current.subtype.data)
    subtype.raw.data = log.and.normalize(subtype.raw.data, subtype.data)
    subtype.raw.data = as.matrix(subtype.raw.data[[1]])
    var2000 = FSbyVar(subtype.raw.data, cut.type="topk", 2000)
    save(var2000, file=paste0("/Users/kyle/Downloads/data/fs/icluster/",subtype,"/var2000"))
    mad2000 = FSbyMAD(subtype.raw.data, cut.type="topk", 2000)
    save(mad2000, file=paste0("/Users/kyle/Downloads/data/fs/icluster/",subtype,"/mad2000"))
    dip2000 = FSbyDip(subtype.raw.data, cut.type = "topk", value=2000)
    save(dip2000, file=paste0("/Users/kyle/Downloads/data/fs/icluster/", subtype,"/dip2000"))
    
    # survival data
    survival.file.path = file.path("/Users/kyle/Downloads/data", subtype, 'survival')
    survival.data = read.table(survival.file.path, header = TRUE)
    patient.names = colnames(subtype.raw.data)
    patient.names.in.file = as.character(survival.data[, 1])
    survival.data$PatientID = toupper(gsub('-', '\\.', substring(patient.names.in.file, 1, 12)))
    stopifnot(all(patient.names %in% survival.data$PatientID))
    indices = match(patient.names, survival.data$PatientID)
    ordered.surv = survival.data[indices,]
  }
}
createdata.fs.icluster.unnorm <- function() {
  for (i in 1:length(SUBTYPES.DATA)) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    subtype.raw.data = get.raw.data(subtype, only.primary=current.subtype.data$only.primary)
    subtype.raw.data = set.omics.list.attr(subtype.raw.data, current.subtype.data)
    subtype.raw.data = log.and.normalize(subtype.raw.data, subtype.data, normalize=F)
    subtype.raw.data = as.matrix(subtype.raw.data[[1]])
    var2000_unnorm = FSbyVar(subtype.raw.data, cut.type="topk", 2000)
    save(var2000_unnorm, file=paste0("/Users/kyle/Downloads/data/fs/icluster/",subtype,"/var2000_unnorm"))
    mad2000_unnorm = FSbyMAD(subtype.raw.data, cut.type="topk", 2000)
    save(mad2000_unnorm, file=paste0("/Users/kyle/Downloads/data/fs/icluster/",subtype,"/mad2000_unnorm"))
    coex2000_unnorm = FSbyCoEX(subtype.raw.data, dataname=subtype, cut.type = "topk", value=2000)
    save(coex2000_unnorm, file=paste0("/Users/kyle/Downloads/data/fs/icluster/", subtype,"/coex2000_unnorm"))
    dip2000_unnorm = FSbyDip(subtype.raw.data, cut.type = "topk", value=2000)
    save(dip2000_unnorm, file=paste0("/Users/kyle/Downloads/data/fs/icluster/", subtype,"/dip2000_unnorm"))
    
    # survival data
    survival.file.path = file.path("/Users/kyle/Downloads/data", subtype, 'survival')
    survival.data = read.table(survival.file.path, header = TRUE)
    patient.names = colnames(subtype.raw.data)
    patient.names.in.file = as.character(survival.data[, 1])
    survival.data$PatientID = toupper(gsub('-', '\\.', substring(patient.names.in.file, 1, 12)))
    stopifnot(all(patient.names %in% survival.data$PatientID))
    indices = match(patient.names, survival.data$PatientID)
    ordered.surv = survival.data[indices,]
    
    cox5_unnorm = FSbyCox(subtype.raw.data, ordered.surv$Survival, ordered.surv$Death, cutoff = 0.05)
    save(cox5_unnorm, file=paste0("/Users/kyle/Downloads/data/fs/icluster/",subtype,"/cox5_unnorm"))
  }
}

tune.icluster.fs <- function() {
  for (i in 1:length(SUBTYPES.DATA)) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    wd = paste0("/Users/kyle/downloads/data/fs/", subtype)
    setwd(wd)
    
    for (k in 1:length(FS)) {
      fs=FS[k]  
      if (fs=="var") {
        load("var2000")
        cur.iteration.data=var2000[[1]]
      } else if (fs=="mad") {
        load("mad2000")
        cur.iteration.data=mad2000[[1]]
      } else if (fs=="med") {
        load("med2000")
        cur.iteration.data=med2000[[1]]
      } else if (fs=="coex") {
        load("coex2000")
        cur.iteration.data=coex2000[[1]]
      } else if (fs=="dip") {
        load("dip2000")
        cur.iteration.data=dip2000[[1]]
      }
      cur.iteration.data=t(cur.iteration.data)
      tune = iClusterPlus::tune.iClusterBayes(cpus=8, dt1=cur.iteration.data, K=1:8)
      tune_wd = paste0("/Users/kyle/Desktop/result_fs/",fs,"/", subtype,"/icluster/tune")
      save(tune, file=tune_wd)
      
      mat = matrix(0, nrow=8, ncol=2)
      for (i in 1:length(tune$fit)) {
        mat[i,1] = tune$fit[[i]]$BIC
        mat[i,2] = tune$fit[[i]]$dev.ratio
      }
      mat_wd = paste0("/Users/kyle/Desktop/result_fs/",fs,"/", subtype,"/icluster/mat")
      save(mat, file=mat_wd)
    }
  }
}
run.icluster.fs <- function() {
  result_icluster = list()
  for (k in 1:length(FS)) {
    fs=FS[k]
    for (i in 1:length(SUBTYPES.DATA)) {
      current.subtype.data = SUBTYPES.DATA[[i]]
      subtype = current.subtype.data$name
      wd = paste0("/Users/kyle/Downloads/data/fs/", subtype)
      setwd(wd)
      if (fs=="var") {
        load("var2000")
        cur.iteration.data=var2000
        if (subtype=="aml") {
          K=5
        } else if (subtype=="breast") { 
          K=3
        } else if (subtype=="colon") {
          K=4
        } else if (subtype=="gbm") {
          K=3
        }
      } else if (fs=="mad") {
        load("mad2000")
        cur.iteration.data=mad2000
        if (subtype=="aml") {
          K=4
        } else if (subtype=="breast") { 
          K=3
        } else if (subtype=="colon") {
          K=4
        } else if (subtype=="gbm") {
          K=3
        }
      } else if (fs=="med") {
        load("med2000")
        cur.iteration.data=med2000
        if (subtype=="aml") {
          K=5
        } else if (subtype=="breast") { 
          K=3
        } else if (subtype=="colon") {
          K=4
        } else if (subtype=="gbm") {
          K=3
        }
      } else if (fs=="dip") {
        load("dip2000")
        cur.iteration.data=dip2000
        if (subtype=="aml") {
          K=3
        } else if (subtype=="breast") { 
          K=3
        } else if (subtype=="colon") {
          K=3
        } else if (subtype=="gbm") {
          K=3
        }
      }
      cur.iteration.data = log.and.norm(cur.iteration.data, subtype=subtype)
      cur.iteration.data=t(cur.iteration.data[[1]])
      
      start = Sys.time()
      result = iClusterPlus::iClusterBayes(dt1=cur.iteration.data, K=K)
      result$clustering.name = paste0("clustering_", subtype,"_icluster")
      time.taken = as.numeric(Sys.time() - start, units='secs')
      result_icluster[[i]]=list(clustering.name=result$clustering.name, clustering=result$clusters,timing=time.taken)
    }
    save(result_icluster, file=paste0("/Users/kyle/Desktop/result_fs/",fs,"/result_icluster.RData"))
  }
}

tune.icluster.fs()
run.icluster.fs()


##################################### NMF, PINS, SNF, NEMO ###################################################
#### SNF & NEMO w/o FS
ALGORITHM.NAMES = c('snf', 'nemo')
ALGORITHM.DISPLAY.NAMES = as.list(c('SNF', 'NEMO'))
run.benchmark.sn()
# SNF patient attribute
for (i in c(1,3,5)) {
  attr(result_sn[[i]]$clustering, "names") = names(result_sn[[i+1]]$clustering)
}
View(result_sn)

#### NMF & PINS
ALGORITHM.NAMES = c('nmf', 'pins')
ALGORITHM.DISPLAY.NAMES = as.list(c('NMF', 'PINS'))
run.benchmark.np()
# NMF patient attribute
#for (i in c(1,3,5)) attr(result_np[[i]]$clustering, "names") = names(result_np[[i+1]]$clustering)

#result_nmf/pins_unnorm=result_np <= run.nmf, run.pins 에서 normalize=F 해서 run.benchmark.np() 돌린거
result_nmf_unnorm=list()
result_pins_unnorm=list()
for (i in 1:3) {
  result_nmf_unnorm[[i]] = result_np_unnorm[[2*i-1]]
  result_pins_unnorm[[i]] = result_np_unnorm[[2*i]]
}

#### PINS, SNF, NEMO with FS
createdata.fs <- function() {
  for (i in 1:length(SUBTYPES.DATA)) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    subtype.raw.data = get.raw.data(subtype, only.primary=current.subtype.data$only.primary)
    subtype.raw.data = set.omics.list.attr(subtype.raw.data, current.subtype.data)
    subtype.raw.data = log.and.normalize(subtype.raw.data, subtype.data)
    subtype.raw.data = as.matrix(subtype.raw.data[[1]])
    start1 = Sys.time()
    var2000 = list(FSbyVar(subtype.raw.data, cut.type="topk", 2000))
    time.taken1 = as.numeric(Sys.time() - start1, units='secs')
    save(var2000, file=paste0("/Users/kyle/Downloads/data/fs/",subtype,"/var2000"))
    
    start2 = Sys.time()
    mad2000 = list(FSbyMAD(subtype.raw.data, cut.type="topk", 2000)) 
    time.taken2 = as.numeric(Sys.time() - start2, units='secs')
    save(mad2000, file=paste0("/Users/kyle/Downloads/data/fs/",subtype,"/mad2000"))
    
    start3 = Sys.time()
    dip2000 = FSbyDip(subtype.raw.data, cut.type = "topk", 2000)
    time.taken3 = as.numeric(Sys.time() - start3, units='secs')
    save(dip2000, file=paste0("/Users/kyle/Downloads/data/fs/", subtype,"/dip2000"))
    
    # survival data
    survival.file.path = file.path("/Users/kyle/Downloads/data", subtype, 'survival')
    survival.data = read.table(survival.file.path, header = TRUE)
    patient.names = colnames(subtype.raw.data)
    patient.names.in.file = as.character(survival.data[, 1])
    survival.data$PatientID = toupper(gsub('-', '\\.', substring(patient.names.in.file, 1, 12)))
    stopifnot(all(patient.names %in% survival.data$PatientID))
    indices = match(patient.names, survival.data$PatientID)
    ordered.surv = survival.data[indices,]
    
    cox5 = list(FSbyCox(subtype.raw.data, ordered.surv$Survival, ordered.surv$Death, cutoff = 0.05))
    save(cox5, file=paste0("/Users/kyle/Downloads/data/fs/",subtype,"/cox5"))
  }
}
createdata.fs.unnorm <- function() {
  for (i in 1:length(SUBTYPES.DATA)) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    subtype.raw.data = get.raw.data(subtype, only.primary=current.subtype.data$only.primary)
    subtype.raw.data = set.omics.list.attr(subtype.raw.data, current.subtype.data)
    subtype.raw.data = log.and.normalize(subtype.raw.data, subtype.data, normalize=FALSE)
    subtype.raw.data = as.matrix(subtype.raw.data[[1]])
    var2000 = list(FSbyVar(subtype.raw.data, cut.type="topk", 2000)) 
    save(var2000, file=paste0("/Users/kyle/Downloads/data/fs/",subtype,"/var2000_unnorm"))
    mad2000 = list(FSbyMAD(subtype.raw.data, cut.type="topk", 2000)) 
    save(mad2000, file=paste0("/Users/kyle/Downloads/data/fs/",subtype,"/mad2000_unnorm"))
    #pca9 = list(FSbyPCA(subtype.raw.data, PC_percent = 0.9, scale=TRUE))
    #save(pca9, file=paste0("/Users/kyle/Downloads/data/fs/",subtype,"/pca9_unnorm"))
    coex2000 = FSbyCoEX(subtype.raw.data, dataname=subtype, cut.type = "topk", 2000)
    save(coex2000, file=paste0("/Users/kyle/Downloads/data/fs/", subtype,"/coex2000_unnorm"))
    dip2000 = FSbyDip(subtype.raw.data, cut.type = "topk", 2000)
    save(dip2000, file=paste0("/Users/kyle/Downloads/data/fs/", subtype,"/dip2000_unnorm"))
    
    # survival data
    survival.file.path = file.path("/Users/kyle/Downloads/data", subtype, 'survival')
    survival.data = read.table(survival.file.path, header = TRUE)
    patient.names = colnames(subtype.raw.data)
    patient.names.in.file = as.character(survival.data[, 1])
    survival.data$PatientID = toupper(gsub('-', '\\.', substring(patient.names.in.file, 1, 12)))
    stopifnot(all(patient.names %in% survival.data$PatientID))
    indices = match(patient.names, survival.data$PatientID)
    ordered.surv = survival.data[indices,]
    
    cox5 = list(FSbyCox(subtype.raw.data, ordered.surv$Survival, ordered.surv$Death, cutoff = 0.05))
    save(cox5, file=paste0("/Users/kyle/Downloads/data/fs/",subtype,"/cox5_unnorm"))
  }
}
createdata.fs()
createdata.fs.unnorm()

ALGORITHM.NAMES = c('nmf', 'pins', 'snf', 'nemo')
run.benchmark.fs <- function() {
  for (k in 1:length(FS)) {
    fs=FS[k]
    j=1
    result_nmf=list()
    result_pins=list()
    result_snf=list()
    result_nemo=list()
    
    for (i in 1:length(SUBTYPES.DATA)) {
      current.subtype.data = SUBTYPES.DATA[[i]]
      subtype = current.subtype.data$name
      wd = paste0("/Users/suzuk/Downloads/data/fs/", subtype)
      setwd(wd)
      
      if (fs=="var") {
        load("var2000")
        cur.iteration.data=var2000
      } else if (fs=="mad") {
        load("mad2000")
        cur.iteration.data=mad2000
      } else if (fs=="med") {
        load("med2000")
        cur.iteration.data=med2000
      } else if (fs=="dip") {
        load("dip2000")
        cur.iteration.data=dip2000
      }
      
      for (algorithm.name in ALGORITHM.NAMES) {
        set.seed(42)
        print(paste('subtype', subtype, 'running algorithm', algorithm.name, 'fs', fs))
        algorithm.func.name = paste0('run.', algorithm.name, '.fs')
        algorithm.func = get(algorithm.func.name)
        clustering.name = paste0('clustering_', subtype,'_', algorithm.name)
        if (algorithm.name!="nmf") {
          cur.iteration.data = log.and.norm(cur.iteration.data, subtype=subtype)
          algorithm.ret = algorithm.func(cur.iteration.data, current.subtype.data)
        } else {
          cur.iteration.data = log.and.norm(cur.iteration.data, subtype=subtype, normalize = F)
          algorithm.ret = algorithm.func(cur.iteration.data, current.subtype.data, fs=fs)
        }
        
        if (algorithm.name=="pins") {
          result_pins[[j]]=list(clustering.name=clustering.name, clustering=algorithm.ret$clustering, timing=algorithm.ret$timing)
        } else if (algorithm.name=="snf") {
          result_snf[[j]]=list(clustering.name=clustering.name, clustering=algorithm.ret$clustering, timing=algorithm.ret$timing)
        } else if (algorithm.name=="nemo") {
          result_nemo[[j]]=list(clustering.name=clustering.name, clustering=algorithm.ret$clustering, timing=algorithm.ret$timing)
        } else if (algorithm.name=="nmf") {
          result_nmf[[j]]=list(clustering.name=clustering.name, clustering=algorithm.ret$clustering, timing=algorithm.ret$timing)
        } 
      }
      j=j+1
    }
    save(result_nmf, file=paste0("/Users/kyle/Desktop/result_fs/",fs,"/result_nmf.RData"))
    save(result_pins, file=paste0("/Users/kyle/Desktop/result_fs/",fs,"/result_pins.RData"))
    save(result_snf, file=paste0("/Users/kyle/Desktop/result_fs/",fs,"/result_snf.RData"))
    save(result_nemo, file=paste0("/Users/kyle/Desktop/result_fs/",fs,"/result_nemo.RData"))
  }
}
run.benchmark.fs()

##################################### Evaluation ################################################################
ALGORITHM.NAMES = c('cc', 'pins', 'pins_unnorm', 'nmf', 'icluster', 'icluster_unnorm', 'snf', 'nemo')
ALGORITHM.DISPLAY.NAMES = as.list(c('CC', 'PINS', 'PINS_UNNORM', 'NMF', 'iClusterBayes', 
                                    'iClusterBayes_UNNORM', 'SNF', 'NEMO'))

for (i in 1:3) result_pins_unnorm[[i]]$clustering.name = paste0(result_pins_unnorm[[i]]$clustering.name,'_unnorm')
for (i in 1:3) result_icluster_unnorm[[i]]$clustering.name = paste0(result_icluster_unnorm[[i]]$clustering.name,'_unnorm')

#### 데이터 결합 (w/o FS)
result_all=list()
for (i in 1:3) result_all[[8*(i-1)+1]] = result_cc[[i]]
for (i in 1:3) result_all[[8*(i-1)+2]] = result_pins[[i]]
for (i in 1:3) result_all[[8*(i-1)+3]] = result_pins_unnorm[[i]]
for (i in 1:3) result_all[[8*(i-1)+4]] = result_nmf_unnorm[[i]]
for (i in 1:3) result_all[[8*(i-1)+5]] = result_icluster[[i]]
for (i in 1:3) result_all[[8*(i-1)+6]] = result_icluster_unnorm[[i]]
for (i in 1:3) result_all[[8*(i-1)+7]] = result_sn[[2*i-1]]
for (i in 1:3) result_all[[8*i]] = result_sn[[2*i]]
#save(result_all, file="/Users/kyle/Desktop/result_final/result_all.RData")
load("/Users/kyle/Desktop/result_final/result_all.RData")

all.result = analyze.benchmark(result_all)
#result_time = benchmark.omics.time(all.result)
#result_num.clusters = benchmark.omics.num.clusters(all.result)
#benchmark.omics.surv(all.result)
perform.all.analyses(all.result)

#### 데이터 결합 (with FS)
for (k in 1:length(FS)) {
  fs=FS[k]
  wd=paste0("/Users/kyle/Desktop/result_fs/",fs)
  setwd(wd)
  load("result_cc.RData")
  load("result_pins.RData")
  load("result_pins_unnorm.RData")
  load("result_nmf_unnorm.RData")
  load("result_icluster.RData")
  load("result_icluster_unnorm.RData")
  load("result_snf.RData")
  load("result_nemo.RData")
  
  result_all=list()
  for (i in 0:2) {
    result_all[[(8*i)+1]] = result_cc[[i+1]]
    result_all[[(8*i)+2]] = result_pins[[i+1]]
    result_all[[(8*i)+3]] = result_pins_unnorm[[i+1]]
    result_all[[(8*i)+3]]$clustering.name = paste0(result_all[[(8*i)+3]]$clustering.name,"_unnorm")
    if (fs!="pca") {
      result_all[[(8*i)+4]] = result_nmf_unnorm[[i+1]]
    } else result_all[[(8*i)+4]] = NULL
    result_all[[(8*i)+5]] = result_icluster[[i+1]]
    result_all[[(8*i)+6]] = result_icluster_unnorm[[i+1]]
    result_all[[(8*i)+7]] = result_snf[[i+1]]
    result_all[[(8*i)+8]] = result_nemo[[i+1]]
  }
  
  # patient attribute
  if (fs!="pca") {
    for (i in c(2:7)) {
      attr(result_all[[i]]$clustering, "names") = names(result_all[[1]]$clustering)
    }
    for (i in c(10:15)) {
      attr(result_all[[i]]$clustering, "names") = names(result_all[[9]]$clustering)
    }
    for (i in c(18:23)) {
      attr(result_all[[i]]$clustering, "names") = names(result_all[[17]]$clustering)
    }
  } else if (fs=="pca") {
    for (i in c(2,3,5:7)) {
      attr(result_all[[i]]$clustering, "names") = names(result_all[[1]]$clustering)
    }
    for (i in c(10,11,13:15)) {
      attr(result_all[[i]]$clustering, "names") = names(result_all[[9]]$clustering)
    }
    for (i in c(18,19,21:23)) {
      attr(result_all[[i]]$clustering, "names") = names(result_all[[17]]$clustering)
    }
  }
  save(result_all, file=paste0("/Users/kyle/Desktop/result_fs/",fs,"/result_all.RData"))
}
load("/Users/kyle/Desktop/result_fs/var/result_all.RData")
#load("/Users/kyle/Desktop/result_fs/mad/result_all.RData")
View(result_all)

analyze.benchmark.fs <- function(result_data) {
  all.clusterings = list()
  all.timings = list()
  j=1
  for (i in 1:length(SUBTYPES.DATA)) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    #subtype.raw.data = get.raw.data(subtype,only.primary=current.subtype.data$only.primary)
    
    all.clusterings[[subtype]] = list()
    all.timings[[subtype]] = list()
    
    for (algorithm.name in ALGORITHM.NAMES) {
      all.clusterings[[subtype]][[algorithm.name]] = list()
      all.timings[[subtype]][[algorithm.name]] = list()
      
      clustering = result_data[[j]]$clustering
      timing = result_data[[j]]$timing
      
      if (length(clustering)<1) {
        all.clusterings[[subtype]][[algorithm.name]] = paste("PINS cannot be used with NMF")
        all.timings[[subtype]][[algorithm.name]] = paste("PINS cannot be used with NMF")
      } else {
        all.clusterings[[subtype]][[algorithm.name]] = clustering
        all.timings[[subtype]][[algorithm.name]] = timing 
      }
      j=j+1
    }
  }
  return(list(all.clusterings=all.clusterings, all.timings=all.timings))
}
for (k in 1:length(FS)) {
  fs=FS[k]
  wd=paste0("/Users/kyle/Desktop/result_fs/",fs)
  setwd(wd)
  load("result_all.RData")
  
  surv.fs.dir <- function() {
    return(paste0('/Users/kyle/Desktop/result_fs/',fs,"/surv"))
  }
  
  all.result.var = analyze.benchmark.fs(result_all) 
  benchmark.omics.surv.fs(all.result.fs)
}


##################################### GOLD STANDARD #########################################################
SUBTYPES.DATA = list(
  list(name='breast', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='BRCA'),
  list(name='colon', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='COAD'))

brca=read.csv("/Users/kyle/Downloads/data/brca/mRNA.csv",row.names=1)
brca=as.matrix(brca)
coad=read.csv("/Users/kyle/Downloads/data/coad/mRNA.csv",row.names=1)
coad=as.matrix(coad)
run.cc.gold <- function(plot) {
  result_cc_gold = list()
  j=1
  for (i in 1:length(SUBTYPES.DATA)) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    subtype.raw.data = read.csv(paste0("/Users/kyle/Downloads/data/",subtype,"/mRNA.csv"),row.names=1)
    subtype.raw.data = as.matrix(subtype.raw.data)
    
    cur.iteration.data = subtype.raw.data
    wd = paste0("/Users/kyle/Desktop/result_final/",subtype)
    setwd(wd)
    
    if (subtype=="brca") {
      k=3
    } else if (subtype=="coad") { 
      k=4
    }
    
    algorithm.ret = ExecuteCC(clusterNum=k, d=cur.iteration.data, plot=plot, maxK=8, clusterAlg="hc", distance="pearson")
    algorithm.ret$clustering.name = paste("clusteirng_", subtype, "_cc")
    
    result_cc_gold[[j]]=list(clustering.name=algorithm.ret$clustering.name, clustering=algorithm.ret$clustering, 
                             timing=algorithm.ret$timing)
    j=j+1
  } 
  setwd("/Users/kyle/Desktop/result_final")
  save(result_cc_gold, file="result_cc_gold.RData")
}
run.cc.gold(plot=NULL)

run.benchmark()
run.benchmark.gold.nmf <- function() {
  j=1
  result_nmf_gold=list()
  
  for (i in 1:length(SUBTYPES.DATA)) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    subtype.raw.data = get.raw.data(subtype, only.primary=current.subtype.data$only.primary)
    subtype.raw.data = set.omics.list.attr(subtype.raw.data, current.subtype.data)
    cur.iteration.data = subtype.raw.data
    
    print(paste('subtype', subtype, 'running algorithm nmf'))
    clustering.name = paste0('clustering_', subtype,'_nmf')
    algorithm.ret = run.nmf(cur.iteration.data, current.subtype.data)
    
    result_nmf_gold[[j]]=list(clustering.name=clustering.name, clustering=algorithm.ret$clustering, timing=algorithm.ret$timing)
    j=j+1
  }
  save(result_nmf_gold, file=paste0("/Users/kyle/Desktop/result_fs/",fs,"/result_nmf_gold.RData"))
}
run.benchmark.gold.nmf()


##################################### RUN #######################################################
#########  10/27
SUBTYPES.DATA = list(
  list(name='aml', only.primary=F, is.rna.seq=T, is.mirna.seq=T, display.name='AML'),
  list(name='gbm', only.primary=T, is.rna.seq=F, is.mirna.seq=F, display.name='GBM'),
  list(name='breast', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='BIC'),
  list(name='colon', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='COAD'))

ALGORITHM.NAMES = c('cc', 'pins', 'nmf', 'icluster', 'snf', 'nemo')
ALGORITHM.DISPLAY.NAMES = as.list(c('CC', 'PINS', 'NMF', 'iClusterBayes', 'SNF', 'NEMO'))

## without FS
setwd("/Users/kyle/Desktop/result_final")
load("result_all(all).RData")
analyze.benchmark <- function(result_data) {
  all.clusterings = list()
  all.timings = list()
  j=1
  for (i in 1:length(SUBTYPES.DATA)) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    #subtype.raw.data = get.raw.data(subtype,only.primary=current.subtype.data$only.primary)
    
    all.clusterings[[subtype]] = list()
    all.timings[[subtype]] = list()
    
    for (algorithm.name in ALGORITHM.NAMES) {
      all.clusterings[[subtype]][[algorithm.name]] = list()
      all.timings[[subtype]][[algorithm.name]] = list()
      
      clustering = result_data[[j]]$clustering
      timing = result_data[[j]]$timing
      
      #if (!any(is.na(clustering))) {
      #  names(clustering) = colnames(subtype.raw.data[[1]])
      #}
      
      all.clusterings[[subtype]][[algorithm.name]] = clustering
      all.timings[[subtype]][[algorithm.name]] = timing
      j=j+1
    }
  }
  return(list(all.clusterings=all.clusterings, all.timings=all.timings))
}
get.logrank.pvalue <- function(survdiff.res) {
  
  1 - pchisq(survdiff.res$chisq, length(survdiff.res$n) - 1)  
  
}
benchmark.omics.surv <- function(benchmark.results) {
  all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  rownames(all.surv.pvalues) = ALGORITHM.NAMES
  colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
  all.clusterings = benchmark.results$all.clusterings
  for (i in 1:length(all.clusterings)) {
    subtype = colnames(all.surv.pvalues)[i]
    subtype.clusterings = all.clusterings[[subtype]]
    for (j in 1:length(subtype.clusterings)) {
      surv.path = file.path("/Users/kyle/Desktop/result_final/surv", paste(subtype, ALGORITHM.NAMES[j], 'surv', sep='_'))
      load(surv.path)
      pvalue = get.logrank.pvalue(check.survival(clustering, subtype))
      all.surv.pvalues[j, i] = pvalue
    }
  }
  return(all.surv.pvalues)
}
benchmark.omics.time <- function(benchmark.results) {
  
  all.alg.times = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  
  rownames(all.alg.times) = ALGORITHM.NAMES
  
  colnames(all.alg.times) = sapply(SUBTYPES.DATA, function(x) x$name)
  
  all.timings = benchmark.results$all.timings
  
  for (i in 1:length(all.timings)) {
    
    subtype = colnames(all.alg.times)[i]
    
    subtype.timings = all.timings[[subtype]]
    
    for (j in 1:length(subtype.timings)) {
      
      timing = subtype.timings[[ALGORITHM.NAMES[j]]]
      
      all.alg.times[j, i] = timing
      
    }
    
  }
  
  return(all.alg.times)
  
}
benchmark.omics.num.clusters <- function(benchmark.results) {
  
  num.clusters = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  
  rownames(num.clusters) = ALGORITHM.NAMES
  
  colnames(num.clusters) = sapply(SUBTYPES.DATA, function(x) x$name)
  
  all.clusterings = benchmark.results$all.clusterings
  
  for (i in 1:length(all.clusterings)) {
    
    subtype = colnames(num.clusters)[i]
    
    subtype.clusterings = all.clusterings[[subtype]]
    
    for (j in 1:length(subtype.clusterings)) {
      
      clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
      
      num.clusters[j, i] = max(clustering)
      
    }
    
  }
  
  return(num.clusters)
  
}
check.survival <- function(groups, subtype, survival.file.path) {
  
  if (missing(survival.file.path)) {
    
    survival.file.path = get.subtype.survival.path(subtype)
    
  }
  
  survival.data = read.table(survival.file.path, header = TRUE)
  
  patient.names = names(groups)
  
  patient.names.in.file = as.character(survival.data[, 1])
  
  patient.names.in.file = toupper(gsub('-', '\\.', substring(patient.names.in.file, 1, 12)))
  
  
  
  stopifnot(all(patient.names %in% patient.names.in.file))
  
  
  
  indices = match(patient.names, patient.names.in.file)
  
  ordered.survival.data = survival.data[indices,]
  
  ordered.survival.data["cluster"] <- groups
  
  ordered.survival.data$Survival[is.na(ordered.survival.data$Survival)] = 0
  
  ordered.survival.data$Death[is.na(ordered.survival.data$Death)] = 0
  
  return(survdiff(Surv(Survival, Death) ~ cluster, data=ordered.survival.data))
  
}
get.subtype.survival.path <- function(subtype) {
  survival.file.path=file.path("/Users/kyle/Downloads/data",subtype,'survival')
}
benchmark.omics.surv <- function(benchmark.results) {
  all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  rownames(all.surv.pvalues) = ALGORITHM.NAMES
  colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
  all.clusterings = benchmark.results$all.clusterings
  for (i in 1:length(all.clusterings)) {
    subtype = colnames(all.surv.pvalues)[i]
    subtype.clusterings = all.clusterings[[subtype]]
    for (j in 1:length(subtype.clusterings)) {
      surv.path = file.path("/Users/kyle/Desktop/result_final/surv", paste(subtype, ALGORITHM.NAMES[j], 'surv', sep='_'))
      
      if (file.exists(surv.path)) {
        load(surv.path)
        pvalue = empirical.surv.ret$pvalue
      } else {
        clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
        if (length(table(clustering)) > 1) {
          #pvalue = -log10(get.logrank.pvalue(check.survival(clustering, subtype)))
          empirical.surv.ret = get.empirical.surv(clustering, subtype)
          save(empirical.surv.ret, file=surv.path)
          pvalue = empirical.surv.ret$pvalue
        } else {
          pvalue = NA
        }  
      }
      all.surv.pvalues[j, i] = pvalue
    }
  }
  return(all.surv.pvalues)
}

all.result = analyze.benchmark(result_all)
benchmark.omics.surv(all.result)
perform.all.analyses(all.result)
benchmark.sil(result_all,fs="no")

## with FS (log, normalize 안 한 상태에서 돌려야함)
FS=c("var","mad","med","dip")
SUBTYPES.DATA = list(
  list(name='aml', only.primary=F, is.rna.seq=T, is.mirna.seq=T, display.name='AML'),
  list(name='gbm', only.primary=T, is.rna.seq=F, is.mirna.seq=F, display.name='GBM'))
createdata.fs()

SUBTYPES.DATA = list(
  list(name='breast', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='BIC'),
  list(name='colon', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='COAD'))
createdata.fs.gold()

# survival for all data
benchmark.omics.surv.fs <- function(benchmark.results) {
  all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  rownames(all.surv.pvalues) = ALGORITHM.NAMES
  colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
  all.clusterings = benchmark.results$all.clusterings
  for (i in 1:length(all.clusterings)) {
    subtype = colnames(all.surv.pvalues)[i]
    subtype.clusterings = all.clusterings[[subtype]]
    for (j in 1:length(subtype.clusterings)) {
      surv.path = file.path(surv.fs.dir(), paste(subtype, ALGORITHM.NAMES[j], 'surv', sep='_'))
      print(surv.path)
      if (file.exists(surv.path)) {
        load(surv.path)
        pvalue = -log10(empirical.surv.ret$pvalue)
      } else {
        clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
        if (length(table(clustering)) > 1) {
          #pvalue = -log10(get.logrank.pvalue(check.survival(clustering, subtype)))
          empirical.surv.ret = get.empirical.surv(clustering, subtype)
          save(empirical.surv.ret, file=surv.path)
          pvalue = -log10(empirical.surv.ret$pvalue)
        } else {
          pvalue = NA
        }  
      }
      all.surv.pvalues[j, i] = pvalue
    }
  }
  return(all.surv.pvalues)
}
benchmark.omics.surv.fs.re <- function(benchmark.results,fs) {
  all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  rownames(all.surv.pvalues) = ALGORITHM.NAMES
  colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
  all.clusterings = benchmark.results$all.clusterings
  for (i in 1:length(all.clusterings)) {
    subtype = colnames(all.surv.pvalues)[i]
    subtype.clusterings = all.clusterings[[subtype]]
    for (j in 1:length(subtype.clusterings)) {
      surv.path = file.path(surv.fs.dir(), paste(subtype, ALGORITHM.NAMES[j], 'pval', sep='_'))
      print(surv.path)
      if (file.exists(surv.path)) {
        load(surv.path)
        pvalue = pvalue
      } else {
        clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
        if (length(table(clustering)) > 1) {
          pvalue = -log10(get.logrank.pvalue(check.survival(clustering, subtype)))
          empirical.surv.ret = get.empirical.surv(clustering, subtype)
          save(pvalue, file=surv.path)
          #pvalue = empirical.surv.ret$pvalue
        } else {
          pvalue = NA
        }  
      }
      all.surv.pvalues[j, i] = pvalue
    }
  }
  write.csv(all.surv.pvalues, file=paste0("/Users/kyle/Desktop/result_fs/",fs,"/tables/pval.csv"))
}
analyze.benchmark.fs <- function(result_data) {
  all.clusterings = list()
  all.timings = list()
  j=1
  for (i in 1:length(SUBTYPES.DATA)) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    #subtype.raw.data = get.raw.data(subtype,only.primary=current.subtype.data$only.primary)
    
    all.clusterings[[subtype]] = list()
    all.timings[[subtype]] = list()
    
    for (algorithm.name in ALGORITHM.NAMES) {
      all.clusterings[[subtype]][[algorithm.name]] = list()
      all.timings[[subtype]][[algorithm.name]] = list()
      
      clustering = result_data[[j]]$clustering
      timing = result_data[[j]]$timing
      
      all.clusterings[[subtype]][[algorithm.name]] = clustering
      all.timings[[subtype]][[algorithm.name]] = timing 
      j=j+1
    }
  }
  return(list(all.clusterings=all.clusterings, all.timings=all.timings))
}
perform.all.analyses.fs <- function(benchmark.ret, fs) {
  
  for (i in 1:3) {
    cur.func = list(benchmark.omics.time, benchmark.omics.num.clusters,benchmark.omics.surv.fs)[[i]]
    
    benchmark.data = cur.func(benchmark.ret)
    
    displayed.benchmark.data = benchmark.data      
    colnames(displayed.benchmark.data)[1:ncol(displayed.benchmark.data)] = 
      sapply(as.list(colnames(displayed.benchmark.data)[1:ncol(displayed.benchmark.data)]), 
             subtype.to.display.name)
    rownames(displayed.benchmark.data) = ALGORITHM.DISPLAY.NAMES
    print.matrix.latex.format(displayed.benchmark.data)
    
    table.name = c('runtime', 'num_cluster', 'survival')[i]
    
    file.name = deparse(substitute(benchmark.ret))
    file = file.path("/Users/kyle/Desktop/result_fs", paste0(fs,'/tables/',table.name, '.csv'))
    
    write.csv(displayed.benchmark.data, file=file)
    
    print('------------------------')
  }
}
FS=c("var","mad", "med", "dip")
for (k in 1:length(FS)) {
  fs=FS[k]
  wd=paste0("/Users/kyle/Desktop/result_fs/",fs)
  setwd(wd)
  load("result_all.RData")
  
  surv.fs.dir <- function() {
    return(paste0("/Users/kyle/Desktop/result_fs/",fs,"/surv"))
  }
  
  all.result.fs = analyze.benchmark.fs(result_all) 
  benchmark.omics.surv.fs(all.result.fs)
  perform.all.analyses.fs(all.result.fs, fs=fs)
  benchmark.sil(result_all,fs=fs)
}


##################################### result for fs (학회용 발표) ###################
res = read.csv("/Users/kyle/Desktop/result_fs/fs_res.csv",header=T)
View(res)
res$comb = paste0(res$fs,"_",res$clustering)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggrepel)

res_aml = res[res$data=="AML",]
ggplot(res_aml, aes(x = fs, y=log.pval)) +
  geom_point(aes(colour=clustering)) +
  geom_text_repel(res_aml[which(res_aml$log.pval>1.3),], mapping=aes(label=clustering),hjust=-0.2, vjust=1, size=3)+
  geom_hline(yintercept=2.25)
labs(title = "-log(pval) for AML")
ggplot(res_aml, aes(x = fs, y=silhouette)) +
  geom_point(aes(colour=clustering)) +
  geom_text(res_aml[which(res_aml$silhouette>0.15),], mapping=aes(label=clustering),hjust=-0.2, vjust=1, size=3)+
  labs(title = "silhoutte coef for AML")

res_gbm = res[res$data=="GBM",]
ggplot(res_gbm, aes(x = fs, y=log.pval)) +
  geom_point(aes(colour=clustering)) +
  geom_text(res_gbm[which(res_gbm$log.pval>1.3),], mapping=aes(label=clustering),hjust=-0.2, vjust=1, size=3)+
  labs(title = "-log(pval) for GBM")
ggplot(res_gbm, aes(x = fs, y=silhouette)) +
  geom_point(aes(colour=clustering)) +
  geom_text(res_gbm[which(res_gbm$silhouette>0.15),], mapping=aes(label=clustering),hjust=-0.2, vjust=1, size=3)+
  labs(title = "silhouette for GBM")

res_bic = res[res$data=="BIC",]
ggplot(res_bic, aes(x = fs, y=log.pval)) +
  geom_point(aes(colour=clustering)) +
  geom_text_repel(res_bic[which(res_bic$log.pval>1.3),], mapping=aes(label=clustering),hjust=-0.2, vjust=1, size=3)+
  labs(title = "-log(pval) for BIC")
ggplot(res_bic, aes(x = fs, y=silhouette)) +
  geom_point(aes(colour=clustering)) +
  geom_text(res_bic[which(res_bic$silhouette>0.15),], mapping=aes(label=clustering),hjust=-0.2, vjust=1, size=3)+
  labs(title = "silhouette for BIC")

res_coad = res[res$data=="COAD",]
ggplot(res_coad, aes(x = fs, y=log.pval)) +
  geom_point(aes(colour=clustering)) +
  geom_text(res_coad[which(res_coad$log.pval>1.3),], mapping=aes(label=clustering),hjust=-0.2, vjust=1, size=3)+
  labs(title = "-log(pval) for COAD")
ggplot(res_coad, aes(x = fs, y=silhouette)) +
  geom_point(aes(colour=clustering)) +
  geom_text(res_coad[which(res_coad$silhouette>0.14),], mapping=aes(label=clustering),hjust=-0.2, vjust=1, size=3)+
  labs(title = "silhouette for COAD")

##################################### KM Curve (w/o FS) ######################
library(ggplot2) 
library(survminer)
check.curve <- function(groups, subtype, survival.file.path) {
  if (missing(survival.file.path)) {
    survival.file.path = get.subtype.survival.path(subtype)
  }
  survival.data = read.table(survival.file.path, header = TRUE)
  patient.names = names(groups)
  patient.names.in.file = as.character(survival.data[, 1])
  patient.names.in.file = toupper(gsub('-', '\\.', substring(patient.names.in.file, 1, 12)))
  
  stopifnot(all(patient.names %in% patient.names.in.file))
  
  indices = match(patient.names, patient.names.in.file)
  ordered.survival.data = survival.data[indices,]
  ordered.survival.data["cluster"] <- groups
  ordered.survival.data$Survival[is.na(ordered.survival.data$Survival)] = 0
  ordered.survival.data$Death[is.na(ordered.survival.data$Death)] = 0
  return(survfit(Surv(Survival, Death) ~ cluster, data=ordered.survival.data))
}

load("/Users/kyle/Desktop/result_final/result_all(all).RData")
all.result = analyze.benchmark.fs(result_all)
all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
rownames(all.surv.pvalues) = ALGORITHM.NAMES
colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
all.clusterings = all.result$all.clusterings
for (i in 1:length(all.clusterings)) {
  subtype = colnames(all.surv.pvalues)[i]
  subtype.clusterings = all.clusterings[[subtype]]
  curve=list()
  for (j in 1:length(subtype.clusterings)) {
    clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
    fit = check.curve(clustering, subtype)
    pvalue = round(get.logrank.pvalue(check.survival(clustering, subtype)), digits = 4)
    curve[[j]] = ggsurvplot(fit, pval=pvalue, pval.coord = c(0, 0.05), legend.title="", font.legend = list(size=10),
                            title=paste0("KM Curve for ",subtype, "_", ALGORITHM.NAMES[j])) +
      guides(colour = guide_legend(nrow = 2))
  }
  arrange_ggsurvplots(curve, ncol=3, nrow=2)
}

######## KM Curve (FS2000) 
check.curve <- function(groups, subtype, survival.file.path) {
  if (missing(survival.file.path)) {
    survival.file.path = get.subtype.survival.path(subtype)
  }
  survival.data = read.table(survival.file.path, header = TRUE)
  patient.names = names(groups)
  patient.names.in.file = as.character(survival.data[, 1])
  patient.names.in.file = toupper(gsub('-', '\\.', substring(patient.names.in.file, 1, 12)))
  
  stopifnot(all(patient.names %in% patient.names.in.file))
  
  indices = match(patient.names, patient.names.in.file)
  ordered.survival.data = survival.data[indices,]
  ordered.survival.data["cluster"] <- groups
  ordered.survival.data$Survival[is.na(ordered.survival.data$Survival)] = 0
  ordered.survival.data$Death[is.na(ordered.survival.data$Death)] = 0
  return(survfit(Surv(Survival, Death) ~ cluster, data=ordered.survival.data))
}

FS=c("var","mad", "med", "dip")

k=2 # change
fs=FS[k]
wd=paste0("/Users/kyle/Desktop/result_fs/",fs)
setwd(wd)
load("result_all.RData")

all.result.fs = analyze.benchmark.fs(result_all) 
all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
rownames(all.surv.pvalues) = ALGORITHM.NAMES
colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
all.clusterings = all.result.fs$all.clusterings

i=4 # change
subtype = colnames(all.surv.pvalues)[i]
subtype.clusterings = all.clusterings[[subtype]]
curve = list()
for (j in 1:length(subtype.clusterings)) {
  clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
  fit = check.curve(clustering, subtype)
  pvalue = round(get.logrank.pvalue(check.survival(clustering, subtype)), digits = 3)
  curve[[j]] = ggsurvplot(fit, pval=pvalue, pval.coord = c(0, 0.05), legend.title="", font.legend = list(size=10),
                          title=paste0("KM Curve for ",subtype, "_", fs, "_", ALGORITHM.NAMES[j])) +
    guides(colour = guide_legend(nrow = 1))
}
arrange_ggsurvplots(curve, ncol=3, nrow=2)


######## KM Curve (FS500)
k=4 # change
fs=FS[k]
wd=paste0("/Users/kyle/Desktop/result_fs_500/",fs)
setwd(wd)
load("result_all.RData")

all.result.fs = analyze.benchmark.fs(result_all) 
all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
rownames(all.surv.pvalues) = ALGORITHM.NAMES
colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
all.clusterings = all.result.fs$all.clusterings

i=4 # change
subtype = colnames(all.surv.pvalues)[i]
subtype.clusterings = all.clusterings[[subtype]]
curve = list()
for (j in 1:length(subtype.clusterings)) {
  clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
  fit = check.curve(clustering, subtype)
  pvalue = round(get.logrank.pvalue(check.survival(clustering, subtype)), digits = 3)
  curve[[j]] = ggsurvplot(fit, pval=pvalue, pval.coord = c(0, 0.05), legend.title="", font.legend = list(size=10),
                          title=paste0("KM Curve for ",subtype, "_", fs, "_", ALGORITHM.NAMES[j])) +
    guides(colour = guide_legend(nrow = 1))
}
arrange_ggsurvplots(curve, ncol=3, nrow=2)

##################################### CLUSTERING INDEX ########################################
library(clusterCrit)
library(pdfCluster)
library(aricode)

#### breast
load("/Users/kyle/Downloads/data/fs/breast/truelabel")
View(truelabel)

unique(truelabel$Label)
part1 = truelabel$Label
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
part1 = as.integer(part1)

# FS2000
for (k in 1:length(FS)) {
  fs=FS[k]
  load(paste0("/Users/kyle/Downloads/data/fs/breast/",fs,"2000"))
  if (fs=="var") {
    data=var2000[[1]]
  } else if (fs=="mad") {
    data=mad2000[[1]]
  } else if (fs=="med") {
    data=med2000[[1]]
  } else if (fs=="dip") {
    data=dip2000[[1]]
  }
  wd=paste0("/Users/kyle/Desktop/result_fs/",fs)
  setwd(wd)
  load("result_all.RData")
  
  all.result.fs = analyze.benchmark.fs(result_all) 
  all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  rownames(all.surv.pvalues) = ALGORITHM.NAMES
  colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
  all.clusterings = all.result.fs$all.clusterings
  subtype = colnames(all.surv.pvalues)[3]
  subtype.clusterings = all.clusterings[[subtype]]
  
  extindex = list()
  for (j in 1:length(subtype.clusterings)) {
    extindex[[j]] = list()
    
    clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
    extindex[[j]]$ARI = adj.rand.index(part1, clustering)
    extindex[[j]]$NMI = NMI(part1, clustering)
  }
  save(extindex, file=paste0("/Users/kyle/Desktop/result_fs/",fs,"/breast/extindex")) 
}

# FS500
for (k in 1:length(FS)) {
  fs=FS[k]
  load(paste0("/Users/kyle/Downloads/data/fs/breast/",fs,"500"))
  if (fs=="var") {
    data=var500[[1]]
  } else if (fs=="mad") {
    data=mad500[[1]]
  } else if (fs=="med") {
    data=med500[[1]]
  } else if (fs=="dip") {
    data=dip500[[1]]
  }
  wd=paste0("/Users/kyle/Desktop/result_fs_500/",fs)
  setwd(wd)
  load("result_all.RData")
  
  all.result.fs = analyze.benchmark.fs(result_all) 
  all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  rownames(all.surv.pvalues) = ALGORITHM.NAMES
  colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
  all.clusterings = all.result.fs$all.clusterings
  subtype = colnames(all.surv.pvalues)[3]
  subtype.clusterings = all.clusterings[[subtype]]
  for (j in 1:length(subtype.clusterings)) {
    intindex[[j]] = list()
    extindex[[j]] = list()
    
    clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
    intindex[[j]] = intCriteria(data, clustering, c("Dunn","Silhouette"))
    #extindex[[j]] = extCriteria(part1, clustering, c("Jaccard","Rand"))
    extindex[[j]]$ARI = adj.rand.index(part1, clustering)
    extindex[[j]]$NMI = NMI(part1, clustering)
  }
  save(intindex, file=paste0("/Users/kyle/Desktop/result_fs_500/",fs,"/breast/intindex"))
  save(extindex, file=paste0("/Users/kyle/Desktop/result_fs_500/",fs,"/breast/extindex")) 
}


#### colon
load("/Users/kyle/Downloads/data/fs/colon/truelabel")
unique(truelabel)
part2 = truelabel$Label
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

# FS2000
for (k in 1:length(FS)) {
  fs=FS[k]
  load(paste0("/Users/kyle/Downloads/data/fs/colon/",fs,"2000"))
  if (fs=="var") {
    data=var2000[[1]]
  } else if (fs=="mad") {
    data=mad2000[[1]]
  } else if (fs=="med") {
    data=med2000[[1]]
  } else if (fs=="dip") {
    data=dip2000[[1]]
  }
  wd=paste0("/Users/kyle/Desktop/result_fs/",fs)
  setwd(wd)
  load("result_all.RData")
  
  all.result.fs = analyze.benchmark.fs(result_all) 
  all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  rownames(all.surv.pvalues) = ALGORITHM.NAMES
  colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
  all.clusterings = all.result.fs$all.clusterings
  subtype = colnames(all.surv.pvalues)[4]
  subtype.clusterings = all.clusterings[[subtype]]
  
  intindex = list()
  extindex = list()
  for (j in 1:length(subtype.clusterings)) {
    #intindex[[j]] = list()
    #extindex[[j]] = list()
    
    clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
    intindex[[j]] = intCriteria(data, clustering, "Silhouette")
    extindex[[j]]$ARI = adj.rand.index(part2, clustering)
    extindex[[j]]$NMI = NMI(part2, clustering)
  }
  save(intindex, file=paste0("/Users/kyle/Desktop/result_fs/",fs,"/colon/intindex"))
  save(extindex, file=paste0("/Users/kyle/Desktop/result_fs/",fs,"/colon/extindex")) 
}

# FS500
for (k in 1:length(FS)) {
  fs=FS[k]
  load(paste0("/Users/kyle/Downloads/data/fs/colon/",fs,"500"))
  if (fs=="var") {
    data=var500[[1]]
  } else if (fs=="mad") {
    data=mad500[[1]]
  } else if (fs=="med") {
    data=med500[[1]]
  } else if (fs=="dip") {
    data=dip500[[1]]
  }
  wd=paste0("/Users/kyle/Desktop/result_fs_500/",fs)
  setwd(wd)
  load("result_all.RData")
  
  all.result.fs = analyze.benchmark.fs(result_all) 
  all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  rownames(all.surv.pvalues) = ALGORITHM.NAMES
  colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
  all.clusterings = all.result.fs$all.clusterings
  subtype = colnames(all.surv.pvalues)[4]
  subtype.clusterings = all.clusterings[[subtype]]
  
  extindex = list()
  for (j in 1:length(subtype.clusterings)) {
    extindex[[j]] = list()
    
    clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
    #intindex[[j]] = intCriteria(data, clustering, c("Dunn","Silhouette"))
    extindex[[j]]$ARI = adj.rand.index(part2, clustering)
    extindex[[j]]$NMI = NMI(part2, clustering)
  }
  save(extindex, file=paste0("/Users/kyle/Desktop/result_fs_500/",fs,"/colon/extindex")) 
}

#### Results ######
SUBTYPES.DATA = list(
  list(name='aml', only.primary=F, is.rna.seq=T, is.mirna.seq=T, display.name='AML'),
  list(name='gbm', only.primary=T, is.rna.seq=F, is.mirna.seq=F, display.name='GBM'),
  list(name='breast', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='BIC'),
  list(name='colon', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='COAD'))

ALGORITHM.NAMES = c('cc', 'pins', 'nmf', 'icluster', 'snf', 'nemo')
ALGORITHM.DISPLAY.NAMES = as.list(c('CC', 'PINS', 'NMF', 'iClusterBayes', 'SNF', 'NEMO'))

## without FS
setwd("/home/leejw01/박지윤/result_final")
load("result_all(all).RData")
analyze.benchmark <- function(result_data) {
  all.clusterings = list()
  all.timings = list()
  j=1
  for (i in 1:length(SUBTYPES.DATA)) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    #subtype.raw.data = get.raw.data(subtype,only.primary=current.subtype.data$only.primary)
    
    all.clusterings[[subtype]] = list()
    all.timings[[subtype]] = list()
    
    for (algorithm.name in ALGORITHM.NAMES) {
      all.clusterings[[subtype]][[algorithm.name]] = list()
      all.timings[[subtype]][[algorithm.name]] = list()
      
      clustering = result_data[[j]]$clustering
      timing = result_data[[j]]$timing
      
      #if (!any(is.na(clustering))) {
      #  names(clustering) = colnames(subtype.raw.data[[1]])
      #}
      
      all.clusterings[[subtype]][[algorithm.name]] = clustering
      all.timings[[subtype]][[algorithm.name]] = timing
      j=j+1
    }
  }
  return(list(all.clusterings=all.clusterings, all.timings=all.timings))
}
get.logrank.pvalue <- function(survdiff.res) {
  
  1 - pchisq(survdiff.res$chisq, length(survdiff.res$n) - 1)  
  
}
benchmark.omics.surv <- function(benchmark.results) {
  all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  rownames(all.surv.pvalues) = ALGORITHM.NAMES
  colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
  all.clusterings = benchmark.results$all.clusterings
  for (i in 1:length(all.clusterings)) {
    subtype = colnames(all.surv.pvalues)[i]
    subtype.clusterings = all.clusterings[[subtype]]
    for (j in 1:length(subtype.clusterings)) {
      surv.path = file.path("/home/leejw01/박지윤/result_final/surv", paste(subtype, ALGORITHM.NAMES[j], 'surv', sep='_'))
      load(surv.path)
      pvalue = get.logrank.pvalue(check.survival(clustering, subtype))
      all.surv.pvalues[j, i] = pvalue
    }
  }
  return(all.surv.pvalues)
}
benchmark.omics.time <- function(benchmark.results) {
  
  all.alg.times = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  
  rownames(all.alg.times) = ALGORITHM.NAMES
  
  colnames(all.alg.times) = sapply(SUBTYPES.DATA, function(x) x$name)
  
  all.timings = benchmark.results$all.timings
  
  for (i in 1:length(all.timings)) {
    
    subtype = colnames(all.alg.times)[i]
    
    subtype.timings = all.timings[[subtype]]
    
    for (j in 1:length(subtype.timings)) {
      
      timing = subtype.timings[[ALGORITHM.NAMES[j]]]
      
      all.alg.times[j, i] = timing
      
    }
    
  }
  
  return(all.alg.times)
  
}
benchmark.omics.num.clusters <- function(benchmark.results) {
  
  num.clusters = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  
  rownames(num.clusters) = ALGORITHM.NAMES
  
  colnames(num.clusters) = sapply(SUBTYPES.DATA, function(x) x$name)
  
  all.clusterings = benchmark.results$all.clusterings
  
  for (i in 1:length(all.clusterings)) {
    
    subtype = colnames(num.clusters)[i]
    
    subtype.clusterings = all.clusterings[[subtype]]
    
    for (j in 1:length(subtype.clusterings)) {
      
      clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
      
      num.clusters[j, i] = max(clustering)
      
    }
    
  }
  
  return(num.clusters)
  
}
check.survival <- function(groups, subtype, survival.file.path) {
  
  if (missing(survival.file.path)) {
    
    survival.file.path = get.subtype.survival.path(subtype)
    
  }
  
  survival.data = read.table(survival.file.path, header = TRUE)
  
  patient.names = names(groups)
  
  patient.names.in.file = as.character(survival.data[, 1])
  
  patient.names.in.file = toupper(gsub('-', '\\.', substring(patient.names.in.file, 1, 12)))
  
  
  
  stopifnot(all(patient.names %in% patient.names.in.file))
  
  
  
  indices = match(patient.names, patient.names.in.file)
  
  ordered.survival.data = survival.data[indices,]
  
  ordered.survival.data["cluster"] <- groups
  
  ordered.survival.data$Survival[is.na(ordered.survival.data$Survival)] = 0
  
  ordered.survival.data$Death[is.na(ordered.survival.data$Death)] = 0
  
  return(survdiff(Surv(Survival, Death) ~ cluster, data=ordered.survival.data))
  
}
get.subtype.survival.path <- function(subtype) {
  survival.file.path=file.path("/home/leejw01/박지윤/data",subtype,'survival')
}
benchmark.omics.surv <- function(benchmark.results) {
  all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  rownames(all.surv.pvalues) = ALGORITHM.NAMES
  colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
  all.clusterings = benchmark.results$all.clusterings
  for (i in 1:length(all.clusterings)) {
    subtype = colnames(all.surv.pvalues)[i]
    subtype.clusterings = all.clusterings[[subtype]]
    for (j in 1:length(subtype.clusterings)) {
      surv.path = file.path("/home/leejw01/박지윤/result_final/surv", paste(subtype, ALGORITHM.NAMES[j], 'surv', sep='_'))
      
      if (file.exists(surv.path)) {
        load(surv.path)
        pvalue = empirical.surv.ret$pvalue
      } else {
        clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
        if (length(table(clustering)) > 1) {
          #pvalue = -log10(get.logrank.pvalue(check.survival(clustering, subtype)))
          empirical.surv.ret = get.empirical.surv(clustering, subtype)
          save(empirical.surv.ret, file=surv.path)
          pvalue = empirical.surv.ret$pvalue
        } else {
          pvalue = NA
        }  
      }
      all.surv.pvalues[j, i] = pvalue
    }
  }
  return(all.surv.pvalues)
}

all.result = analyze.benchmark(result_all)
benchmark.omics.surv.fs(all.result)
perform.all.analyses(all.result)
benchmark.sil(result_all,fs="no")

## with FS (log, normalize 안 한 상태에서 돌려야함)
FS=c("var","mad","med","dip")
SUBTYPES.DATA = list(
  list(name='aml', only.primary=F, is.rna.seq=T, is.mirna.seq=T, display.name='AML'),
  list(name='gbm', only.primary=T, is.rna.seq=F, is.mirna.seq=F, display.name='GBM'))
createdata.fs()

SUBTYPES.DATA = list(
  list(name='breast', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='BIC'),
  list(name='colon', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='COAD'))
createdata.fs.gold()

# survival for all data
benchmark.omics.surv.fs <- function(benchmark.results) {
  all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  rownames(all.surv.pvalues) = ALGORITHM.NAMES
  colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
  all.clusterings = benchmark.results$all.clusterings
  for (i in 1:length(all.clusterings)) {
    subtype = colnames(all.surv.pvalues)[i]
    subtype.clusterings = all.clusterings[[subtype]]
    for (j in 1:length(subtype.clusterings)) {
      clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
      pvalue = get.logrank.pvalue(check.survival(clustering, subtype))
      all.surv.pvalues[j, i] = pvalue
    }
  }
  return(all.surv.pvalues)
}
benchmark.omics.surv.fs.re <- function(benchmark.results,fs) {
  all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  rownames(all.surv.pvalues) = ALGORITHM.NAMES
  colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
  all.clusterings = benchmark.results$all.clusterings
  for (i in 1:length(all.clusterings)) {
    subtype = colnames(all.surv.pvalues)[i]
    subtype.clusterings = all.clusterings[[subtype]]
    for (j in 1:length(subtype.clusterings)) {
      surv.path = file.path(surv.fs.dir(), paste(subtype, ALGORITHM.NAMES[j], 'pval', sep='_'))
      print(surv.path)
      if (file.exists(surv.path)) {
        load(surv.path)
        pvalue = pvalue
      } else {
        clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
        if (length(table(clustering)) > 1) {
          pvalue = -log10(get.logrank.pvalue(check.survival(clustering, subtype)))
          empirical.surv.ret = get.empirical.surv(clustering, subtype)
          save(pvalue, file=surv.path)
          #pvalue = empirical.surv.ret$pvalue
        } else {
          pvalue = NA
        }  
      }
      all.surv.pvalues[j, i] = pvalue
    }
  }
  write.csv(all.surv.pvalues, file=paste0("/home/leejw01/박지윤/result_fs/",fs,"/tables/pval.csv"))
}
analyze.benchmark.fs <- function(result_data) {
  all.clusterings = list()
  all.timings = list()
  j=1
  for (i in 1:length(SUBTYPES.DATA)) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    #subtype.raw.data = get.raw.data(subtype,only.primary=current.subtype.data$only.primary)
    
    all.clusterings[[subtype]] = list()
    all.timings[[subtype]] = list()
    
    for (algorithm.name in ALGORITHM.NAMES) {
      all.clusterings[[subtype]][[algorithm.name]] = list()
      all.timings[[subtype]][[algorithm.name]] = list()
      
      clustering = result_data[[j]]$clustering
      timing = result_data[[j]]$timing
      
      all.clusterings[[subtype]][[algorithm.name]] = clustering
      all.timings[[subtype]][[algorithm.name]] = timing 
      j=j+1
    }
  }
  return(list(all.clusterings=all.clusterings, all.timings=all.timings))
}
perform.all.analyses.fs <- function(benchmark.ret, fs) {
  
  for (i in 1:3) {
    cur.func = list(benchmark.omics.time, benchmark.omics.num.clusters,benchmark.omics.surv.fs)[[i]]
    
    benchmark.data = cur.func(benchmark.ret)
    
    displayed.benchmark.data = benchmark.data      
    colnames(displayed.benchmark.data)[1:ncol(displayed.benchmark.data)] = 
      sapply(as.list(colnames(displayed.benchmark.data)[1:ncol(displayed.benchmark.data)]), 
             subtype.to.display.name)
    rownames(displayed.benchmark.data) = ALGORITHM.DISPLAY.NAMES
    print.matrix.latex.format(displayed.benchmark.data)
    
    table.name = c('runtime', 'num_cluster', 'survival')[i]
    
    file.name = deparse(substitute(benchmark.ret))
    file = file.path("/home/leejw01/박지윤/result_fs", paste0(fs,'/tables/',table.name, '.csv'))
    
    write.csv(displayed.benchmark.data, file=file)
    
    print('------------------------')
  }
}
FS=c("var","mad", "med", "dip")
for (k in 1:length(FS)) {
  fs=FS[k]
  wd=paste0("/home/leejw01/박지윤/result_fs/",fs)
  setwd(wd)
  load("result_all.RData")
  
  all.result.fs = analyze.benchmark.fs(result_all) 
  #benchmark.omics.surv.fs(all.result.fs)
  perform.all.analyses.fs(all.result.fs, fs=fs)
  #benchmark.sil(result_all,fs=fs)
}


######## result for fs (학회용 발표) ###################
res = read.csv("/home/leejw01/박지윤/result_fs/fs_res.csv",header=T)
View(res)
res$comb = paste0(res$fs,"_",res$clustering)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggrepel)

res_aml = res[res$data=="AML",]
ggplot(res_aml, aes(x = fs, y=log.pval)) +
  geom_point(aes(colour=clustering)) +
  geom_text_repel(res_aml[which(res_aml$log.pval>1.3),], mapping=aes(label=clustering),hjust=-0.2, vjust=1, size=3)+
  geom_hline(yintercept=2.25)
labs(title = "-log(pval) for AML")
ggplot(res_aml, aes(x = fs, y=silhouette)) +
  geom_point(aes(colour=clustering)) +
  geom_text(res_aml[which(res_aml$silhouette>0.15),], mapping=aes(label=clustering),hjust=-0.2, vjust=1, size=3)+
  labs(title = "silhoutte coef for AML")

res_gbm = res[res$data=="GBM",]
ggplot(res_gbm, aes(x = fs, y=log.pval)) +
  geom_point(aes(colour=clustering)) +
  geom_text(res_gbm[which(res_gbm$log.pval>1.3),], mapping=aes(label=clustering),hjust=-0.2, vjust=1, size=3)+
  labs(title = "-log(pval) for GBM")
ggplot(res_gbm, aes(x = fs, y=silhouette)) +
  geom_point(aes(colour=clustering)) +
  geom_text(res_gbm[which(res_gbm$silhouette>0.15),], mapping=aes(label=clustering),hjust=-0.2, vjust=1, size=3)+
  labs(title = "silhouette for GBM")

res_bic = res[res$data=="BIC",]
ggplot(res_bic, aes(x = fs, y=log.pval)) +
  geom_point(aes(colour=clustering)) +
  geom_text_repel(res_bic[which(res_bic$log.pval>1.3),], mapping=aes(label=clustering),hjust=-0.2, vjust=1, size=3)+
  labs(title = "-log(pval) for BIC")
ggplot(res_bic, aes(x = fs, y=silhouette)) +
  geom_point(aes(colour=clustering)) +
  geom_text(res_bic[which(res_bic$silhouette>0.15),], mapping=aes(label=clustering),hjust=-0.2, vjust=1, size=3)+
  labs(title = "silhouette for BIC")

res_coad = res[res$data=="COAD",]
ggplot(res_coad, aes(x = fs, y=log.pval)) +
  geom_point(aes(colour=clustering)) +
  geom_text(res_coad[which(res_coad$log.pval>1.3),], mapping=aes(label=clustering),hjust=-0.2, vjust=1, size=3)+
  labs(title = "-log(pval) for COAD")
ggplot(res_coad, aes(x = fs, y=silhouette)) +
  geom_point(aes(colour=clustering)) +
  geom_text(res_coad[which(res_coad$silhouette>0.14),], mapping=aes(label=clustering),hjust=-0.2, vjust=1, size=3)+
  labs(title = "silhouette for COAD")

############### KM Curve (w/o FS) ######################
library(ggplot2) 
library(survminer)
check.curve <- function(groups, subtype, survival.file.path) {
  if (missing(survival.file.path)) {
    survival.file.path = get.subtype.survival.path(subtype)
  }
  survival.data = read.table(survival.file.path, header = TRUE)
  patient.names = names(groups)
  patient.names.in.file = as.character(survival.data[, 1])
  patient.names.in.file = toupper(gsub('-', '\\.', substring(patient.names.in.file, 1, 12)))
  
  stopifnot(all(patient.names %in% patient.names.in.file))
  
  indices = match(patient.names, patient.names.in.file)
  ordered.survival.data = survival.data[indices,]
  ordered.survival.data["cluster"] <- groups
  ordered.survival.data$Survival[is.na(ordered.survival.data$Survival)] = 0
  ordered.survival.data$Death[is.na(ordered.survival.data$Death)] = 0
  return(survfit(Surv(Survival, Death) ~ cluster, data=ordered.survival.data))
}

load("/home/leejw01/박지윤/result_final/result_all(all).RData")
all.result = analyze.benchmark.fs(result_all)
all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
rownames(all.surv.pvalues) = ALGORITHM.NAMES
colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
all.clusterings = all.result$all.clusterings
for (i in 1:length(all.clusterings)) {
  subtype = colnames(all.surv.pvalues)[i]
  subtype.clusterings = all.clusterings[[subtype]]
  curve=list()
  for (j in 1:length(subtype.clusterings)) {
    clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
    fit = check.curve(clustering, subtype)
    pvalue = round(get.logrank.pvalue(check.survival(clustering, subtype)), digits = 4)
    curve[[j]] = ggsurvplot(fit, pval=pvalue, pval.coord = c(0, 0.05), legend.title="", font.legend = list(size=10),
                            title=paste0("KM Curve for ",subtype, "_", ALGORITHM.NAMES[j])) +
      guides(colour = guide_legend(nrow = 2))
  }
  arrange_ggsurvplots(curve, ncol=3, nrow=2)
}

benchmark.sil(result_all,fs="no")

######## KM Curve (FS2000) 
check.curve <- function(groups, subtype, survival.file.path) {
  if (missing(survival.file.path)) {
    survival.file.path = get.subtype.survival.path(subtype)
  }
  survival.data = read.table(survival.file.path, header = TRUE)
  patient.names = names(groups)
  patient.names.in.file = as.character(survival.data[, 1])
  patient.names.in.file = toupper(gsub('-', '\\.', substring(patient.names.in.file, 1, 12)))
  
  stopifnot(all(patient.names %in% patient.names.in.file))
  
  indices = match(patient.names, patient.names.in.file)
  ordered.survival.data = survival.data[indices,]
  ordered.survival.data["cluster"] <- groups
  ordered.survival.data$Survival[is.na(ordered.survival.data$Survival)] = 0
  ordered.survival.data$Death[is.na(ordered.survival.data$Death)] = 0
  return(survfit(Surv(Survival, Death) ~ cluster, data=ordered.survival.data))
}

FS=c("var","mad", "med", "dip")

k=4 # change
fs=FS[k]
wd=paste0("/home/leejw01/박지윤/result_fs/",fs)
setwd(wd)
load("result_all.RData")

all.result.fs = analyze.benchmark.fs(result_all) 
all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
rownames(all.surv.pvalues) = ALGORITHM.NAMES
colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
all.clusterings = all.result.fs$all.clusterings

i=1 # change
subtype = colnames(all.surv.pvalues)[i]
subtype.clusterings = all.clusterings[[subtype]]
curve = list()
for (j in 1:length(subtype.clusterings)) {
  clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
  fit = check.curve(clustering, subtype)
  pvalue = round(get.logrank.pvalue(check.survival(clustering, subtype)), digits = 3)
  curve[[j]] = ggsurvplot(fit, pval=pvalue, pval.coord = c(0, 0.05), legend.title="", font.legend = list(size=10),
                          title=paste0("KM Curve for ",subtype, "_", fs, "_", ALGORITHM.NAMES[j])) +
    guides(colour = guide_legend(nrow = 1))
}
arrange_ggsurvplots(curve, ncol=3, nrow=2)


######## KM Curve (FS500)
k=4 # change
fs=FS[k]
wd=paste0("/home/leejw01/박지윤/result_fs_500/",fs)
setwd(wd)
load("result_all.RData")

all.result.fs = analyze.benchmark.fs(result_all) 
all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
rownames(all.surv.pvalues) = ALGORITHM.NAMES
colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
all.clusterings = all.result.fs$all.clusterings

i=4 # change
subtype = colnames(all.surv.pvalues)[i]
subtype.clusterings = all.clusterings[[subtype]]
curve = list()
for (j in 1:length(subtype.clusterings)) {
  clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
  fit = check.curve(clustering, subtype)
  pvalue = round(get.logrank.pvalue(check.survival(clustering, subtype)), digits = 3)
  curve[[j]] = ggsurvplot(fit, pval=pvalue, pval.coord = c(0, 0.05), legend.title="", font.legend = list(size=10),
                          title=paste0("KM Curve for ",subtype, "_", fs, "_", ALGORITHM.NAMES[j])) +
    guides(colour = guide_legend(nrow = 1))
}
arrange_ggsurvplots(curve, ncol=3, nrow=2)

##################################### CLUSTERING INDEX ########################################
library(clusterCrit)
library(pdfCluster)
library(aricode)
FS=c("var","mad", "med", "dip")

#### breast
load("/home/leejw01/박지윤/data/fs/breast/truelabel")
unique(truelabel$Label)
part1 = truelabel$Label
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
part1 = as.integer(part1)

# w/o FS
load("/home/leejw01/박지윤/result_final/result_all(all).RData")
all.result = analyze.benchmark.fs(result_all) 
all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
rownames(all.surv.pvalues) = ALGORITHM.NAMES
colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
all.clusterings = all.result$all.clusterings
subtype = colnames(all.surv.pvalues)[3]
subtype.clusterings = all.clusterings[[subtype]]
extindex = list()
for (j in 1:length(subtype.clusterings)) {
  extindex[[j]] = list()
  clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
  extindex[[j]]$ARI = adj.rand.index(part1, clustering)
  extindex[[j]]$NMI = NMI(part1, clustering)
}
save(extindex, file="/home/leejw01/박지윤/result_final/breast/extindex") 

# FS2000
for (k in 1:length(FS)) {
  fs=FS[k]
  load(paste0("/home/leejw01/박지윤/data/fs/breast/",fs,"2000"))
  if (fs=="var") {
    data=var2000[[1]]
  } else if (fs=="mad") {
    data=mad2000[[1]]
  } else if (fs=="med") {
    data=med2000[[1]]
  } else if (fs=="dip") {
    data=dip2000[[1]]
  }
  wd=paste0("/home/leejw01/박지윤/result_fs/",fs)
  setwd(wd)
  load("result_all.RData")
  
  all.result.fs = analyze.benchmark.fs(result_all) 
  all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  rownames(all.surv.pvalues) = ALGORITHM.NAMES
  colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
  all.clusterings = all.result.fs$all.clusterings
  subtype = colnames(all.surv.pvalues)[3]
  subtype.clusterings = all.clusterings[[subtype]]
  
  extindex = list()
  for (j in 1:length(subtype.clusterings)) {
    extindex[[j]] = list()
    
    clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
    extindex[[j]]$ARI = adj.rand.index(part1, clustering)
    extindex[[j]]$NMI = NMI(part1, clustering)
  }
  save(extindex, file=paste0("/home/leejw01/박지윤/result_fs/",fs,"/breast/extindex")) 
}

# FS500
for (k in 1:length(FS)) {
  fs=FS[k]
  load(paste0("/home/leejw01/박지윤/data/fs/breast/",fs,"500"))
  if (fs=="var") {
    data=var500[[1]]
  } else if (fs=="mad") {
    data=mad500[[1]]
  } else if (fs=="med") {
    data=med500[[1]]
  } else if (fs=="dip") {
    data=dip500[[1]]
  }
  wd=paste0("/home/leejw01/박지윤/result_fs_500/",fs)
  setwd(wd)
  load("result_all.RData")
  
  all.result.fs = analyze.benchmark.fs(result_all) 
  all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  rownames(all.surv.pvalues) = ALGORITHM.NAMES
  colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
  all.clusterings = all.result.fs$all.clusterings
  subtype = colnames(all.surv.pvalues)[3]
  subtype.clusterings = all.clusterings[[subtype]]
  
  intindex = list()
  extindex = list()
  for (j in 1:length(subtype.clusterings)) {
    intindex[[j]] = list()
    extindex[[j]] = list()
    
    clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
    intindex[[j]] = intCriteria(data, clustering, c("Dunn","Silhouette"))
    extindex[[j]]$ARI = adj.rand.index(part1, clustering)
    extindex[[j]]$NMI = NMI(part1, clustering)
  }
  save(intindex, file=paste0("/home/leejw01/박지윤/result_fs_500/",fs,"/breast/intindex"))
  save(extindex, file=paste0("/home/leejw01/박지윤/result_fs_500/",fs,"/breast/extindex")) 
}

#### colon
load("/home/leejw01/박지윤/data/fs/colon/truelabel")
unique(truelabel)
part2 = truelabel$Label
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

# w/o FS
load("/home/leejw01/박지윤/result_final/result_all(all).RData")
all.result = analyze.benchmark.fs(result_all) 
all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
rownames(all.surv.pvalues) = ALGORITHM.NAMES
colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
all.clusterings = all.result$all.clusterings
subtype = colnames(all.surv.pvalues)[4]
subtype.clusterings = all.clusterings[[subtype]]
extindex = list()
for (j in 1:length(subtype.clusterings)) {
  extindex[[j]] = list()
  clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
  extindex[[j]]$ARI = adj.rand.index(part2, clustering)
  extindex[[j]]$NMI = NMI(part2, clustering)
}
save(extindex, file="/home/leejw01/박지윤/result_final/colon/extindex") 

# FS2000
for (k in 1:length(FS)) {
  fs=FS[k]
  load(paste0("/home/leejw01/박지윤/data/fs/colon/",fs,"2000"))
  if (fs=="var") {
    data=var2000[[1]]
  } else if (fs=="mad") {
    data=mad2000[[1]]
  } else if (fs=="med") {
    data=med2000[[1]]
  } else if (fs=="dip") {
    data=dip2000[[1]]
  }
  wd=paste0("/home/leejw01/박지윤/result_fs/",fs)
  setwd(wd)
  load("result_all.RData")
  
  all.result.fs = analyze.benchmark.fs(result_all) 
  all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  rownames(all.surv.pvalues) = ALGORITHM.NAMES
  colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
  all.clusterings = all.result.fs$all.clusterings
  subtype = colnames(all.surv.pvalues)[4]
  subtype.clusterings = all.clusterings[[subtype]]
  
  extindex=list()
  for (j in 1:length(subtype.clusterings)) {
    extindex[[j]] = list()
    
    clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
    extindex[[j]]$ARI = adj.rand.index(part2, clustering)
    extindex[[j]]$NMI = NMI(part2, clustering)
  }
  save(extindex, file=paste0("/home/leejw01/박지윤/result_fs/",fs,"/colon/extindex")) 
}

# FS500
for (k in 1:length(FS)) {
  fs=FS[k]
  load(paste0("/home/leejw01/박지윤/data/fs/colon/",fs,"500"))
  if (fs=="var") {
    data=var500[[1]]
  } else if (fs=="mad") {
    data=mad500[[1]]
  } else if (fs=="med") {
    data=med500[[1]]
  } else if (fs=="dip") {
    data=dip500[[1]]
  }
  wd=paste0("/home/leejw01/박지윤/result_fs_500/",fs)
  setwd(wd)
  load("result_all.RData")
  
  all.result.fs = analyze.benchmark.fs(result_all) 
  all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
  rownames(all.surv.pvalues) = ALGORITHM.NAMES
  colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
  all.clusterings = all.result.fs$all.clusterings
  subtype = colnames(all.surv.pvalues)[4]
  subtype.clusterings = all.clusterings[[subtype]]
  
  extindex = list()
  for (j in 1:length(subtype.clusterings)) {
    extindex[[j]] = list()
    clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
    extindex[[j]]$ARI = adj.rand.index(part2, clustering)
    extindex[[j]]$NMI = NMI(part2, clustering)
  }
  save(extindex, file=paste0("/home/leejw01/박지윤/result_fs_500/",fs,"/colon/extindex")) 
}

load("/home/leejw01/박지윤/result_fs/med/breast/extindex")
load("/home/leejw01/박지윤/result_fs/med/colon/extindex")
View(extindex)

############################# silhouette and dunn
load("/home/leejw01/박지윤/result_final/result_all(all).RData")
all.result = analyze.benchmark.fs(result_all) 
all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
rownames(all.surv.pvalues) = ALGORITHM.NAMES
colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
all.clusterings = all.result$all.clusterings

for (i in 1:length(SUBTYPES.DATA)) {
  current.subtype.data = SUBTYPES.DATA[[i]]
  subtype = current.subtype.data$name
  subtype.clusterings = all.clusterings[[subtype]]
  data = get.raw.data(subtype, only.primary=current.subtype.data$only.primary)
  data = set.omics.list.attr(data, current.subtype.data)
  data = log.and.normalize(data, subtype.data)
  data = as.matrix(data[[1]])
  
  intindex = list()
  for (j in 1:length(subtype.clusterings)) {
    clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
    intindex[[j]] = intCriteria(data, clustering, c("Dunn","Silhouette"))
  }
  save(intindex, file=paste0("/home/leejw01/박지윤/result_final/tables/",subtype,"_intindex"))
}


for (k in 1:length(FS)) {
  fs=FS[k]
  load(paste0("/home/leejw01/박지윤/data/fs/colon/",fs,"500"))
  if (fs=="var") {
    data=var500[[1]]
  } else if (fs=="mad") {
    data=mad500[[1]]
  } else if (fs=="med") {
    data=med500[[1]]
  } else if (fs=="dip") {
    data=dip500[[1]]
  }
  wd=paste0("/home/leejw01/박지윤/result_fs_500/",fs)
  setwd(wd)
  load("result_all.RData")
  
  for (i in 1:length(colnames(all.surv.pvalues))) {
    all.result.fs = analyze.benchmark.fs(result_all) 
    all.surv.pvalues = matrix(1, ncol=length(SUBTYPES.DATA), nrow=length(ALGORITHM.NAMES))
    rownames(all.surv.pvalues) = ALGORITHM.NAMES
    colnames(all.surv.pvalues) = sapply(SUBTYPES.DATA, function(x) x$name)
    all.clusterings = all.result.fs$all.clusterings
    subtype = colnames(all.surv.pvalues)[i]
    subtype.clusterings = all.clusterings[[subtype]]
    
    intindex = list()
    for (j in 1:length(subtype.clusterings)) {
      clustering = subtype.clusterings[[ALGORITHM.NAMES[j]]]
      intindex[[j]] = intCriteria(data, clustering, c("Dunn","Silhouette"))
    } 
  }
  save(intindex, file=paste0("/home/leejw01/박지윤/result_fs_500/",fs,subtype,"/intindex"))
}

load("/home/leejw01/박지윤/result_fs/med/breast/extindex")
load("/home/leejw01/박지윤/result_fs/med/colon/extindex")
View(extindex)

