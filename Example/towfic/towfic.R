# Downloads the non-normalized data from GSE40566
# Using the covariate annotation from GSE61901
loadtowfic = function(downloaddata=TRUE)
{
  geoaccession="GSE40566"
  rawfn = "http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE40566&format=file&file=GSE40566%5Fnon%5Fnormalized%2Etxt%2Egz"	
  
  #createsampleannotation() # ad.hoc function to create sampleannotation.txt based on GSE61901 and GSE40566
  sampleannotation = read.table("sampleannotation.txt", sep="\t",
                                header=TRUE, stringsAsFactors=FALSE)
  
  if(downloaddata)
  {
    temp = tempfile()    
    download.file(url=rawfn, destfile=temp, mode = "wb")
    datamatrix_raw = read.table(temp, sep="\t", header=TRUE, 
                                stringsAsFactors=FALSE, skip=0, strip.white=TRUE, fill=TRUE)
    unlink(temp)
  }else{
    datamatrix_raw = read.table(paste("not_in_github/",geoaccession,"_non-normalized.txt", sep=""), sep="\t", header=TRUE, 
                                stringsAsFactors=FALSE, skip=0, strip.white=TRUE, fill=TRUE)
  }
  
  rownames(datamatrix_raw) =  datamatrix_raw[,1]
  datamatrix_raw = datamatrix_raw[,-1]
  datamatrix_raw = as.matrix(datamatrix_raw)
  
  sampleannotation = sampleannotation[ match(sampleannotation$title, colnames(datamatrix_raw) ) , ]
  return(list(sampleannotation=sampleannotation, data=datamatrix_raw))
}	

# ad hoc function to fetch the covariate labels which are described in GSE61901 but not in GSE40566
# GSE40566 has a "relation" tag that connects to the sample in GSE61901 which has covaraite label.
# I dont know a method of getting the sample annotation from GEO whitout downloading the whole matrix (including the GPL).
createsampleannotation = function()
{
  require("GEOquery")
  GSE61901eset=getGEO("GSE61901")[[1]]
  GSE40566eset=getGEO("GSE40566")[[1]]
  
  #t(pData(GSE61901eset)[1,])
  GSE61901sa = pData(GSE61901eset)[, c("characteristics_ch1", "characteristics_ch1.3", "characteristics_ch1.2")]
  names(GSE61901sa) = c("batch", "covariate", "array_strip_address")
  GSE61901sa[,1]=gsub("batch: ", "", GSE61901sa[,1] )
  GSE61901sa[,2]=gsub("batchcovar: ", "", GSE61901sa[,2] )
  GSE61901sa[,3]=gsub("array_address: ", "", GSE61901sa[,3] )
  GSE61901sa[,"array_hyb_address"]=gsub("_[12]", "", GSE61901sa[,3] )
  # some renaming for conveniance
  
  GSE61901sa$covariate = make.names(GSE61901sa$covariate)
  GSE61901sa$covariate[GSE61901sa$covariate=="GA.DP"] = "DP"
  GSE61901sa$covariate[GSE61901sa$covariate=="GA.Q"] = "N"
  GSE61901sa$covariate[GSE61901sa$covariate=="Medium"] = "M"
  GSE61901sa$covariate[GSE61901sa$covariate=="GA.RS"] = "RS"
  
  #t(pData(GSE40566eset)[1,])
  GSE40566sa  = pData(GSE40566eset)[, c("title", "relation")]
  GSE40566sa$title = make.names(GSE40566sa$title)
  GSE40566sa[,"relation"]=gsub("Reanalyzed by: ", "", GSE40566sa[,"relation"] )
  GSE40566sa = cbind(GSE40566sa, GSE61901sa[GSE40566sa$relation,])
  write.table(GSE40566sa, file="sampleannotation.txt", sep="\t", col.names=NA, quote=FALSE)
}

ret = loadtowfic()
sampleannotation = ret[["sampleannotation"]]
rawdata = ret[["data"]]
library(sva)
library(limma)

qnormdata = normalizeBetweenArrays(rawdata, method="quantile") 


########################################################################
## Start here to read the data instead of download & process the data ##
########################################################################

adj_pval <- function(data,batchind,treat,noise){
  num_batch <- max(batchind)
  size <- length(batchind)
  batch <- mat.or.vec(size,num_batch)
  for (i in 1:num_batch){
    batch[,i] <- ifelse(batchind==i,1,0)
  }
  h1 <- treat%*%solve(t(treat)%*%treat)%*%t(treat)
  h12 <- batch%*%solve(t(batch)%*%(diag(size)-h1)%*%batch)%*%t(batch)%*%(diag(size)-h1)
  reduce <- diag(size)-h12
  m=reduce%*%t(reduce)
  v=eigen(m)$vectors
  d=eigen(m)$'values'
  n=sum(d)*noise
  d=ifelse(d<0.0001,n,d)
  k=v%*%diag(d)%*%t(v)
  s <- solve(t(chol(k)))
  x11 <- as.matrix(cbind(1,treat))
  rex <- s%*%x11
  reb <- t(s%*%t(data))
  mod <- model.matrix(~rex-1)
  mod0 <- model.matrix(~rex[,1]-1)
  adj_p <- f.pvalue(reb,mod,mod0)
  return(adj_p)
}

adj_pval <- compiler::cmpfun(adj_pval)

towfic <- readRDS('towfic.rds')

qnormdata <- towfic[['data']]
sampleannotation <- towfic[['annotation']]

# combat adjust
library(sva)
combatdata <- as.matrix(ComBat(dat=qnormdata,
                             batch=sampleannotation$batch,
                             mod=model.matrix(~as.factor(sampleannotation$covariate))))

## Compare DP vs N 
select <- which(sampleannotation$covariate=='DP'|sampleannotation$covariate=='N')
ebdata <- combatdata[,select]
sample1 <- sampleannotation[select,]
t <- ifelse(sample1$covariate=='DP',1,0)
b <- factor(sample1$batch,labels = 1:17)
b <- as.numeric(b)
mod <- model.matrix(~t)
mod0 <- model.matrix(~1,data = sample1)
p1 <- as.numeric(f.pvalue(ebdata,mod,mod0))
p2 <- as.numeric(adj_pval(ebdata,b,t,0.01))
result <- data.frame(cbind(p1,p2))
length(which(result$p1<0.05&result$p2>0.05))
length(which(result$p1<0.05))
length(which(result$p2<0.05))
plot(result$p1,result$p2)
abline(a=0,b=1,col='red')
## Histogram lines plot based on data
x=hist(result$p1)
y=hist(result$p2)
plot(x$breaks[-length(x$breaks)],x$counts,type='l',ylab='Frequency',xlab='p-value',main = '(a) Example 2')
lines(y$breaks[-length(y$breaks)],y$counts,lty=2,lwd=1.5)
legend('topright',legend = c('ComBat','ComBat3.0'),lty=c(1,2),lwd=c(1.5,1.5),cex = 1,xjust = 0.5,yjust = 0.5)

result1=apply(result,2,function(x) p.adjust(x,method = "BH"))
result1 <- data.frame(result1)
length(which(result1$p1<0.05&result1$p2>0.05))
length(which(result1$p1<0.05))
length(which(result1$p2<0.05))

## Try to do simulation based on this data

#########################################################
## The no-reference-batch version of simbatch function ##
#########################################################

simbatch=function(sample,treated,intercept,mean,var,alpha,beta,res,eff){
  library(sva)
  library(invgamma)
  library(BiocParallel)
  library(Matrix)
  batch <- length(sample)
  size <- sum(sample)
  gene <- length(res)
  output <- mat.or.vec(gene,11)
  k <- cumsum(sample)
  x2 <- mat.or.vec(size,batch-1)
  x1 <- rep(c(0,1),times=c(sample[1]-treated[1],treated[1])) ## Treatment indicator
  for (i in 1:(batch-1)){
    x1 <- c(x1,rep(c(0,1),times=c(sample[i+1]-treated[i+1],treated[i+1])))
    x2[(k[i]+1):k[i+1],i] <- 1
  }
  x22 <- as.matrix(cbind(rep(c(1,0),times=c(sample[1],size-sample[1])),x2)) ## For simulation only
  x11 <- as.matrix(cbind(1,x1))
  xx <- as.matrix(cbind(x11,x2))
  x <- as.matrix(cbind(x11,x22))
  h1 <- x1%*%solve(t(x1)%*%x1)%*%t(x1)
  h12 <- x22%*%solve(t(x22)%*%(diag(size)-h1)%*%x22)%*%t(x22)%*%(diag(size)-h1)
  reduce <- diag(size)-h12
  batchind <- rep(1:batch,times=sample)
  mod <- model.matrix(~x1)
  mod0 <- model.matrix(~1,data=data.frame(x))
  varbatch <- t(vinvrgamma(gene,shape=alpha,rate=beta))[rep(1:batch,times=sample),]
  resid <- vrnorm(size,mean = 0, sd = res)
  meanbatch <- t(vrnorm(gene,mean,var^0.5))
  null <- mat.or.vec(batch,gene)
  y0 <- t(x%*%rbind(intercept,eff,null)+resid)
  y <- t(x%*%rbind(intercept,eff,meanbatch)+resid*varbatch)
  cy <- t(reduce%*%t(y))
  output[,1] <- eff
  output[,2] <- f.pvalue(y,mod,mod0) ## Raw Data T-Test
  mod1 <- model.matrix(~xx-1)
  mod10 <- model.matrix(~xx[,-2]-1)
  output[,3] <- f.pvalue(y,mod1,mod10) ## One step approach
  output[,4] <- f.pvalue(y0,mod,mod0) ## Benchmark approach, no batch effect
  output[,5] <- f.pvalue(cy,mod,mod0) ## Two step approach without adjustment
  m=reduce%*%t(reduce)
  v=eigen(m)$vectors
  d=eigen(m)$'values'
  d=ifelse(d<0.01,0.01,d) ## t is the threshold, typically chosen as 0.01
  k=v%*%diag(d)%*%t(v)
  s <- solve(t(chol(k)))
  ecy <- t(s%*%t(cy))
  rx11 <- s%*%x11
  mod2 <- model.matrix(~rx11-1)
  mod20 <- model.matrix(~rx11[,1]-1)
  output[,6] <- f.pvalue(ecy,mod2,mod20) ## Two step approach with adjustment 
  cbmod <-ComBatp(y,batch = batchind,mod = mod) ## Power would be reduced if choose a reference batch
  eb <- cbmod$adjdata
  output[,7] <- f.pvalue(eb,mod,mod0) ## ComBat without adjustment
  ## Try the correlation matrix formed by reduce 
  m=reduce%*%t(reduce)
  v=eigen(m)$vectors
  d=eigen(m)$'values'
  noise=sum(d)*0.1
  d=ifelse(d<0.0001,noise,d) ## t is the threshold, typically chosen as 0.1
  k=v%*%diag(d)%*%t(v)
  s <- solve(t(chol(k)))
  rex <- s%*%x11
  reb <- t(s%*%t(eb))
  mod4 <- model.matrix(~rex-1)
  mod40 <- model.matrix(~rex[,1]-1)
  output[,8] <- f.pvalue(reb,mod4,mod40)
  ## Try the correlation matrix formed by reduce
  d=eigen(m)$'values'
  noise=sum(d)*0.01
  d=ifelse(d<0.0001,noise,d) ## t is the threshold, typically chosen as 0.01
  k=v%*%diag(d)%*%t(v)
  s <- solve(t(chol(k)))
  rex <- s%*%x11
  reb <- t(s%*%t(eb))
  mod4 <- model.matrix(~rex-1)
  mod40 <- model.matrix(~rex[,1]-1)
  output[,9] <- f.pvalue(reb,mod4,mod40)
  ## Try the correlation matrix formed by reduce
  d=eigen(m)$'values'
  noise=sum(d)*0.001
  d=ifelse(d<0.0001,noise,d) ## t is the threshold, typically chosen as 0.001
  k=v%*%diag(d)%*%t(v)
  s <- solve(t(chol(k)))
  rex <- s%*%x11
  reb <- t(s%*%t(eb))
  mod4 <- model.matrix(~rex-1)
  mod40 <- model.matrix(~rex[,1]-1)
  output[,10] <- f.pvalue(reb,mod4,mod40)
  ## Try the correlation matrix formed by reduce
  d=eigen(m)$'values'
  noise=sum(d)*0.00001
  d=ifelse(d<0.0001,noise,d) ## t is the threshold, typically chosen as 0.1
  k=v%*%diag(d)%*%t(v)
  s <- solve(t(chol(k)))
  rex <- s%*%x11
  reb <- t(s%*%t(eb))
  mod4 <- model.matrix(~rex-1)
  mod40 <- model.matrix(~rex[,1]-1)
  output[,11] <- f.pvalue(reb,mod4,mod40)
  colnames(output) <- c('TEff','T','1step','Bench','2step','2step_A','EB','EB_A1','EB_A2','EB_A3','EB_A4')
  ## result <- apply(output,2,function(x) length(x[x<0.05])/length(x))
  return(output)
}

simbatch <- compiler::cmpfun(simbatch)

setwd('../')
source("helper.R")
source("Combat_parout.R")
library(sva)
library(invgamma)
library(tidyverse)
library(BiocParallel)
vrnorm <- Vectorize(rnorm,c('mean','sd'))
vinvrgamma <- Vectorize(rinvgamma,c('shape','rate'))
qnormdata1 <- qnormdata[,select]
eb <- ComBatp(dat=qnormdata1, batch=sample1$batch, mod=model.matrix(~as.factor(sample1$covariate)))
batch <- as.numeric(table(b))
control <- numeric(17)
control[c(2,3,6,7,8,9,11,12)] <- c(2,1,1,1,1,1,2,2)
treat <- batch-control
out <- simbatch(batch,treat,intercept =5+rgamma(46547,3,1.8),mean = rep(0,17),var = c(4,16,9,16,9,4,8,10,12,18,16,14,20,20,15,18,25),alpha = rep(150,17),beta = rep(150,17),res = rgamma(46547,1,12),eff = c(rep(0,46447),rep(c(-1,1),each=50)))
out <- data.frame(out)
## Histogram lines plot based on simulation
x=hist(out$EB)
y=hist(out$Bench)
z=hist(out$EB_A2)
plot(x$breaks[-length(x$breaks)],x$counts,type='l',ylab='Frequency',xlab='p-value',main = '(a) Example 2')
lines(y$breaks[-length(y$breaks)],y$counts,lty=2,lwd=1.5)
lines(z$breaks[-length(z$breaks)],z$counts,lty=3,lwd=1.5)
legend('topright',legend = c('Benchmark','ComBat','ComBat3.0'),lty=c(2,1,3),lwd=c(1.5,1.5,1.5),cex = 1,xjust = 0.5,yjust = 0.5)







