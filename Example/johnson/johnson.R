loadjohnsondata = function()
{  
  datamatrix = as.matrix(read.table("dataExample2.txt", sep="\t", header=TRUE))
  sampleannotation = read.table("sampleInfoExample2.txt", 
                                sep="\t", header=TRUE, stringsAsFactors=FALSE)
  rownames(sampleannotation)=sampleannotation$ArrayName
  sampleannotation$Batch=factor(as.character(sampleannotation$Batch)) 
  datamatrix=datamatrix[,sampleannotation$Type!="WT"]
  sampleannotation=sampleannotation[sampleannotation$Type!="WT",]
  sampleannotation$Type=factor(sampleannotation$Type)
  
  log2data = datamatrix
  # flooring to 1
  log2data[log2data<1]=1
  # take out data with to much low/missing values.
  negativeprobesfilter =( rowSums(log2data>1) >= (0.9*ncol(log2data)) )
  log2data = log2data[negativeprobesfilter,]
  # quantilenormalize
  log2data=normalizeBetweenArrays(log2(log2data), method="quantile")
  
  return(list(sampleannotation=sampleannotation, rawdata=datamatrix, normdata=log2data))
}

##########################################################
## Make the function to do the adjustment for real data ##
##########################################################

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

library(sva)
library(limma)
tmp = loadjohnsondata()
sampleannotation = tmp[["sampleannotation"]]
normdata = tmp[["normdata"]]
combatdata = as.matrix(ComBat(
  dat=normdata, 
  batch=sampleannotation$Batch, 
  mod=model.matrix(~as.factor(sampleannotation$"Type"))))
t <- ifelse(sampleannotation$Type=='R',1,0)
b <- as.numeric(sampleannotation$Batch)
mod <- model.matrix(~t)
mod0 <- model.matrix(~1,data = sampleannotation)
p1 <- as.numeric(f.pvalue(combatdata,mod,mod0))
p2 <- as.numeric(adj_pval(combatdata,b,t,0.01))
result <- data.frame(cbind(p1,p2))
length(which(result$p1<0.05&result$p2>0.05))
length(which(result$p1<0.05))
length(which(result$p2<0.05))
plot(result$p1,result$p2)
abline(a=0,b=1,col='red')
## Histogram lines plot based on data
x=hist(result$p1)
y=hist(result$p2)
plot(x$breaks[-length(x$breaks)],x$counts,type='l',ylab='Frequency',xlab='p-value',main = '(b) Example 3')
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
## eb <- ComBatp(dat=normdata, batch=b, mod=model.matrix(~t))

## Only report the FPR and TPR based on p values
out1 <- simbatch(c(8,7,15),c(6,3,9),intercept = rgamma(52275,12,1.8),mean = rep(0,3),var = c(0.5,0.3,0.15),alpha = rep(100,3),beta = rep(100,3),res = rgamma(52275,0.4,2),eff = c(rep(0,48275),rep(c(-2,-1,1,2),each=1000)))
out1 <- data.frame(out1)
by(out1[,-1],out1$TEff,function(x) apply(x,2,function(x) length(x[x<0.05])))

## Histogram lines plot
x1=hist(out1$EB)
y1=hist(out1$Bench)
z1=hist(out1$EB_A2)
plot(x1$breaks[-length(x1$breaks)],x1$counts,type='l',ylab='Frequency',xlab='p-value',main = '(b) Example 3')
lines(y1$breaks[-length(y1$breaks)],y1$counts,lty=2,lwd=1.5)
lines(z1$breaks[-length(z1$breaks)],z1$counts,lty=3,lwd=1.5)
legend('topright',legend = c('Benchmark','ComBat','ComBat3.0'),lty=c(2,1,3),lwd=c(1.5,1.5,1.5),cex = 1,xjust = 0.5,yjust = 0.5)


## Likewise, only report FDR and TPR based on q values
out11 <- apply(out1[,-1],2,function(x) p.adjust(x,method = "BH"))
out11 <- cbind(out1[,1],out11)
out11 <- data.frame(out11)
by(out11[,-1],out11$V1,function(x) apply(x,2,function(x) length(x[x<0.05])))

