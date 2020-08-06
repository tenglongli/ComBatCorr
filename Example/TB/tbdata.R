setwd('../')
source("helper.R")
source("Combat_parout.R")
library(sva)
library(invgamma)
library(tidyverse)
library(BiocParallel)
vrnorm <- Vectorize(rnorm,c('mean','sd'))
vinvrgamma <- Vectorize(rinvgamma,c('shape','rate'))


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
  output <- mat.or.vec(gene,15)
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
  output[,1] <- eff
  output[,2] <- f.pvalue(y,mod,mod0) ## Raw Data T-Test
  mod1 <- model.matrix(~xx-1)
  mod10 <- model.matrix(~xx[,-2]-1)
  output[,3] <- f.pvalue(y,mod1,mod10) ## One step approach
  output[,4] <- f.pvalue(y0,mod,mod0) ## Benchmark approach, no batch effect
  cbmod <-ComBatp(y,batch = batchind,mod = mod) ## Power would be reduced if choose a reference batch
  eb <- cbmod$adjdata
  output[,5] <- f.pvalue(eb,mod,mod0) ## ComBat without adjustment
  ## Try the correlation matrix formed by reduce 
  ## Choose delta as 10%
  m=reduce%*%t(reduce)
  v=eigen(m)$vectors
  d=eigen(m)$'values'
  noise=sum(d)*0.1
  d=ifelse(d<0.0001,noise,d) 
  k=v%*%diag(d)%*%t(v)
  s <- solve(t(chol(k)))
  rex <- s%*%x11
  reb <- t(s%*%t(eb))
  mod4 <- model.matrix(~rex-1)
  mod40 <- model.matrix(~rex[,1]-1)
  output[,6] <- f.pvalue(reb,mod4,mod40)
  ## Choose delta as 5%
  d=eigen(m)$'values'
  noise=sum(d)*0.05
  d=ifelse(d<0.0001,noise,d) 
  k=v%*%diag(d)%*%t(v)
  s <- solve(t(chol(k)))
  rex <- s%*%x11
  reb <- t(s%*%t(eb))
  mod4 <- model.matrix(~rex-1)
  mod40 <- model.matrix(~rex[,1]-1)
  output[,7] <- f.pvalue(reb,mod4,mod40)
  ## Choose delta as 2%
  d=eigen(m)$'values'
  noise=sum(d)*0.02
  d=ifelse(d<0.0001,noise,d) 
  k=v%*%diag(d)%*%t(v)
  s <- solve(t(chol(k)))
  rex <- s%*%x11
  reb <- t(s%*%t(eb))
  mod4 <- model.matrix(~rex-1)
  mod40 <- model.matrix(~rex[,1]-1)
  output[,8] <- f.pvalue(reb,mod4,mod40)
  ## Choose delta as 1%
  d=eigen(m)$'values'
  noise=sum(d)*0.01
  d=ifelse(d<0.0001,noise,d) 
  k=v%*%diag(d)%*%t(v)
  s <- solve(t(chol(k)))
  rex <- s%*%x11
  reb <- t(s%*%t(eb))
  mod4 <- model.matrix(~rex-1)
  mod40 <- model.matrix(~rex[,1]-1)
  output[,9] <- f.pvalue(reb,mod4,mod40)
  ## Choose delta as 0.5%
  d=eigen(m)$'values'
  noise=sum(d)*0.005
  d=ifelse(d<0.0001,noise,d) 
  k=v%*%diag(d)%*%t(v)
  s <- solve(t(chol(k)))
  rex <- s%*%x11
  reb <- t(s%*%t(eb))
  mod4 <- model.matrix(~rex-1)
  mod40 <- model.matrix(~rex[,1]-1)
  output[,10] <- f.pvalue(reb,mod4,mod40)
  ## Choose delta as 0.1%
  d=eigen(m)$'values'
  noise=sum(d)*0.001
  d=ifelse(d<0.0001,noise,d) 
  k=v%*%diag(d)%*%t(v)
  s <- solve(t(chol(k)))
  rex <- s%*%x11
  reb <- t(s%*%t(eb))
  mod4 <- model.matrix(~rex-1)
  mod40 <- model.matrix(~rex[,1]-1)
  output[,11] <- f.pvalue(reb,mod4,mod40)
  ## Choose delta as 0.01%
  d=eigen(m)$'values'
  noise=sum(d)*0.0001
  d=ifelse(d<0.0001,noise,d) 
  k=v%*%diag(d)%*%t(v)
  s <- solve(t(chol(k)))
  rex <- s%*%x11
  reb <- t(s%*%t(eb))
  mod4 <- model.matrix(~rex-1)
  mod40 <- model.matrix(~rex[,1]-1)
  output[,12] <- f.pvalue(reb,mod4,mod40)
  ## Choose delta as 0.001%
  d=eigen(m)$'values'
  noise=sum(d)*0.00001
  d=ifelse(d<0.0001,noise,d) 
  k=v%*%diag(d)%*%t(v)
  s <- solve(t(chol(k)))
  rex <- s%*%x11
  reb <- t(s%*%t(eb))
  mod4 <- model.matrix(~rex-1)
  mod40 <- model.matrix(~rex[,1]-1)
  output[,13] <- f.pvalue(reb,mod4,mod40)
  ## Choose delta as 0.0001%
  d=eigen(m)$'values'
  noise=sum(d)*0.000001
  d=ifelse(d<0.0001,noise,d) 
  k=v%*%diag(d)%*%t(v)
  s <- solve(t(chol(k)))
  rex <- s%*%x11
  reb <- t(s%*%t(eb))
  mod4 <- model.matrix(~rex-1)
  mod40 <- model.matrix(~rex[,1]-1)
  output[,14] <- f.pvalue(reb,mod4,mod40)
  ## Choose delta as 0.00001%
  d=eigen(m)$'values'
  noise=sum(d)*0.0000001
  d=ifelse(d<0.0001,noise,d) 
  k=v%*%diag(d)%*%t(v)
  s <- solve(t(chol(k)))
  rex <- s%*%x11
  reb <- t(s%*%t(eb))
  mod4 <- model.matrix(~rex-1)
  mod40 <- model.matrix(~rex[,1]-1)
  output[,15] <- f.pvalue(reb,mod4,mod40)
  colnames(output) <- c('TEff','T','1step','Bench','EB',paste0('EB',c(0.1,0.05,0.02,0.01,0.005,0.001,0.0001,0.00001,0.000001,0.0000001)))
  ## result <- apply(output,2,function(x) length(x[x<0.05])/length(x))
  return(output)
}

simbatch <- compiler::cmpfun(simbatch)

setwd('Data')

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


## Example: TB Batch Data
## Simulation based on the TB data
## out <- simbatch(c(181,399,19),c(77,95,0),intercept = rgamma(24000,0.2,0.00006),mean = c(-0.5,0.25,-0.9),var = c(0.27,0.07,0.69),alpha = c(8,100,1.8),beta = c(4,120,0.4),res = rgamma(24000,0.2,0.0002),eff = c(rep(0,14000),rep(c(-2000,-1000,1000,2000),each=2500)))
out2 <- readRDS('tbbatchsim.rds')
## Process the result based on p values
processp <- function(out){
  out <- data.frame(out)
  result <- matrix(unlist(by(out[,-1],out$TEff,function(x) apply(x,2,function(x) length(x[x<0.05])))),nrow = 5,byrow = TRUE)
  ## Report FPR and TPR based on p values
  ## Report FPR
  output <- mat.or.vec(2,14)
  rownames(output) <- c('FPR','TPR')
  colnames(output) <- colnames(out)[-1]
  output[1,] <- apply(result,2,function(x) round(x[3]/14000,digits = 3))
  ## Report TPR
  output[2,] <- apply(result,2,function(x) round(sum(x[-3])/10000,digits = 3))
  return(output)
}
st <- processp(out2)
out2 <- data.frame(out2)
## Histogram lines plot based on simulation
x2=hist(out2$EB)
y2=hist(out2$Bench)
z2=hist(out2$EB_A3)
plot(x2$breaks[-length(x2$breaks)],x2$counts,type='l',ylab='Frequency',xlab='p-value',main = '(c) Example 4')
lines(y2$breaks[-length(y2$breaks)],y2$counts,lty=2,lwd=1.5)
lines(z2$breaks[-length(z2$breaks)],z2$counts,lty=3,lwd=1.5)
legend('topright',legend = c('Benchmark','ComBat','ComBat3.0'),lty=c(2,1,3),lwd=c(1.5,1.5,1.5),cex = 1,xjust = 0.5,yjust = 0.5)
## Making graph about the value of delta
tpr <- st[2,5:14]
fpr <- st[1,5:14]
delta <-factor(c('10%','5%','2%','1%','0.5%','0.1%','0.01%','0.001%','0.0001%','0.00001%'),levels=c('10%','5%','2%','1%','0.5%','0.1%','0.01%','0.001%','0.0001%','0.00001%'))
dat <- cbind.data.frame(delta,tpr)
dat1 <- cbind.data.frame(delta,fpr)
library(tidyverse)
ggplot(data=dat,aes(x=delta,y=tpr,group=1))+geom_point()+geom_line()+labs(x=expression(delta),y='TPR')
ggplot(data=dat1,aes(x=delta,y=fpr,group=1))+geom_point()+geom_line()+labs(x=expression(delta),y='FPR')

## Read data example
library(SummarizedExperiment)
source('ComBat.R')
source('../helper.R')
dat <- readRDS('tbbatch.rds')
x1 <- assay(dat)
x2 <- data.frame(colData(dat))
s <- which(x2$Label!='Active'&x2$Dataset!='Brazil')
## s <- which(x2$Label!='Active')
x3 <- x2[s,]
x4 <- x1[,s]
b <- as.numeric(factor(x3$SequencingBatch,labels = 1:3))
t <- ifelse(x3$Label=='Progressor',1,0)
mod <- model.matrix(~t)
eb <- ComBat1(x4,batch = b,mod = mod)
data1 <- eb$adjdata
mod0 <- model.matrix(~1,data=x3)
p1 <- as.numeric(f.pvalue(data1,mod,mod0))
p2 <- as.numeric(adj_pval(data1,b,t,0.001))
result <- data.frame(cbind(p1,p2))
## k <- which(is.nan(p1))
## Some p-values are NaN, due to the fact that the expression of those genes are 0 across all the samples
## So we have to drop those genes in the final result also the analysis in the paper
## result <- result[-k,]
length(which(result$p1<0.05&result$p2>0.05))
length(which(result$p1<0.05))
length(which(result$p2<0.05))
plot(result$p1,result$p2)
abline(a=0,b=1,col='red')
## q-values for probe
result1 <- apply(result,2,function(x) p.adjust(x,method = "BH"))
result1 <- data.frame(result1)
colnames(result1) <- c('q1','q2')
length(which(result1$q1<0.05&result1$q2>0.05))
length(which(result1$q1<0.05))
length(which(result1$q2<0.05))




