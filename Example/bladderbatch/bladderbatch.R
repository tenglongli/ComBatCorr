## Try to use the following function to calculate the correlation among residuals
## Also need to source ComBat function that outputs the parameters
source("helper.R")
source("Combat_parout.R")
library(sva)
library(invgamma)
library(tidyverse)
library(BiocParallel)
vrnorm <- Vectorize(rnorm,c('mean','sd'))
vinvrgamma <- Vectorize(rinvgamma,c('shape','rate'))

##############################################################################################################
## A simple version with gene-specific mean and variance effects
## Will simulate those effects within the function
## The mean effect will be generated from normal distribution with mean and var specified for each batch
## The variance effect will be generated from inverse-gamma distribution with alpha and beta
## specified for each batch
## res is the vector of residual variances (one for each gene)
## eff is the known treatment effect


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
  output <- mat.or.vec(gene,13)
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
  colnames(output) <- c('TEff','T','1step','Bench','EB',paste0('EB',c(0.1,0.05,0.02,0.01,0.005,0.001,0.0001,0.00001)))
  ## result <- apply(output,2,function(x) length(x[x<0.05])/length(x))
  return(output)
}

simbatch <- compiler::cmpfun(simbatch)


##################################################
## Simple simulation for mean batch effect only ##
##################################################

simbatch1=function(sample,treated,intercept,mean,var,res,eff){
  library(sva)
  library(BiocParallel)
  library(Matrix)
  batch <- length(sample)
  size <- sum(sample)
  gene <- length(res)
  output <- mat.or.vec(gene,13)
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
  resid <- vrnorm(size,mean = 0, sd = res)
  meanbatch <- t(vrnorm(gene,mean,var^0.5))
  null <- mat.or.vec(batch,gene)
  y0 <- t(x%*%rbind(intercept,eff,null)+resid)
  y <- t(x%*%rbind(intercept,eff,meanbatch)+resid)
  output[,1] <- eff
  output[,2] <- f.pvalue(y,mod,mod0) ## Raw Data T-Test
  mod1 <- model.matrix(~xx-1)
  mod10 <- model.matrix(~xx[,-2]-1)
  output[,3] <- f.pvalue(y,mod1,mod10) ## One step approach
  output[,4] <- f.pvalue(y0,mod,mod0) ## Benchmark approach, no batch effect
  cbmod <-ComBatp(y,batch = batchind,mod = mod,mean.only = TRUE) ## Power would be reduced if choose a reference batch
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
  colnames(output) <- c('TEff','T','1step','Bench','EB',paste0('EB',c(0.1,0.05,0.02,0.01,0.005,0.001,0.0001,0.00001)))
  ## result <- apply(output,2,function(x) length(x[x<0.05])/length(x))
  return(output)
}

simbatch1 <- compiler::cmpfun(simbatch1)

##########################################################
## The simple simulation for variance batch effect only ##
##########################################################

simbatch2=function(sample,treated,intercept,alpha,beta,res,eff){
  library(sva)
  library(invgamma)
  library(BiocParallel)
  library(Matrix)
  batch <- length(sample)
  size <- sum(sample)
  gene <- length(res)
  output <- mat.or.vec(gene,13)
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
  null <- mat.or.vec(batch,gene)
  y0 <- t(x%*%rbind(intercept,eff,null)+resid)
  y <- t(x%*%rbind(intercept,eff,null)+resid*varbatch)
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
  colnames(output) <- c('TEff','T','1step','Bench','EB',paste0('EB',c(0.1,0.05,0.02,0.01,0.005,0.001,0.0001,0.00001)))
  ## result <- apply(output,2,function(x) length(x[x<0.05])/length(x))
  return(output)
}

simbatch2 <- compiler::cmpfun(simbatch2)

## Simulation design: 3*3
unbalanced <- vector("list",9)
## For unbalanced design
## Just mean batch effect only
## No mean batch effect & No variance batch effect
out <- simbatch1(c(11,18,4,5,19),c(11,14,0,0,15),intercept=3+rgamma(20000,4.5,1.5),mean=rep(0,5),var=rep(0,5),res=rgamma(20000,4,10),eff=c(rep(0,18000),rep(c(-2,-1,1,2),each=500)))
unbalanced[[1]] <- out
## Small mean batch effect & No variance batch effect
out <- simbatch1(c(11,18,4,5,19),c(11,14,0,0,15),intercept=3+rgamma(20000,4.5,1.5),mean=c(-0.04,0.15,-0.15,-0.1,-0.08),var=c(0.15,0.35,0.82,0.46,0.12),res=rgamma(20000,4,10),eff=c(rep(0,18000),rep(c(-2,-1,1,2),each=500)))
unbalanced[[2]] <- out
## Large mean batch effect & No variance batch effect
out <- simbatch1(c(11,18,4,5,19),c(11,14,0,0,15),intercept=3+rgamma(20000,4.5,1.5),mean=c(-0.4,1.5,-1.5,-1,-0.8),var=c(0.15,0.35,0.82,0.46,0.12),res=rgamma(20000,4,10),eff=c(rep(0,18000),rep(c(-2,-1,1,2),each=500)))
unbalanced[[3]] <- out
## Just variance batch effect only
## No mean batch effect & Small variance batch effect 
out <- simbatch2(c(11,18,4,5,19),c(11,14,0,0,15),intercept=3+rgamma(20000,4.5,1.5),alpha=c(60,100,56,30,100),beta=c(60,100,50,30,100),res=rgamma(20000,4,10),eff=c(rep(0,18000),rep(c(-2,-1,1,2),each=500)))
unbalanced[[4]] <- out
## No mean batch effect & Large variance batch effect
out <- simbatch2(c(11,18,4,5,19),c(11,14,0,0,15),intercept=3+rgamma(20000,4.5,1.5),alpha=c(100,120,100,60,40),beta=c(100,40,60,100,120),res=rgamma(20000,4,10),eff=c(rep(0,18000),rep(c(-2,-1,1,2),each=500)))
unbalanced[[5]] <- out
## Have both mean batch effect and variance batch effect
## Small mean batch effect & Small variance batch effect
out <- simbatch(c(11,18,4,5,19),c(11,14,0,0,15),intercept=3+rgamma(20000,4.5,1.5),mean=c(-0.04,0.15,-0.15,-0.1,-0.08),var=c(0.15,0.35,0.82,0.46,0.12),alpha=c(60,100,56,30,100),beta=c(60,100,50,30,100),res=rgamma(20000,4,10),eff=c(rep(0,18000),rep(c(-2,-1,1,2),each=500)))
unbalanced[[6]] <- out
## Large mean batch effect & Small variance batch effect
out <- simbatch(c(11,18,4,5,19),c(11,14,0,0,15),intercept=3+rgamma(20000,4.5,1.5),mean=c(-0.4,1.5,-1.5,-1,-0.8),var=c(0.15,0.35,0.82,0.46,0.12),alpha=c(60,100,56,30,100),beta=c(60,100,50,30,100),res=rgamma(20000,4,10),eff=c(rep(0,18000),rep(c(-2,-1,1,2),each=500)))
unbalanced[[7]] <- out
## Small mean batch effect & Large variance batch effect
out <- simbatch(c(11,18,4,5,19),c(11,14,0,0,15),intercept=3+rgamma(20000,4.5,1.5),mean=c(-0.04,0.15,-0.15,-0.1,-0.08),var=c(0.15,0.35,0.82,0.46,0.12),alpha=c(100,120,100,60,40),beta=c(100,40,60,100,120),res=rgamma(20000,4,10),eff=c(rep(0,18000),rep(c(-2,-1,1,2),each=500)))
unbalanced[[8]] <- out
## Large mean batch effect & Large variance batch effect 
out <- simbatch(c(11,18,4,5,19),c(11,14,0,0,15),intercept=3+rgamma(20000,4.5,1.5),mean=c(-0.4,1.5,-1.5,-1,-0.8),var=c(0.15,0.35,0.82,0.46,0.12),alpha=c(100,120,100,60,40),beta=c(100,40,60,100,120),res=rgamma(20000,4,10),eff=c(rep(0,18000),rep(c(-2,-1,1,2),each=500)))
unbalanced[[9]] <- out
saveRDS(unbalanced,'unbalanced.rds')

balanced <- vector("list",9)
## For balanced design
## Just mean batch effect only
## No mean batch effect & No variance batch effect
out <- simbatch1(c(12,18,4,6,20),c(6,9,2,3,10),intercept=3+rgamma(20000,4.5,1.5),mean=rep(0,5),var=rep(0,5),res=rgamma(20000,4,10),eff=c(rep(0,18000),rep(c(-2,-1,1,2),each=500)))
balanced[[1]] <- out
## Small mean batch effect & No variance batch effect
out <- simbatch1(c(12,18,4,6,20),c(6,9,2,3,10),intercept=3+rgamma(20000,4.5,1.5),mean=c(-0.04,0.15,-0.15,-0.1,-0.08),var=c(0.15,0.35,0.82,0.46,0.12),res=rgamma(20000,4,10),eff=c(rep(0,18000),rep(c(-2,-1,1,2),each=500)))
balanced[[2]] <- out
## Large mean batch effect & No variance batch effect
out <- simbatch1(c(12,18,4,6,20),c(6,9,2,3,10),intercept=3+rgamma(20000,4.5,1.5),mean=c(-0.4,1.5,-1.5,-1,-0.8),var=c(0.15,0.35,0.82,0.46,0.12),res=rgamma(20000,4,10),eff=c(rep(0,18000),rep(c(-2,-1,1,2),each=500)))
balanced[[3]] <- out
## Just variance batch effect only
## No mean batch effect & Small variance batch effect 
out <- simbatch2(c(12,18,4,6,20),c(6,9,2,3,10),intercept=3+rgamma(20000,4.5,1.5),alpha=c(60,100,56,30,100),beta=c(60,100,50,30,100),res=rgamma(20000,4,10),eff=c(rep(0,18000),rep(c(-2,-1,1,2),each=500)))
balanced[[4]] <- out
## No mean batch effect & Large variance batch effect
out <- simbatch2(c(12,18,4,6,20),c(6,9,2,3,10),intercept=3+rgamma(20000,4.5,1.5),alpha=c(100,120,100,60,40),beta=c(100,40,60,100,120),res=rgamma(20000,4,10),eff=c(rep(0,18000),rep(c(-2,-1,1,2),each=500)))
balanced[[5]] <- out
## Have both mean batch effect and variance batch effect
## Small mean batch effect & Small variance batch effect
out <- simbatch(c(12,18,4,6,20),c(6,9,2,3,10),intercept=3+rgamma(20000,4.5,1.5),mean=c(-0.04,0.15,-0.15,-0.1,-0.08),var=c(0.15,0.35,0.82,0.46,0.12),alpha=c(60,100,56,30,100),beta=c(60,100,50,30,100),res=rgamma(20000,4,10),eff=c(rep(0,18000),rep(c(-2,-1,1,2),each=500)))
balanced[[6]] <- out
## Large mean batch effect & Small variance batch effect
out <- simbatch(c(12,18,4,6,20),c(6,9,2,3,10),intercept=3+rgamma(20000,4.5,1.5),mean=c(-0.4,1.5,-1.5,-1,-0.8),var=c(0.15,0.35,0.82,0.46,0.12),alpha=c(60,100,56,30,100),beta=c(60,100,50,30,100),res=rgamma(20000,4,10),eff=c(rep(0,18000),rep(c(-2,-1,1,2),each=500)))
balanced[[7]] <- out
## Small mean batch effect & Large variance batch effect
out <- simbatch(c(12,18,4,6,20),c(6,9,2,3,10),intercept=3+rgamma(20000,4.5,1.5),mean=c(-0.04,0.15,-0.15,-0.1,-0.08),var=c(0.15,0.35,0.82,0.46,0.12),alpha=c(100,120,100,60,40),beta=c(100,40,60,100,120),res=rgamma(20000,4,10),eff=c(rep(0,18000),rep(c(-2,-1,1,2),each=500)))
balanced[[8]] <- out
## Large mean batch effect & Large variance batch effect 
out <- simbatch(c(12,18,4,6,20),c(6,9,2,3,10),intercept=3+rgamma(20000,4.5,1.5),mean=c(-0.4,1.5,-1.5,-1,-0.8),var=c(0.15,0.35,0.82,0.46,0.12),alpha=c(100,120,100,60,40),beta=c(100,40,60,100,120),res=rgamma(20000,4,10),eff=c(rep(0,18000),rep(c(-2,-1,1,2),each=500)))
balanced[[9]] <- out
saveRDS(balanced,'balanced.rds')

## Postprocessing
## Process the result based on p values
processp <- function(out){
  out <- data.frame(out)
  result <- matrix(unlist(by(out[,-1],out$TEff,function(x) apply(x,2,function(x) length(x[x<0.05])))),nrow = 5,byrow = TRUE)
  ## Report FPR and TPR based on p values
  ## Report FPR
  output <- mat.or.vec(2,12)
  rownames(output) <- c('FPR','TPR')
  colnames(output) <- colnames(out)[-1]
  output[1,] <- apply(result,2,function(x) round(x[3]/18000,digits = 3))
  ## Report TPR
  output[2,] <- apply(result,2,function(x) round(sum(x[-3])/2000,digits = 3))
  return(output)
}
## Process the result based on q values
processq <- function(out){
  output <- mat.or.vec(2,12)
  rownames(output) <- c('FDR','TPR')
  colnames(output) <- colnames(out)[-1]
  out1 <- apply(out[,-1],2,function(x) p.adjust(x,method = "BH"))
  out1 <- cbind(out[,1],out1)
  out1 <- data.frame(out1)
  result <- matrix(unlist(by(out1[,-1],out1[,1],function(x) apply(x,2,function(x) length(x[x<0.05])))),nrow = 5,byrow = TRUE)
  ## Report FDR
  output[1,] <- apply(result,2,function(x) x[3]/sum(x))
  ## Report TPR
  output[2,] <- apply(result,2,function(x) sum(x[-3])/2000)
  return(output)
}


## Read and process the output
setwd('Output/')
unbalanced <- readRDS('unbalanced.rds')
balanced <- readRDS('balanced.rds')
u <- lapply(unbalanced,processp)
b <- lapply(balanced,processp)
## Making graph about the value of delta
tpr_u <- apply(do.call("rbind",lapply(u,function(x) x[2,5:12])),2,mean)
fpr_u <- apply(do.call("rbind",lapply(u,function(x) x[1,5:12])),2,mean)
tpr_b <- apply(do.call("rbind",lapply(b,function(x) x[2,5:12])),2,mean)
fpr_b <- apply(do.call("rbind",lapply(b,function(x) x[1,5:12])),2,mean)
delta <-factor(c('10%','5%','2%','1%','0.5%','0.1%','0.01%','0.001%'),levels=c('10%','5%','2%','1%','0.5%','0.1%','0.01%','0.001%'))
dat <- cbind.data.frame(delta,unbalanced=tpr_u,balanced=tpr_b)
dat1 <- cbind.data.frame(delta,unbalanced=fpr_u,balanced=fpr_b)
library(tidyverse)
dat <- dat%>%gather(key = 'design',value = 'value',-delta)
dat1 <- dat1%>%gather(key = 'design',value = 'value',-delta)
ggplot(data=dat,aes(x=delta,y=value,group=design))+geom_point(aes(shape=design))+geom_line(aes(linetype=design))+labs(x=expression(delta),y='TPR')
ggplot(data=dat1,aes(x=delta,y=value,group=design))+geom_point(aes(shape=design))+geom_line(aes(linetype=design))+labs(x=expression(delta),y='FPR')



## Generate Line-Histogram Plot for small mean and variance batch effect
out <- simbatch(c(11,18,4,5,19),c(11,14,0,0,15),intercept=3+rgamma(20000,4.5,1.5),mean=c(-0.04,0.15,-0.15,-0.1,-0.08),var=c(0.15,0.35,0.82,0.46,0.12),alpha=c(60,100,56,30,100),beta=c(60,100,50,30,100),res=rgamma(20000,4,10),eff=c(rep(0,18000),rep(c(-2,-1,1,2),each=500)))
out <- data.frame(out)
x=hist(out$EB)
y=hist(out$Bench)
z=hist(out$EB0.01)
## QQ Plot with reference line at Benchmark
par(mfrow=c(1,3))
qqplot(out$EB,out$Bench,xlab = "ComBat", ylab = "Benchmark",main = '(a) QQ plot of ComBat versus Benchmark')
abline(a=0,b=1,lty=2)
qqplot(out$EB0.01,out$Bench,xlab="ComBat3.0",ylab='Benchmark',main = '(b) QQ plot of ComBat3.0 versus Benchmark')
abline(a=0,b=1,lty=2)
plot(x$breaks[-length(x$breaks)],x$counts,type='l',ylab='Frequency',xlab='p-value',main = '(c) Distributions of p-values for the three approaches')
lines(y$breaks[-length(y$breaks)],y$counts,lty=2,lwd=1.5)
lines(z$breaks[-length(z$breaks)],z$counts,lty=3,lwd=1.5)
legend('topright',legend = c('Benchmark','ComBat','ComBat3.0'),lty=c(2,1,3),lwd=c(1.5,1.5,1.5),cex = 1,xjust = 0.5,yjust = 0.5)
par(mfrow=c(1,1))







