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

## Simple examples
## 2 batches: each batch has 8 samples, each sample consists of 10 gene expressions
## group-batch design (treated/control): batch 1: 1/7 and batch 2: 8/0
## Example 1: Show ComBat will exaggerate significance for null genes
## Treatment effect is fixed at 0
par(mfrow=c(1,2))
gene <- 10000
batch.design <- cbind(rep(c(1,0),each=8),rep(c(0,1),each=8))
group.design <- c(rep(c(1,0),times=c(1,7)),rep(c(1,0),times=c(8,0)))
background <- rnorm(gene)
meanbatch <- rbind(rep(1,gene),rep(4,gene))
varbatch <- matrix(1,16,gene)
resid <- vrnorm(16,mean = rep(0,gene),sd=rep(1,gene)) 
design <- cbind(1,group.design,batch.design)
effect <- rbind(background,0,meanbatch)
null <- rbind(background,0,mat.or.vec(2,gene))
y0 <- t(design%*%null+resid)
y <- t(design%*%effect+resid)
batch.ind <- rep(c(1,2),each=8)
mod <- model.matrix(~as.factor(group.design))
mod0 <- model.matrix(~1,data=data.frame(design))
y1 <- ComBat(y,batch = batch.ind,mod = mod,mean.only = TRUE)
p1 <- f.pvalue(y1,mod,mod0)
p0 <- f.pvalue(y0,mod,mod0)
qqplot(p1,p0,xlab = "ComBat", ylab = "Benchmark",main = '(a) No group effect and variance batch effect \n strong mean batch effect')
abline(a=0,b=1,lty=2)
## Use boxplots to show the comparison of groups
## boxplot(p1,p0)

## Example 2: Show ComBat will lose power for differentially expressed genes
## Treatment effect fixed at 1
gene <- 10000
batch.design <- cbind(rep(c(1,0),each=8),rep(c(0,1),each=8))
group.design <- c(rep(c(1,0),times=c(1,7)),rep(c(1,0),times=c(8,0)))
background <- rnorm(gene)
meanbatch <- rbind(rep(1,gene),rep(4,gene))
varbatch <- rbind(rep(0.5,gene),rep(4,gene))[rep(1:2,each=8),]
resid <- vrnorm(16,mean = rep(0,gene),sd=rep(1,gene)) 
design <- cbind(1,group.design,batch.design)
effect <- rbind(background,1,meanbatch)
null <- rbind(background,1,mat.or.vec(2,gene))
y0 <- t(design%*%null+resid)
y <- t(design%*%effect+resid*varbatch)
batch.ind <- rep(c(1,2),each=8)
mod <- model.matrix(~as.factor(group.design))
mod0 <- model.matrix(~1,data=data.frame(design))
y1 <- ComBat(y,batch = batch.ind,mod = mod)
p1 <- f.pvalue(y1,mod,mod0)
p0 <- f.pvalue(y0,mod,mod0)
qqplot(p1,p0,xlab = "ComBat", ylab = "Benchmark",main = '(b) Moderate group effect \n strong mean and variance batch effect')
abline(a=0,b=1,lty=2)
par(mfrow=c(1,1))
## Use boxplots to show the comparison of groups
## boxplot(p1,p0)









