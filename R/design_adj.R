#' A function for correcting p-values for unbalanced design & ComBat in sva
#' 
#' This function adjusts for unbalanced design matrix to calculate correct f-statistics/p-values
#' Each row of a data matrix comparing the nested models
#' defined by the design matrices for the alternative (mod) and and null (mod0) cases.
#' The columns of mod0 must be a subset of the columns of mod.
#' The user needs to provide group and batch indicators in order to adjust for unbalanced group-batch design  
#' 
#' @param dat The transformed data matrix with the variables in rows and samples in columns
#' @param mod The model matrix being used to fit the data
#' @param mod0 The null model being compared when fitting the data
#' @param treat The treatment indicator, i.e., the group design 
#' @param batch The batch indicator, i.e., the batch design 
#' @param ref.batch The reference batch, the default is null (no reference batch)
#' @param noise The noise input, if not set will use the default value 1/n
#' @param output The desired output, the default is adjusted p-values (could also be sample covariance matrix)
#'
#' @return result could be a vector of F-statistic p-values one for each row of dat or sample covariance matrix
#' 
#' @examples 
#' 
#' library(bladderbatch)
#' data(bladderdata)
#' dat <- bladderEset[1:50,]
#' 
#' pheno = pData(dat)
#' edata = exprs(dat)
#' mod = model.matrix(~as.factor(cancer), data=pheno)
#' mod0 = model.matrix(~1,data=pheno)
#' 
#' 
#' pValues = f.pvalue(edata,mod,mod0)
#' qValues = p.adjust(pValues,method="BH")
#'
#' batch = pheno$batch
#' treat = ifelse(pheno$cancer=='Cancer',1,0)
#' d = ComBat(dat=edata, batch=batch, mod=mod)
#' mod = model.matrix(~as.factor(treat), data=pheno)
#' mod0 =  model.matrix(~1, data=pheno)
#' adjusted_pValues = design_adj(d,mod,mod0,treat = treat,batch = batch)
#' adjusted_pValues1 = design_adj(d,mod,mod0,treat = treat,batch = batch, ref = 1)
#' cov_matrix = design_adj(d,mod,mod0,treat = treat,batch = batch,output = "matrix")
#'
#' @export
#' 

design_adj <- function(dat,mod,mod0,treat,batch,ref.batch=NULL,noise=NULL,output="pvalue"){
  if (output=="pvalue"|output=="matrix"){
  n <- dim(dat)[2]
  m <- dim(dat)[1]
  df1 <- dim(mod)[2]
  df0 <- dim(mod0)[2]
  p <- rep(0,m)
  Id <- diag(n)
  noise <- ifelse(is.null(noise),1/n,noise)
  nbatch <- max(batch)
  if (is.null(ref.batch)){
    batchm <- mat.or.vec(n,nbatch)
    for (i in 1:nbatch){
        batchm[,i] <- ifelse(batch==i,1,0)
    }
    h1 <- treat%*%solve(t(treat)%*%treat)%*%t(treat)
    h12 <- batchm%*%solve(t(batchm)%*%(Id-h1)%*%batchm)%*%t(batchm)%*%(Id-h1)
  } else {
    batchm <- mat.or.vec(n,nbatch-1)
    batchv <- unique(batch)
    batchv <- batchv[batchv!=ref.batch]
    treat1 <- as.matrix(cbind(1,treat))
    for (i in 1:(nbatch-1)){
      batchm[,i] <- ifelse(batch==batchv[i],1,0)
    }
    h1 <- treat1%*%solve(t(treat1)%*%treat1)%*%t(treat1)
    h12 <- batchm%*%solve(t(batchm)%*%(Id-h1)%*%batchm)%*%t(batchm)%*%(Id-h1)
  }
  reduce <- Id-h12
  mt <- reduce%*%t(reduce)
  v <- eigen(mt)$vectors
  d <- eigen(mt)$'values'
  err <- sum(d)*noise
  d <- ifelse(d<0.0001,err,d)
  k <- v%*%diag(d)%*%t(v)
  if (output=="matrix"){
    return(k)
  } else {
  s <- solve(t(chol(k)))
  dat <-  t(s%*%t(dat))
  mod <- s%*%mod
  mod0 <- s%*%mod0
  resid <- dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*% t(mod))
  rss1 <- rowSums(resid*resid)
  rm(resid)
  resid0 <- dat %*% (Id - mod0 %*% solve(t(mod0) %*% mod0) %*% t(mod0))
  rss0 <- rowSums(resid0*resid0)
  rm(resid0)
  fstats <- ((rss0 - rss1)/(df1-df0))/(rss1/(n-df1))
  p <-  1-pf(fstats,df1=(df1-df0),df2=(n-df1))
  return(p)
  }
  } else {
    stop("Need to specify the output type; Must be either pvalue or matrix")
  }
}



