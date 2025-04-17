# Load required libraries
rm(list = ls())
library(parallel)
library(MASS)
library("mvtnorm")
library("mnormt")
library(ICSNP)
library("ICS")
library("SpatialNP")
library("MNM")
library(pcaPP)
library(Matrix)
library(matrixcalc)
library(matrixStats)
library(mixedCCA)
library(sparsepca)
###some CCA packages###
library(nscancor)
library(ccaPP)
library(CCA) #contain dataset 'nutrimouse'
library(latentcor)
library(robustHD)
library(corpcor) #make.positive.definite
library(psych)
library(corpcor)
#One replication#
Func_onesim <- function(i,X1,X2,split_ratio=4/5,BICtype=1){
  ##
  p1 <- dim(X1)[2];  
  p2 <- dim(X2)[2];
  n <- dim(X1)[1];
  d = p1 + p2;
  #
  X1 <- scale(X1, center = T,scale = T)
  X2 <- scale(X2, center = T,scale = T)
  X <- cbind(X1,X2)
  #train and test data
  set.seed(i)
  print(i)
  idx_n <- sample(1:n,ceiling(n*split_ratio));
  trainX <- X[idx_n,]
  testX <- X[-idx_n,]
  #
  trainX1 <- trainX[,1:p1]
  trainX2 <- trainX[,(p1+1):d]
  #
  testX1 <- testX[,1:p1]
  testX2 <- testX[,(p1+1):d]
  ##In mixedCCA
  kendallcca1 <- mixedCCA(trainX1, trainX2, 
                          type1 = "continuous", 
                          type2 = "continuous", 
                          BICtype = BICtype, 
                          nlamseq = 10)
  vx1<-kendallcca1$w1
  res_sv1_ken <- sum(vx1!=0)
  print(paste0('sv1_ken',res_sv1_ken))
  vy1<-kendallcca1$w2
  res_sv2_ken <- sum(vy1!=0)
  #
  if (vx1[1]<0) vx1<--vx1
  if (vy1[1]<0) vy1<--vy1
  #estimated kendall tau canonical vector
  kv <- c(vx1,vy1)
  ###########################
  ###sample covariance#######
  ###########################
  covx <- cov(trainX); 
  #make it semi-positive definite#
  covx <- as.matrix(Matrix::nearPD(covx)$mat)
  #make it symmetric#
  covx = (covx + t(covx))/2
  ###
  eSigmax<-covx[(1:p1),(1:p1)]
  eSigmay<-covx[(p1+1):d,(p1+1):d]
  eSigmaxy<-covx[(1:p1),(p1+1):d]
  ###
  lam<-lambdaseq_generate(nlamseq = 20,lam.eps = 0.01,eSigmax,eSigmay,eSigmaxy)
  lam1<-lam[[1]]
  lam2<-lam[[2]]
  ###
  n_train <- dim(trainX)[1]
  res <- find_w12bic(n_train,eSigmax,eSigmay,eSigmaxy,BICtype=BICtype,w1init=vx1,w2init=vy1,lamseq1=lam1,lamseq2=lam2)
  #
  evx<-res$w1
  evy<-res$w2
  #
  if (evx[1]<0) evx<--evx
  if (evy[1]<0) evy<--evy
  #
  evx <- evx/sqrt(sum(evx^2))
  res_sv1_cov <- sum(evx!=0)
  print(paste0('sv1_cov',res_sv1_cov))
  evy<-evy/sqrt(sum(evy^2))
  res_sv2_cov <- sum(evy!=0)
  #estimated canonical vector of sample covariance
  ev <- c(evx,evy)
  ########### Spatial sign covariance####
  scovx <- SCov(trainX)*d
  ##make it nearPD#
  scovx <- as.matrix(Matrix::nearPD(scovx)$mat)
  #make it symmetric
  scovx = (scovx + t(scovx))/2
  ##
  eSigmax <- scovx[(1:p1),(1:p1)]#*d
  eSigmay <- scovx[(p1+1):d,(p1+1):d]#*d#(d-p1)
  eSigmaxy <- scovx[(1:p1),(p1+1):d]#*d#p1
  #
  lam <- lambdaseq_generate(nlamseq = 20,lam.eps = 0.01,eSigmax,eSigmay,eSigmaxy)
  lam1<-lam[[1]]
  lam2<-lam[[2]]
  res<-find_w12bic(n_train,eSigmax,eSigmay,eSigmaxy,BICtype=BICtype,w1init=vx1,
                   w2init=vy1,lamseq1=lam1,lamseq2=lam2)
  svx<-res$w1
  svy<-res$w2
  #
  if (svx[1]<0) svx<--svx
  if (svy[1]<0) svy<--svy
  #
  esvx<-svx/sqrt(sum(svx^2)) #
  res_sv1_ss <- sum(esvx!=0)
  print(paste0('sv1_ss',res_sv1_ss))
  esvy<-svy/sqrt(sum(svy^2))
  res_sv2_ss <- sum(esvy!=0)
  sv <- c(esvx,esvy)
  ##compute canonical correlation of three methods##
  rhof<-function(u,v,S1,S2,S12)
  {t(u)%*%S12%*%v/sqrt(t(u)%*%S1%*%u)/sqrt(t(v)%*%S2%*%v)
  }
  #compute sample covariance of test data# 
  covx_test <- cov(testX);
  covx_test <- as.matrix(Matrix::nearPD(covx_test)$mat)
  covx_test = (covx_test + t(covx_test))/2
  #S1, S2 and S12
  eSigmax_cov <- covx_test[(1:p1),(1:p1)]
  eSigmay_cov <- covx_test[(p1+1):d,(p1+1):d]
  eSigmaxy_cov <- covx_test[(1:p1),(p1+1):d]
  #
  rho_cov_cov <- abs(rhof(evx,evy,eSigmax_cov,eSigmay_cov,eSigmaxy_cov))
  #compute spatial sign covariance of test data##
  scovx_test <- SCov(testX)*d
  scovx_test <- as.matrix(Matrix::nearPD(scovx_test)$mat)
  scovx_test = (scovx_test + t(scovx_test))/2
  ##S1, S2 and S12
  eSigmax_ss <- scovx_test[(1:p1),(1:p1)]
  eSigmay_ss<-scovx_test[(p1+1):d,(p1+1):d]
  eSigmaxy_ss<-scovx_test[(1:p1),(p1+1):d]
  rho_ss_ss <- abs(rhof(esvx,esvy,eSigmax_ss,eSigmay_ss,eSigmaxy_ss))
  #compute kendall tau correlation of test sample#
  kendallR_test <- estimateR(testX,type = 'continuous',use.nearPD = T);
  kendallR_test <- kendallR_test$R 
  ##S1,S2 and S12
  eSigmax_ken <- kendallR_test[(1:p1),(1:p1)]
  eSigmay_ken <- kendallR_test[(p1+1):d,(p1+1):d]
  eSigmaxy_ken <- kendallR_test[(1:p1),(p1+1):d]
  rho_ken_ken <- abs(rhof(vx1,vy1,eSigmax_ken,eSigmay_ken,eSigmaxy_ken))
  ##
  print(paste0('rho_ss_ss',round(rho_ss_ss,3),
               'rho_cov_cov',round(rho_cov_cov,3),
               'rho_ken_ken',round(rho_ken_ken,3)))
  ###
  res_cov <- c(rho_cov_cov,res_sv1_cov,res_sv2_cov);
  res_ken <- c(rho_ken_ken,res_sv1_ken,res_sv2_ken);
  res_ss <- c(rho_ss_ss,res_sv1_ss,res_sv2_ss);
  return(c(res_ss,res_cov,res_ken))
}
#we scale data, consider nonzero-variance cols#  
preprocess_for_scca <- function(X) {
  # Center and scale
  X <- scale(X) #default center and scale both TRUE
  
  # 3. Remove near-zero variance features
  nzv <- caret::nearZeroVar(X)
  if (length(nzv) > 0) X <- X[, -nzv]
  
  return(X)
}
###
library(CCA)
data(nutrimouse)
X1=as.matrix(nutrimouse$gene)
X2=as.matrix(nutrimouse$lipid)
# Basic preprocessing
X1 <- preprocess_for_scca(X1)
X2 <- preprocess_for_scca(X2)
library(parallel)
library(doParallel)
tic = proc.time()
nRep <- 500
cl <- makeCluster(16)
registerDoParallel(cl)
Res <- foreach(j=1:nRep,.combine = rbind,
               .packages = c('mixedCCA','SpatialNP'))%dopar%{
                 set.seed(j+1234)
                 res <- Func_onesim(j,X1,X2,split_ratio = 4/5,BICtype = 1) 
                 return(res)
               }
stopCluster(cl)
proc.time() - tic
Res <- as.data.frame(Res)
colnames(Res) = c('rho_ss','sv1_ss','sv2_ss',
                  'rho_cov','sv1_cov','sv2_cov',
                  'rho_ken','sv1_ken','sv2_ken')
apply(Res,2,mean)
apply(Res,2,sd)
