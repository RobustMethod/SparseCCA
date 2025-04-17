# Load required libraries

library(parallel)
library(MASS)
library("mvtnorm")
library("mnormt")
library(ICSNP)
library("ICS")
library("SpatialNP")
library("MNM")
library(Matrix)

library(mixedCCA)

tt1<-proc.time()[1]

# A small toy example

Res1<-matrix(0,6,9)
Res2<-matrix(0,6,18)
Res3<-matrix(0,18,12)

SIM<-100

for (di in 1:2)
{
d<-400*di
p<-d/2
# Generate block diagonal covariance matrix with 5 blocks

B<-5
pb<-p/B

Sigmab<- matrix(0,pb,pb)
for(i in 1:pb){
Sigmab[i,] <- 1:(pb)-i
}
Sigmab <- 0.8^abs(Sigmab)

Sigmax<-as.matrix(bdiag(Sigmab,Sigmab,Sigmab,Sigmab,Sigmab))
Sigmay<-Sigmax

lambda1<-0.9
vx<-c(1,0,0,0,0,1,rep(0,4),1,rep(0,p-11))/sqrt(3)
vy<-vx

vx<-vx/c(sqrt(t(vx)%*%Sigmax%*%vx))

vy<-vy/c(sqrt(t(vy)%*%Sigmay%*%vy))

v<-c(vx,vy)

kd<-50

Lambda<-diag(kd)*0.05
Qx<-eigen(Sigmax)$vector
Lx<-diag(1/sqrt(eigen(Sigmax)$values))
Vx<-(Qx%*%Lx)[,1:kd]

Qy<-eigen(Sigmay)$vector
Ly<-diag(1/sqrt(eigen(Sigmay)$values))
Vy<-(Qy%*%Ly)[,1:kd]


Sigmaxy<-Sigmax%*%vx%*%t(vy)%*%Sigmay*lambda1+Sigmax%*%Vx%*%Lambda%*%t(Vy)%*%Sigmay

Sigma<-matrix(0,d,d)
Sigma[1:p,1:p]<-Sigmax
Sigma[(p+1):d,(p+1):d]<-Sigmay
Sigma[1:p,(p+1):d]<-Sigmaxy
Sigma[(p+1):d,1:p]<-t(Sigmaxy)

eigen(Sigma)$values

for (ni in 1:3)
{

n<-100*ni


for (mi in 1:3)
{


Resultrho<-matrix(0,SIM,3)
ResultL<-matrix(0,SIM,6)
Result<-matrix(0,SIM,12)
for (kk in 1:SIM)
{

X <- rmvnorm(n, rep(0, d), Sigma)
 if (mi == 2) 
{
    X <- rmt(n, rep(0, d), Sigma, df = 3)/sqrt(3)
} 
if (mi == 3) 
{
    xr <- rbinom(n, 1, 0.8)
    X <- (X * xr + X * 10 * (1 - xr)) / sqrt(20.8)

}

covx<-cov(X)

A<-covx
A[1:p,1:p]<-0
A[(p+1):d,(p+1):d]<-0

B<-covx
B[1:p,(p+1):d]<-0
B[(p+1):d,1:p]<-0

k<-6

lambda<-sqrt(log(d)/n)
K<-1
#P<-initial.convex(A, B, lambda, K, nu = 1, epsilon = 0.005, maxiter = 1000, trace = FALSE)$Pi

#init<-eigen(P)$vector[,1]

#if (abs(t(init)%*%A%*%init)>0)
#{
#ev<-rifle(A, B, init, k, eta = 0.01, convergence = 0.001, maxiter = 5000)
#}


##### Kendall tau matrix


X1<-X[,1:p]
X2<-X[,(1+p):d]
kendallcca1 <- mixedCCA(X1, X2, type1 = "continuous", type2 = "continuous", BICtype = 1, nlamseq = 10)

vx1<-kendallcca1$w1
vy1<-kendallcca1$w2

if (vx1[1]<0) vx1<--vx1
if (vy1[1]<0) vy1<--vy1

kv<-c(vx1,vy1)

#ev<-rifle(A, B, kv, k, eta = 0.01, convergence = 0.001, maxiter = 5000)

####### covariance matrix


covx<-cov(X)

eSigmax<-covx[(1:p),(1:p)]
eSigmay<-covx[(p+1):d,(p+1):d]
eSigmaxy<-covx[(1:p),(p+1):d]


lam<-lambdaseq_generate(nlamseq = 20,lam.eps = 0.01,eSigmax,eSigmay,eSigmaxy)
lam1<-lam[[1]]
lam2<-lam[[2]]

res<-find_w12bic(n,eSigmax,eSigmay,eSigmaxy,BICtype=1,w1init=vx1,w2init=vy1,lamseq1=lam1,lamseq2=lam2)

evx<-res$w1
evy<-res$w2

if (evx[1]<0) evx<--evx
if (evy[1]<0) evy<--evy

evx<-evx/sqrt(sum(evx^2))

evy<-evy/sqrt(sum(evy^2))

ev<-c(evx,evy)


########### Spatial sign method
scovx<-SCov(X)

eSigmax<-scovx[(1:p),(1:p)]*p
eSigmay<-scovx[(p+1):d,(p+1):d]*p
eSigmaxy<-scovx[(1:p),(p+1):d]*p


lam<-lambdaseq_generate(nlamseq = 20,lam.eps = 0.01,eSigmax,eSigmay,eSigmaxy)
lam1<-lam[[1]]
lam2<-lam[[2]]

res<-find_w12bic(n,eSigmax,eSigmay,eSigmaxy,BICtype=1,w1init=vx1,w2init=vy1,lamseq1=lam1,lamseq2=lam2)

svx<-res$w1
svy<-res$w2

if (svx[1]<0) svx<--svx
if (svy[1]<0) svy<--svy

esvx<-svx/sqrt(sum(svx^2))

esvy<-svy/sqrt(sum(svy^2))

#c(sqrt(sum((evx-vx)^2)),sqrt(sum((kv[1:p]-vx)^2)),sqrt(sum((esvx-vx)^2)))

#c(sqrt(sum((evy-vy)^2)),sqrt(sum((kv[(1+p):d]-vx)^2)),sqrt(sum((esvy-vy)^2)))

cosd<-function(a,b)
{
 1-sum(a*b)/sqrt(sum(a^2))/sqrt(sum(b^2))
}

c(cosd(evx,vx),cosd(vx1,vx),cosd(esvx,vx))*100
c(cosd(evy,vy),cosd(vy1,vy),cosd(esvy,vy))*100

rhof<-function(u,v,S1,S2,S12)
{
t(u)%*%S12%*%v/sqrt(t(u)%*%S1%*%u)/sqrt(t(v)%*%S2%*%v)
}

Resultrho[kk,1]<-lambda1-abs(rhof(evx,evy,Sigmax,Sigmay,Sigmaxy))
Resultrho[kk,2]<-lambda1-abs(rhof(vx1,vy1,Sigmax,Sigmay,Sigmaxy))
Resultrho[kk,3]<-lambda1-abs(rhof(esvx,esvy,Sigmax,Sigmay,Sigmaxy))


lf<-function(u,eu,S)
{
1-abs(t(eu)%*%S%*%u)/sqrt(t(u)%*%S%*%u)
}

ResultL[kk,1]<-lf(evx,vx,Sigmax)
ResultL[kk,2]<-lf(evy,vy,Sigmay)
ResultL[kk,3]<-lf(vx1,vx,Sigmax)
ResultL[kk,4]<-lf(vy1,vy,Sigmay)
ResultL[kk,5]<-lf(esvx,vx,Sigmax)
ResultL[kk,6]<-lf(esvy,vy,Sigmay)

tpr<-function(u,v)
{
 length(which(u*v!=0))/length(which(v!=0))
}

tnr<-function(u,v)
{
 length(intersect(which(u==0),which(v==0)))/length(which(v==0))
}

Result[kk,1]<-tpr(evx,vx)
Result[kk,2]<-tnr(evx,vx)

Result[kk,3]<-tpr(evy,vy)
Result[kk,4]<-tnr(evy,vy)

Result[kk,5]<-tpr(vx1,vx)
Result[kk,6]<-tnr(vx1,vx)

Result[kk,7]<-tpr(vy1,vy)
Result[kk,8]<-tnr(vy1,vy)

Result[kk,9]<-tpr(esvx,vx)
Result[kk,10]<-tnr(esvx,vx)

Result[kk,11]<-tpr(esvy,vy)
Result[kk,12]<-tnr(esvy,vy)

}

Resultrho[which(is.na(Resultrho)==TRUE)]<-1
ResultL[which(is.na(ResultL)==TRUE)]<-1
Result[which(is.na(Result)==TRUE)]<-1

Res1[ni+3*(di-1),(1:3)+3*(mi-1)]<-colMeans(abs(Resultrho*100))

Res2[ni+3*(di-1),(1:6)+6*(mi-1)]<-colMeans(ResultL)*100

Res3[ni+3*(di-1)+6*(mi-1),]<-colMeans((1-Result)*100)

}
}
}

Res1

Res2

Res3

#write.table(round(Res1,1),"ccarho2.txt",sep="&")
#write.table(round(Res2,1),"ccal2.txt",sep="&")
#write.table(round(Res3,1),"ccasle2.txt",sep="&")

tt2<-proc.time()[1]

tt2-tt1


