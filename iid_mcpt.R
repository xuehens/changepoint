## Multiple Changepoint Simulation
## Author: Xueheng Shi
## Date: 08/29/2019
## This program is designed for parallel simulation
## on Clemson SuperComputer


#rm(list=ls())

library(GA)
library(wbs, lib.loc = "/home/xuehens/software/r_lib/")
library(lars, lib.loc = "/home/xuehens/software/r_lib/")
library(glmnet,lib.loc = "/home/xuehens/software/r_lib/")

#library(wbs)
#library(lars)
#library(glmnet)
#setwd("C:\\Users\\Xueheng Shi\\Documents\\All MCPT")

##MDL.bin uses loc.ind as input which is binary with length N
# Note X[1] cannot be a CPT,
# so loc.ind[1] is initialized to be zero
# loc.ind = 0/1. 
# 1 corresponds a CPT location, 0 means not.

MDL.bin = function(loc.ind){
  loc.ind[1]=0
  N = length(Xt) #length of the series
  m = sum(loc.ind) #Number of CPTs
  
  if(m==0){
    ##Case 1, Zero Changepoint
    mu.hat = mean(Xt)
    sigma.hatsq = sum( (Xt-mu.hat)^2 )/N
    MDL.obj = 0.5*N*log(sigma.hatsq)  
  }
  else{
    tau.vec = loc.ind*(1:N) #convert binary to CPT location
    tau = tau.vec[tau.vec>0] #keep CPT locations only
    tau.ext = c(1,tau,(N+1)) #include CPT boundary 1 and N+1
    
    ## Split Xt to regimes/segments
    seg.len = diff(tau.ext) #length of each segments
    ff = rep(0:m, times=seg.len) ##create factors for segmentation
    Xseg = split(Xt, ff) ##Segmentation list
    mu.seg = unlist(lapply(Xseg,mean), use.names=F)
    mu.hat = rep(mu.seg, seg.len)
    sigma.hatsq = sum( (Xt-mu.hat)^2 )/N
    MDL.obj = 0.5*N*log(sigma.hatsq)+sum(log(diff(tau.ext))/2)+log(m)+sum(log(tau))
  }
  
  return(-MDL.obj)
}



##loc.ind is the binary input of length N
# Note X[1] cannot be a changepoint,
# so loc.ind[1] is always set to be zero
# loc.ind = 0/1. 1 corresponds a CPT location, 0 means not.
AIC.bin = function(loc.ind){
  loc.ind[1]=0
  N = length(Xt) #length of the series
  m = sum(loc.ind) #Number of CPTs
  
  if(m==0){
    ##Case 1, Zero Changepoint
    mu.hat = mean(Xt)
    sigma.hatsq = sum( (Xt-mu.hat)^2 )/N
    AIC.obj = 0.5*N*log(sigma.hatsq) + 4 #2*# of parameters
  }
  else{
    tau.vec = loc.ind*(1:N) #convert binary to CPT location
    tau = tau.vec[tau.vec>0] #keep CPT locations only
    tau.ext = c(1,tau,(N+1)) #include CPT boundary 1 and N+1
    
    ## Split Xt to regimes/segments 
    seg.len = diff(tau.ext) #length of each segments
    ff = rep(0:m, times=seg.len) ##create factors for segmentation
    Xseg = split(Xt, ff) ##Segmentation list
    mu.seg = unlist(lapply(Xseg,mean), use.names=F)
    mu.hat = rep(mu.seg, seg.len)
    sigma.hatsq = sum( (Xt-mu.hat)^2 )/N
    AIC.obj = 0.5*N*log(sigma.hatsq) + 4*m + 4
  }
  return(-AIC.obj)
}

##BIC
BIC.bin = function(loc.ind){
  loc.ind[1]=0
  N = length(Xt) #length of the series
  m = sum(loc.ind) #Number of CPTs
  
  if(m == 0){
    ##Case 1, Zero Changepoint
    mu.hat = mean(Xt)
    sigma.hatsq = sum( (Xt-mu.hat)^2 )/N
    BIC.obj = 0.5*N*log(sigma.hatsq)+ 2*log(N) 
  }
  else{
    tau.vec = loc.ind*(1:N) #convert binary to CPT location
    tau = tau.vec[tau.vec>0] #keep CPT locations only
    tau.ext = c(1,tau,(N+1)) #include CPT boundary 1 and N+1
    
    ## Split Xt to regimes/segments 
    seg.len = diff(tau.ext) #length of each segments
    ff = rep(0:m, times=seg.len) ##create factors for segmentation
    Xseg = split(Xt, ff) ##Segmentation list
    mu.seg = unlist(lapply(Xseg,mean), use.names=F)
    mu.hat = rep(mu.seg, seg.len)
    sigma.hatsq = sum( (Xt-mu.hat)^2 )/N
    BIC.obj = 0.5*N*log(sigma.hatsq) + (2*m + 2)*log(N)
  }
  return(-BIC.obj)
}

##loc.ind is the binary input of length N
# Note X[1] cannot be a changepoint,
# so loc.ind[1] is set to be zero
# loc.ind = 0/1. 1 corresponds a CPT location, 0 is not.
MBIC.bin = function(loc.ind){
  loc.ind[1]=0
  N = length(Xt) #length of the series
  m = sum(loc.ind) #Number of CPTs
  
  if(m==0){
    ##Case 1, Zero Changepoint
    mu.hat = mean(Xt)
    sigma.hatsq = sum( (Xt-mu.hat)^2 )/N
    MBIC.obj = 0.5*N*log(sigma.hatsq)  
  }
  else{
    tau.vec = loc.ind*(1:N) #convert binary to CPT location
    tau = tau.vec[tau.vec>0] #keep CPT locations only
    tau.ext = c(1,tau,(N+1)) #include CPT boundary 1 and N+1
    
    ## Split Xt to regimes/segments 
    seg.len = diff(tau.ext) #length of each segments
    ff = rep(0:m, times=seg.len) ##create factors for segmentation
    Xseg = split(Xt, ff) ##Segmentation list
    mu.seg = unlist(lapply(Xseg,mean), use.names=F)
    mu.hat = rep(mu.seg, seg.len)
    sigma.hatsq = sum( (Xt-mu.hat)^2 )/N
    MBIC.obj = 0.5*N*log(sigma.hatsq) + 1.5*m*log(N) + 0.5*sum(log(seg.len/N))
  }
  return(-MBIC.obj)
}


##LASSO changepoint 
MCPT.LASSO = function(Xt){ 
  ##Create Design matrix
  Nd = length(Xt)
  
  ##differencing to compute phi.hat
  ##differencing is not needed for IID White Noise
  #DXt = diff(Xt)
  #length(DXt)
  
  ##estimating phi
  #acv.dh = acf(DXt, lag.max = 1, type = "covariance", plot = F)$acf
  #phi.dhat = 2*acv.dh[2]/acv.dh[1] + 1
  #Yt = as.vector( Xt - phi.dhat*c(0, Xt[-1]) )
  #plot(Yt)
  
  Yt = Xt ##change name to comply with linear model setting
  Xd = matrix(0, nrow = Nd, ncol=Nd-1)
  Xd[lower.tri(Xd,diag = T)] = 1
  
  ##Use Least Angle Regression to selection best lambda
  lasso.path  = lars(Xd,Yt,intercept=F)
  n = Nd
  betas    = lasso.path$beta
  df       = lasso.path$df
  MSE      = lasso.path$RSS/n
  bic      = log(n)*df + n*log(MSE)
  bestb    = which.min(bic[seq(1, 0.5*Nd )])
  
  beta_lasso = betas[bestb, ]
  
  ##perform adaptive lasso
  if ( sum(abs(beta_lasso)) == 0 ){
    beta_adalasso = beta_lasso
  } else {
    ##If beta_lasso=0, no weight, so convert 0 to 1
    beta_lasso[beta_lasso==0] = 10^(-4) 
    ##Create weight  
    weight.h = abs(1/beta_lasso)
    #adalasso=ic.glmnet(Xd,Yt,crit="bic",penalty.factor=what)
    adalasso = glmnet(x = Xd, y = Yt,
                      nlambda = Nd,
                      intercept = F,
                      ## 'alpha = 1' is lasso penalty, '.. = 0' ridge penalty.
                      alpha = 1, 
                      penalty.factor = weight.h)
    
    ##how the following result comes from???
	##I will insert a reference later
    tLL = adalasso$nulldev - deviance(adalasso)
    k = adalasso$df
    n = adalasso$nobs
    BIC = log(n)*k - tLL
    beta_adalasso = adalasso$beta[ ,which.min(BIC[seq(1:0.5*Nd)])]
  }
  #ncpt = sum( abs(beta_adalasso) != 0 )
  return(beta_adalasso)
}  



#start.time = proc.time()
AIC.comb = NULL
BIC.comb = NULL
MDL.comb = NULL
MBIC.comb = NULL

WBS.comb = NULL
BS.comb = NULL

LASSO.comb = NULL

#Store Xt data for further analysis
#Xt.bind = NULL

kk = 1 ##No. of iterations on each CPU
for(i in 1:kk){
  Xt = rnorm(100)
  #Xt = rnorm(100)+c(rep(0,25), rep(1,25), rep(0,25), rep(1,25))
  #Xt = rnorm(100)+c(rep(0,25), rep(2,25), rep(0,25), rep(2,25))
  #Xt = rnorm(100)+c(rep(0,50), rep(1,50))
  #Xt = rnorm(100)+c(rep(0,50), rep(2,50))

  N = length(Xt)
  
  #Run AICGA
  AICGA = ga(type="binary", fitness = AIC.bin,
             nBits = N, maxiter = 30000, run = 3000,
             popSize = 200, monitor = T)
  AICsol = (AICGA@solution)[1, ] ##can be 2 sols due to X[1], take 1st sol
  AICsol[1] = 0
  
  mAIC = sum(AICsol)
  AIC.loc = AICsol*(1:N)
  tAIC = ifelse(mAIC==0, "NA", toString(AIC.loc[AIC.loc>0]) )
  AIC.comb = rbind(AIC.comb, cbind("AIC", mAIC, tAIC))
  
  
  #Run BICGA
  BICGA = ga(type="binary", fitness = BIC.bin,
             nBits = N, maxiter = 30000, run = 3000,
             popSize = 200, monitor = T)
  BICsol = (BICGA@solution)[1, ] ##can be 2 sols since X[1] is free, take 1st sol
  BICsol[1] = 0 
  
  mBIC = sum(BICsol)
  BIC.loc = BICsol*(1:N)
  tBIC = ifelse(mBIC==0, "NA", toString(BIC.loc[BIC.loc>0]) )
  BIC.comb = rbind(BIC.comb, cbind("BIC", mBIC, tBIC))
  
  
  #Run MBICGA
  MBICGA = ga(type="binary", fitness = MBIC.bin,
              nBits = N, maxiter = 30000, run = 3000,
              popSize = 200, monitor = T)
  MBICsol = (MBICGA@solution)[1, ] ##can be 2 sols since X[1] is free, take 1st sol
  MBICsol[1] = 0 
  
  mMBIC = sum(MBICsol)
  MBIC.loc =MBICsol*(1:N)
  tMBIC = ifelse(mMBIC==0, "NA", toString(MBIC.loc[MBIC.loc>0]) )
  MBIC.comb = rbind(MBIC.comb, cbind("MBIC", mMBIC, tMBIC))
  
  
  #MDL
  MDLGA = ga(type="binary", fitness = MDL.bin,
             nBits = N, maxiter = 30000, run=3000,
             popSize = 200, monitor = T)
  MDLsol = (MDLGA@solution)[1,] ##can be 2 sols
  MDLsol[1]=0
  
  mMDL = sum(MDLsol)
  MDL.loc=MDLsol*(1:N)
  tMDL = ifelse(sum(MDLsol)==0, "NA", toString(MDL.loc[MDL.loc>0]) )
  MDL.comb = rbind(MDL.comb, cbind("MDL", mMDL, tMDL))
  
  
  ##wbs
  ww = wbs(Xt, M=5000)
  wbscpt = changepoints(ww)
  ##result of wbs
  mWBS = wbscpt$no.cpt.th       #number of CPT
  tWBS = toString(as.character(sort(unlist(wbscpt$cpt.th))))
  WBS.comb = rbind(WBS.comb, cbind("WBS", mWBS, tWBS))
  
  ##bs
  bs = sbs(Xt)
  bscpt = changepoints(bs)
  mBS = bscpt$no.cpt.th #number of CPT
  tBS = toString(as.character(sort(unlist(bscpt$cpt.th))))
  BS.comb = rbind(BS.comb, cbind("BS", mBS, tBS))
  
  
  ##LASSO MCPT method
  lasso = MCPT.LASSO(Xt)
  mlasso = sum(abs(lasso) !=0 )
  tlasso = ifelse(mlasso==0, "NA", 
                  toString(as.character(which(lasso !=0))))
  LASSO.comb = rbind(LASSO.comb, cbind("LASSO", mlasso, tlasso))
  
}


colnames(AIC.comb) = c("AIC", "#CPT", "tau")
colnames(BIC.comb) = c("BIC", "#CPT", "tau")
colnames(MBIC.comb) = c("MBIC", "#CPT", "tau")
colnames(MDL.comb) = c("MDL", "#CPT", "tau")
colnames(WBS.comb) = c("WBS", "#CPT", "tau")
colnames(BS.comb)  = c("BS",  "#CPT", "tau")
colnames(LASSO.comb) = c("LASSO", "#CPT", "tau")


MCPT.comb = cbind(AIC.comb,
                  BIC.comb,
                  MBIC.comb,
                  MDL.comb,
                  WBS.comb,
                  BS.comb,
                  LASSO.comb)

write.csv(MCPT.comb, "MCPT=0_Run=100_N=100_IIDWN_sigma=1_delta=0.csv")










