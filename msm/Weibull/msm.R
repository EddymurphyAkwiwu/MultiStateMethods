# rm(list = ls())
library(tidyverse)
library(msm)
library(doParallel)
library(MCMCpack)
library(coda)
library(mvtnorm)
library(MHadaptive)
library(mcmcr)
library(MASS)
library(MCMCpack)
library(mvtnorm)
library(actuar)
library(parallel)
library(minqa)


# We source startup.r which loads all required functions to generate the data using BayesTSM package
source('controllers/startup.r')



source('data_gen_func_Weibull_msm.R')  # For Weibull distributed survival times



msm.dat = data
head(msm.dat)

# statetable.msm(state = state, subject = id, data = msm.dat)



#################################################################################
Q = rbind ( c( 0, 0.01, 0 ),
            c( 0, 0, 0.01 ),
            c( 0, 0, 0 ))




msm.mod =  msm( state ~ time, subject = id, data = msm.dat, qmatrix = Q,
                covariates = ~Z.1+Z.2, hessian = T,
                control = list(maxit=50000, fnscale=3000,reltol = 1e-16) )





load('Weibull_TrueCDF_SCen_p_2.RData') # Assumed true CDF simulated in a large population setting

# load('Weibull_TrueCDF_MCen_p_2.RData')







CDF_percentiles.x = CDF_percentiles.y = numeric()


F1 = ecdf(dat$X)
F2 = ecdf(dat$S)

t.x = quantile(dat$X, seq(0, 0.99, 0.01))

t.y = quantile(dat$S, seq(0, 0.99, 0.01))

TrueX = F1(t.x)

TrueY = F2(t.y)


Npersons = nrow(cov.data)

cdf.x = cdf.y = matrix(NA, nr = length(t.x) ,nc = Npersons)



for(j in 1: Npersons){
  
  
  for(i in 1: length(t.x)){
    prob.x       = ppass.msm(msm.mod, tot=t.x[i], covariates=list(Z.1= cov.data[j,'Z.1'],
                                                                  Z.2= cov.data[j,'Z.2']))
    
    prob.y       = ppass.msm(msm.mod, tot=t.y[i], covariates=list(Z.1= cov.data[j,'Z.1'],
                                                                  Z.2= cov.data[j,'Z.2']))
    
    
    CDF_percentiles.x[i] = prob.x[1,2] 
    CDF_percentiles.y[i] = prob.y[2,3]  
    
    
  }
  
  
  cdf.x[,j] = CDF_percentiles.x
  
  cdf.y[,j] = CDF_percentiles.y
  
}

est.x = rowMeans(cdf.x)
est.y = rowMeans(cdf.y)
#####################  plot cdfs ###################



plot(t.x,est.x,type="l",col=1,lty=1,lwd=2,ylim=c(0,1),
     xlab="time since entry in state 1",ylab="CDF",
     main = "Transition from state 1 to 2")

lines(t.x,TrueX,type="l",lwd=2,col='blue',lty=2)

legend("bottomright", legend=c("Model Estimate", "True"),
       col=c("black", "blue"), lty=c(1, 2), lwd=2)


plot(t.y,est.y,type="l",col=1,lty=1,lwd=2,ylim=c(0,1),
     xlab="time since entry in state 2",ylab="CDF",
     main = "Transition from state 2 to 3")

lines(t.y,TrueY,type="l",lwd=2,col='blue',lty=2)

legend("bottomright", legend=c("Model Estimate", "True"),
       col=c("black", "blue"), lty=c(1, 2), lwd=2)
