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
library(expm) 

# We source startup.r which loads all required functions to generate the data using BayesTSM package
source('controllers/startup.r')



source('data_gen_func_Exp_msm-phase.R')  # For exponentially distributed survival times




msm.dat = data
head(msm.dat)
length(unique(msm.dat$id))

# statetable.msm(state = state, subject = id, data = msm.dat)



#################################################################################
Q = rbind ( c( 0, 0.01, 0 ),
            c( 0, 0, 0.01 ),
            c( 0, 0, 0 ))



msm.mod0 =  msm( state ~ time, subject = id, data = msm.dat, qmatrix = Q,
                 covariates = ~Z.1+Z.2, hessian = FALSE )

Q = qmatrix.msm(msm.mod0)$estimate %>% as.matrix()

a = list(0.01,  c(  0.1, 0.01 ) )

b = list(0.01,  c(  0.01, 0.1 ) )

Qphase =  list(a,b)


msm.mod = msm( state ~ time, subject = id, data = msm.dat, qmatrix = Q,
               phase.inits =Qphase, hessian = FALSE,
               phase.states=c(1,2),
               covariates = ~Z.1+Z.2, constraint = list(Z.1=c(1,1,1,2,2,2), Z.2=c(1,1,1,2,2,2)),
               control = list(maxit=50000, fnscale=3000,reltol = 1e-16) )








load('Exp_TrueCDF_SCen_p_2.RData') # Assumed true CDF simulated in a large population setting

# load('Exp_TrueCDF_MCen_p_2.RData



CDF_percentiles.x = CDF_percentiles.y = numeric()


F1 = ecdf(dat$X)
F2 = ecdf(dat$S)

t.x = quantile(dat$X, seq(0, 0.99, 0.01)) %>% as.numeric()

t.y = quantile(dat$S, seq(0, 0.99, 0.01))  %>% as.numeric()

TrueX = F1(t.x)

TrueY = F2(t.y)


Npersons = nrow(cov.data)

cdf.x = cdf.y = matrix(NA, nr = length(t.x) ,nc = Npersons)



for(j in 1: Npersons){
  # tot=t.x[i]
  
  for(i in 1: length(t.x)){

    qmat   = qmatrix.msm(msm.mod,  covariates=list(Z.1= cov.data[j,'Z.1'],
                                                   Z.2= cov.data[j,'Z.2']))[[1]]
    
    
    qmat.x= qmat.y= matrix(nrow=5,ncol=5)   # define new transition intensity matrix 
    
    
    qmat.x[,] = qmat   # Extract the original estimated transition intensity matrix qmat
    
    qmat.x[3:4,] = 0  # set 3rd and 4th row of the original matrix intensity to zero 
    
    
    prob.x = expm(qmat.x*t.x[i])  # Compute the cdf for time t.x[i]
    
    
    qmat.y[,] = qmat   # Extract the original estimated transition intensity matrix qmat
    
    qmat.y[5,] = 0  # set 5th row of the original matrix intensity to zero 
    
    
    prob.y = expm(qmat.y*t.y[i])  # Compute the cdf for time t.y[i]
    
    
    CDF_percentiles.x[i] = prob.x[1,3]  ## State 1 to 2 (here 3 is the phase 1 of disease state 2)
    CDF_percentiles.y[i] = prob.y[3,5]  ## State 2 to 3 (here 5 is the final state (disease state 3))
    
    
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
