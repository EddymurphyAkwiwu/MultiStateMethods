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



source('data_gen_func_Exp_cthmm.R')  # For exponentially distributed survival times


head(dat2) # View data






head(cov.data)   # covariate(s)
# 
######## Intensity or rate or transition  distribution setup ############################
rates.setup=list()
rates.setup$num.states=5
design=read.csv("rate_design.csv",header=T,row.names=1)
design.matrix=(design[1:6,3:dim(design)[2]])
rates.setup$design.matrix=design.matrix
rates.setup$transition.codes=design[1:6,c(1,2)]
rates.setup$param.types=rep(0,times=10)
rates.setup$param.values = c(0.7646372, -1.6795936, 0.8461781, -2.6481639, -1.6541670, -3.8273887,
                             -0.3550389,-0.2508152, 0.250812, 0.8461781)+rnorm(n=10)/5


##############   #Is this part correct???? ############## 
rates.setup$variable.predictors = matrix(c(1,7,1,8,2,7,2,8,3,7,3,8,4,9,4,10,5,9,
                                           5,10,6,9,6,10),ncol=2,
                                         byrow=T,dimnames=list(seq(1:12),c("i","j")))

############## ############## ############## ############## ############## #############
rates.setup$covariate.array = get.covariate.array(design.matrix,cov.data,
                                                  rates.setup$variable.predictors)



rates.setup$deriv.array = get.deriv.array(rates.setup$covariate.array, rates.setup$param.types)

############################# initial distribution setup ########################################

init.setup=list()
init.setup$fixed.dist=c(1,0,0,0,0)  # Everyone starts from state 1

############################# #emission distribution setup ########################################

emission.setup=list()

# No misclassification
emission.setup$fixed.dist = matrix(c(1,0,0,1,0,0,0,1,0,0,1,0,0,0,1),byrow=T,nrow=5)




fit.4state=EM(rates.setup=rates.setup,
              init.setup=init.setup,
              emission.setup=emission.setup,
              the.data=dat2,
              num.subjects=length(dat2),
              num.states=5,
              num.obs.states=3,
              tol=1e-7,
              absorb.state=5,
              maxiter=1000)

out = fit.4state



covariance_4state=get_covariance(par=fit.4state$param,
                                 the.data=dat2,
                                 num.subjects=length(dat2),
                                 num.states=5,
                                 num.obs.states=3,
                                 rates.setup=rates.setup,
                                 emission.setup=emission.setup,
                                 init.setup=init.setup,
                                 DDO.setup=NULL,
                                 do.DDO=F)

#
# # The estimated standard error for the parameters (rates, emission, and initial distribution, in order):

se.est  = round(sqrt(diag(covariance_4state$covariance)),digits=3)

se.est






######################################################
# Get estimated rate matrices for each individual
######################################################

rates.list = get.rate.matrix.list(current.params=out$params[1:10],rate.setup=rates.setup,
                                  do.list=T)



######################################################


load('Exp_TrueCDF_SCen_p_2.RData') # Assumed true CDF simulated in a large population setting

# load('Exp_TrueCDF_MCen_p_2.RData



F1 = ecdf(dat$X)
F2 = ecdf(dat$S)

t.x = quantile(dat$X, seq(0, 0.99, 0.01))

t.y = quantile(dat$S, seq(0, 0.99, 0.01))

TrueX = F1(t.x)

TrueY = F2(t.y)





#############################   CDFs ######################

covar=covariance_4state$covariance[1:10,1:10]
num.transitions = dim(rates.setup$transition.codes)[1]
num.params = length(out$params)
transitions = rates.setup$transition.codes
num.states = 5



cdf.x = cdf.y = matrix(NA, nr = length(t.x) ,nc = n)


CDF_percentiles.x = CDF_percentiles.y = numeric()



for(i in 1: n){
  
  for(j in 1: length(t.x)){
    
    
    
    ##########################   State 1 to State 2 ####################
    
    param.deriv = rates.setup$deriv.array[,,i]
    
    rate.firstpassage.AB.1 = rates.list[[i]]
    rate.firstpassage.AB.1[3:5,] = 0
    
    cdfA.B.1 = sub_dist_times(start=t.x[j],end=t.x[j],states=3,alpha=c(1,0,0,0,0),rate.firstpassage.AB.1, length.out = 1)
    
    
    
    CDF_percentiles.x[j] = cdfA.B.1$dist 
    
   
    
    
    ##########################   State 2 to State 3 ####################
    param.deriv = rates.setup$deriv.array[,,i]
    
    rate.firstpassage.AB.2 = rates.list[[i]]
    rate.firstpassage.AB.2[c(5),] = 0  
    
    cdfA.B.2 = sub_dist_times(start=t.y[j],end=t.y[j],states=5,alpha=c(0,0,1,0,0),rate.firstpassage.AB.2,length.out = 1)
    
    
    CDF_percentiles.y[j] = cdfA.B.2$dist
    
  }
  
  
  
  cdf.x[,i] = CDF_percentiles.x
  
  cdf.y[,i] = CDF_percentiles.y
  

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
