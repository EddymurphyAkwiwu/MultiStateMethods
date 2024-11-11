# rm(list = ls())
library(tidyverse)
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
library(igraph)
library(smms)


# We source startup.r which loads all required functions to generate the data using BayesTSM package
source('controllers/startup.r')



source('data_gen_func_Exp_smms.R')  # For exponentially distributed survival times


smms.dat = data
colnames(smms.dat)[1] <- c("patient") # Current version of smms package requires this

head(smms.dat)


head(cov.data)



#################################################################################

gg = igraph::graph_from_literal("1"--+"2"--+"3")

f_01 = function(param, x, t){dexp(t,exp(param[1]+param[2]*x[1]+param[3]*x[2]))}

f_12 = function(param, x, t){dexp(t,exp(param[4]+param[5]*x[1]+param[6]*x[2]))}


S_01 = function(param, x, t){(1-pexp(t,exp(param[1]+param[2]*x[1]+param[3]*x[2])))}
S_12 = function(param, x, t){(1-pexp(t,exp(param[4]+param[5]*x[1]+param[6]*x[2])))}




startval <- c(-2.0, -0.5, -0.5,  -1.2,-0.5,-0.5) 


mlo.last <- smms(startval,smms.dat,gg, cov.data, mc_cores = 1,hessian_matrix = F,abs_exact = FALSE)
ho.last <- hessian_matrix(mlo.last$opt$par,smms.dat, gg, cov.data ,mc_cores=1)


sdterr =  sqrt(diag(solve(ho.last)))

se = sdterr
se

est = mlo.last $opt$par  # ML estimates
est



load('Exp_TrueCDF_SCen_p_2.RData') # Assumed true CDF simulated in a large population setting

# load('Exp_TrueCDF_MCen_p_2.RData



F1 = ecdf(dat$X)
F2 = ecdf(dat$S)

t.x = quantile(dat$X, seq(0, 0.99, 0.01)) %>% as.numeric()

t.y = quantile(dat$S, seq(0, 0.99, 0.01))%>% as.numeric()

TrueX = F1(t.x)

TrueY = F2(t.y)


# source custom function for computing marginal CIF using pexp() function
source('smms_compute_Exp_marginal_CIF_func.R')


est.x = compute_exp_cdf(est[1], est[2], est[3],  t.x,cov.data )

est.y = compute_exp_cdf(est[4], est[5], est[6], t.y,cov.data )

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
