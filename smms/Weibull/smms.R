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



source('data_gen_func_Weibull_smms.R')  # For exponentially distributed survival times


smms.dat = data
colnames(smms.dat)[1] <- c("patient") # Current version of smms package requires this

head(smms.dat)


head(cov.data)



#################################################################################

gg = igraph::graph_from_literal("1"--+"2"--+"3")


f_01 = function(param, x, t){dweibull(t,exp(param[1]),
                                      exp(param[2]+param[3]*x[1]+param[4]*x[2]  ))}

f_12 = function(param, x, t){dweibull(t,exp(param[5]),
                                      exp(param[6]+param[7]*x[1]+param[8]*x[2]))}



S_01 = function(param, x, t){(1-pweibull(t,exp(param[1]),
                                         exp(param[2]+param[3]*x[1]+param[4]*x[2])))}
S_12 = function(param, x, t){(1-pweibull(t,exp(param[5]),
                                         exp(param[6]+param[7]*x[1]+param[8]*x[2])))}






startval <- c(1.6, 3.0, 0.5, 0.5,  1.2,1.2,0.5,0.5) 


mlo.last <- smms(startval,smms.dat,gg, cov.data, mc_cores = 1,hessian_matrix = F,abs_exact = FALSE)
ho.last <- hessian_matrix(mlo.last$opt$par,smms.dat, gg, cov.data ,mc_cores=1)


sdterr =  sqrt(diag(solve(ho.last)))

# se = sdterr
# se

# Transform the # ML estimates of the Weibull shapes parameter to the original scale (i.e., not the log scale)

est = c(exp(mlo.last $opt$par[1]), mlo.last $opt$par[2:4], exp(mlo.last $opt$par[5]), mlo.last $opt$par[6:8])
est



load('Weibull_TrueCDF_SCen_p_2.RData') # Assumed true CDF simulated in a large population setting

# load('Weibull_TrueCDF_MCen_p_2.RData



F1 = ecdf(dat$X)
F2 = ecdf(dat$S)

t.x = quantile(dat$X, seq(0, 0.99, 0.01)) %>% as.numeric()

t.y = quantile(dat$S, seq(0, 0.99, 0.01))%>% as.numeric()

TrueX = F1(t.x)

TrueY = F2(t.y)


# source custom function for computing marginal CIF using pexp() function
source('smms_compute_Weibull_marginal_CIF_func.R')


est.x = compute_Weibull_cdf(est[1], est[2], est[3], est[4], t.x,cov.data )

est.y = compute_Weibull_cdf(est[5], est[6], est[7], est[8], t.y,cov.data )

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
