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



source('data_gen_func_Weibull_BayesTSM.R')  # For Weibull distributed survival times

head(dat)





# Function bts_survreg is used to estimate the model (parallel processing requires at least 3 free CPUs)
mod.ini          = bts_survreg(d              = dat$d, # censoring indicator, assumes values 1, 2, 3
                           L              = dat$L, # Time of left censoring; assumes time of last visit of d=1
                           R              = dat$R, # Time of right censoring; assumes inf if d=1
                           Z.X            = dat[,c('Z.1','Z.2')], # Covariates of X
                           Z.S            = dat[,c('Z.1','Z.2')], # Covariates of S
                           mc             = 1e3,       # MCMC draws (can be updated later, see below), half will be dropped for burn-in (other values can also be specified using brunin argument)
                           chains         = 3,         # Number of parallel MCMC chains
                           thin           = 5,         # The function returns thinned and unthinned chains; there the thinning interval may be given
                           do.seperate.MH = F,         # Whether the metropolis step should be done jointly (F) or seperately for parameters of X and S
                           prop.sd.X      = 0.01,      # The proposal standard deviation of the normal distribution used for the metropolis step (if do.seperate.MH = T also prop.sd.S should be tuned)
                           beta.prior.X   = 4,         # The degrees of freedom of a t-distribution for prior of model betas and intercept (X)
                           beta.prior.S   = 4,         # The degrees of freedom of a t-distribution for prior of model betas and intercept (S)
                           sig.prior.X    = sqrt(10),  # The sd=tau of a half normal distribution N+(0,tau^2) (X)
                           sig.prior.S    = sqrt(10),  # The sd=tau of a half normal distribution N+(0,tau^2) (S)
                           dist.X         = 'weibull', # Distribution of X
                           dist.S         = 'weibull', # Distribution of S
                           fix.sigma.X    = F,         # Should sigma of X be fixed at its prior value?
                           fix.sigma.S    = F,         # Should sigma of S be fixed at its prior value? (use e.g. to specify exponential distribution)
                           parallel       = T,         # TRUE if chains should be run on seperate CPUs
                           beta.prior     = 't')       # Alternatively a normal prior for beta can be specified with 'norm'



# Second, we run search.prop.sd implementing the heuristic search
s   = search.prop.sd(mod.ini, acc.bounds.X =c(0.22,0.24)) # acc.bounds.X gives acceptable acceptance rate bounds as described in the supplemental material
s$prop.sd.X # The tuned acceptance probability


mc.ini = 1e5

mod          = bts_survreg(d              = dat$d, # censoring indicator, assumes values 1, 2, 3
                         L              = dat$L, # Time of left censoring; assumes time of last visit of d=1
                         R              = dat$R, # Time of right censoring; assumes inf if d=1
                         Z.X            = dat[,c('Z.1','Z.2')], # Covariates of X
                         Z.S            = dat[,c('Z.1','Z.2')], # Covariates of S
                         mc             = mc.ini ,       # MCMC draws (can be updated later, see below), half will be dropped for burn-in (other values can also be specified using brunin argument)
                         chains         = 3,         # Number of parallel MCMC chains
                         thin           = 5,         # The function returns thinned and unthinned chains; there the thinning interval may be given
                         do.seperate.MH = F,         # Whether the metropolis step should be done jointly (F) or seperately for parameters of X and S
                         prop.sd.X      = s$prop.sd.X,# The proposal standard deviation of the normal distribution used for the metropolis step (if do.seperate.MH = T also prop.sd.S should be tuned)
                         beta.prior.X   = 4,         # The degrees of freedom of a t-distribution for prior of model betas and intercept (X)
                         beta.prior.S   = 4,         # The degrees of freedom of a t-distribution for prior of model betas and intercept (S)
                         sig.prior.X    = sqrt(10),  # The sd=tau of a half normal distribution N+(0,tau^2) (X)
                         sig.prior.S    = sqrt(10),  # The sd=tau of a half normal distribution N+(0,tau^2) (S)
                         dist.X         = 'weibull', # Distribution of X
                         dist.S         = 'weibull', # Distribution of S
                         fix.sigma.X    = F,         # Should sigma of X be fixed at its prior value?
                         fix.sigma.S    = F,         # Should sigma of S be fixed at its prior value? (use e.g. to specify exponential distribution)
                         parallel       = T,         # TRUE if chains should be run on seperate CPUs
                         beta.prior     = 't')       # Alternatively a normal prior for beta can be specified with 'norm'


# Updating MCMC chain 

# mod.updated = bts_survreg(prev.run = mod, mc = 1e5) # An update for another 10^5 mc draws
# mod = mod.updated




load('Weibull_TrueCDF_SCen_p_2.RData') # Assumed true CDF simulated in a large population setting

# load('Weibull_TrueCDF_MCen_p_2.RData')


F1 = ecdf(dat$X)
F2 = ecdf(dat$S)

t.x = quantile(dat$X, seq(0, 0.99, 0.01)) %>% as.numeric()

t.y = quantile(dat$S, seq(0, 0.99, 0.01))%>% as.numeric()

TrueX = F1(t.x)

TrueY = F2(t.y)



# Marginal CIF estimates

cdf.x = get.pCIF.p(mod=mod, pst.samples = 1e4, q= t.x)
cdf.y = get.pCIF.p(mod = mod, pst.samples = 1e4, q = t.y)


# Extract median estimates

est.x = cdf.x$med.cdf.x 

est.y = cdf.y$med.cdf.s 

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
