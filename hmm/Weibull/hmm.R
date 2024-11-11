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

library(eddy, include.only = 'fn_log_post') # or install.packages(eddy_0.0.0.9000.tar.gz, repos = NULL, type="source")


# We source startup.r which loads all required functions to generate the data using BayesTSM package
source('controllers/startup.r')

source("mcmc_routine_p2.r") 

source('data_gen_func_Weibull_hmm.R')  # For Weibull distributed survival times


head(data)



load('Weibull_TrueCDF_SCen_p_2.RData')

# load('Weibull_TrueCDF_MCen_p_2.RData')

F1 = ecdf(dat$X)
F2 = ecdf(dat$S)

t.x.p = quantile(dat$X, seq(0, 0.99, 0.01))%>% as.numeric()
t.x.end = length(t.x.p)
t.x = seq(0,ceiling(t.x.p[t.x.end ]), length.out=100)
lx = length(t.x)


t.y.p = quantile(dat$S, seq(0, 0.99, 0.01))%>% as.numeric()
t.y.end = length(t.y.p)
t.y = seq(0,ceiling(t.y.p[t.y.end ]), length.out=100)
ly = length(t.y)


TrueX = F1(t.x)

TrueY = F2(t.y)
#############################################################

# head(data)



n.chains = 3





# n.post.samples = 1e1
# steps = 100
# burnin = 50
# 
n.post.samples = 1e2
steps = 1e5
burnin = 1e4


# *************************************************************************
disc = TRUE  #--> piece-wise time homogeneous solution

# 
init_par = c(-3, 0, 0, 0, # beta(0, x), beta(1, x), beta(2, x), beta(2, x_time)
             -3, 0, 0, 0) # beta(0, t), beta(1, t), beta(2, t), beta(2, t_time)

par_index = list( beta=1:8 )

# Defining the mean and variance for the uninformed Gaussian priors for the MH update
prior_par = data.frame( prior_mean=rep( 0, length(init_par)),
                        prior_sd=rep( 20, length(init_par)))


# Transform the in put data into hmm format


################### transform to hmm format #######################################
censor_times_btsm <- function(tat) {
  min_t = 0
  max_t = floor(max(tat))
  new_time = seq(min_t, max_t, by = 1)
  
  return(new_time)
}




disc_time <- sapply(data$time, floor)

obstrue <- rep(0,nrow(data))

hold <- cbind(data,obstrue,disc_time)

tempRow <- rep(0,ncol(hold))
names(tempRow) <- colnames(hold)

num <- 1
data_added_rows <- NULL
for(i in unique(data$id)){
  
  current <- NULL
  subject <- hold[hold$id==i,,drop=FALSE]
  
  #------------------------------------
  censoredAges <- censor_times_btsm(subject$time)
  
  for(tt in censoredAges ){
    
    # Rounding t, subject$time, & subject$disc_time to make sure we have equality
    t_round = round(tt, digits = 5)
    yrs_round = round(subject$time, digits = 5)
    disc_round = round(subject$disc_time, digits = 5)
    
    # If 't' corresponds to an observed age, then the next row will include the observed clinical visit data.
    if(t_round %in% yrs_round){
      current <- rbind( current, subject[disc_round==round(floor(tt), digits=5),])
    } else{
      
      # Create a CENSORED row for each subject at each discritezed time.
      tempRow['id'] <- i
      tempRow['time'] <- tt
      tempRow['Z.1'] <- subject$Z.1[1]
      tempRow['Z.2'] <- subject$Z.2[1]
      tempRow['state'] <- 99
      tempRow['obstrue'] <- 1
      tempRow['disc_time'] <- tt
      
      current <- rbind( current, tempRow)
      
      # If 't' corresponds to an observed INTEGER years, then the subject was observed some time during this years.  According, the next row will include the observed clinical visit data.  Recall that integer years is simply the floor(years).
      if(t_round %in% disc_round){ current <- rbind( current, subject[disc_round==t_round,])}
    }
    
  }
  #------------------------------------
  
  data_added_rows <- rbind( data_added_rows, current)
  num <- num+1
}
colnames(data_added_rows) <- colnames(hold)
rownames(data_added_rows) <- NULL

data.hmm = data_added_rows


bayestsm_dat = data.hmm

temp_data = as.matrix(bayestsm_dat); rownames(temp_data) = NULL



id = temp_data[, "id"]
y  = temp_data[, "state"]
x  = temp_data[, c("Z.1", "Z.2"), drop=F]
tt  = temp_data[, "time"]

disc_t = rep(-1, length(tt))
disc_t = temp_data[,"disc_time"]



# mcmc_chain = mcmc_routine(y, x, tt, id, init_par, prior_par, par_index,
#                           steps, burnin, disc, disc_t)

# mcmc_out$chain
# dim(mcmc_out$chain)
#############################################################################


out=list()
for(l in 1:n.chains){
  init_par = c(-3-2*runif(1), 0, 0, 0, # beta(0, x), beta(1, x), beta(2, x), beta(2, x_time)
               -3+runif(1), 0, 0, 0) # beta(0, t), beta(1, t), beta(2, t), beta(2, t_time)
  
  
  mcmc_chain = mcmc_routine(y, x, tt, id, init_par, prior_par, par_index,
                            steps, burnin, disc, disc_t)
  out[[l]]=  list(chain = mcmc_chain$chain, accept = mcmc_chain$accept)
  
}

mcmc.par = list()
mcmc.par.all = mcmc.par.all.est =  matrix()
accept.rate = NA

mat = list()

for(j in 1: length(out) ){
  
   chain = out[[j]]$chain
  
  chain.all = out[[j]]$chain
  
  mcmc.par[[j]] = mcmc(chain)
  
  accept.rate[j] = out[[j]]$accept
  ############################# 
  chain2 = chain
  colnames(chain2) = paste('pr',1:ncol(chain2), sep = '')
  mat[[j]] = chain2
  #############################  
  if(j==1){
    mcmc.par.all = chain
    
    mcmc.par.all.est  = chain.all
  }else{
    mcmc.par.all = rbind(mcmc.par.all, chain)
    mcmc.par.all.est = rbind(mcmc.par.all.est, chain.all)
  }
}
accept.rate
par.mcmc.list = mcmc.list(mcmc.par)
##########################################

myData = bayestsm_dat
# 

QQ <- function(x_ik,beta){
  
  betaMat = matrix(beta, ncol = 4, byrow = T)
  q_x  = exp(c(1,x_ik) %*% betaMat[1,])  # Transition from state 1 to state 2
  q_t  = exp(c(1,x_ik) %*% betaMat[2,])  # Transition from state 2 to state 3
  
  qmat = matrix(c( 0,q_x,  0,
                   0,  0,q_t,
                   0,  0,  0), nrow=3, byrow = T)
  diag(qmat) = -rowSums(qmat)
  
  return(qmat)
}




prob_evolution_expm <- function(state_from, state_to, z1, z2, par) {
  
  t_seq = NULL
  if(state_to == 2) {
    
    t_seq = t.x
  } else {
    
    t_seq = t.y
  }
  
  probEvo = rep(NA,length(t_seq))
  
  P = rep(0,3)
  P[state_from] = 1
  index = 1
  
  for(t in 2:length(t_seq)) {
    beta = par
    # Q_mat = Q(c(z1, z2, t_seq[t-1]), par)
    Q_mat = QQ(c(z1, z2, floor(t_seq[t-1])), par)
    new_P = expm((t_seq[t] - t_seq[t-1]) * Q_mat, method = 'Pade')
    P = P %*% new_P
    
    if(state_to == 2) {
      probEvo[index] = P[state_to] + P[state_to+1]
    } else {
      probEvo[index] = P[state_to]
    }
    
    index = index + 1
  }
  
  return(probEvo[-lx])
}

#
unique_ids = 1:n
cdf.x =  matrix(nrow =n, ncol = length(t.x)-1)
cdf.y = matrix(nrow =n, ncol = length(t.y)-1)

compute.cdf = function(para){
  par = para
  # par = mcmc.samples[para,]
  
  for(i in 1:nrow(cdf.x)) {
    # print(i)
    subDat = myData[myData[,'id'] == unique_ids[i], ]
    z_1 = subDat[1,'Z.1']
    z_2 = subDat[1,'Z.2']
    
    cdf.x[i,] = prob_evolution_expm(1, 2, z_1, z_2, par)
    
    cdf.y[i,] = prob_evolution_expm(2, 3, z_1, z_2, par)
    
  }
  
  cdx = colMeans(cdf.x)
  cdy = colMeans(cdf.y)
  return(c(cdx,cdy))
  
}




mcmc.samples = mcmc.par.all.est[sample(1:nrow(mcmc.par.all.est),n.post.samples),]



cdf.all = apply(mcmc.samples,1,compute.cdf)


p.x = rbind(0,cdf.all[1:((dim(cdf.all )[1])/2),])

p.y = rbind(0,cdf.all[(((dim(cdf.all )[1])/2)+1): dim(cdf.all )[1],])
#
cdf.x.med = apply(p.x,1,median,na.rm=T)
cdf.x.Lower = apply(p.x,1,quantile, c(0.025, 0.975),na.rm=T)[1,]
cdf.x.Upper = apply(p.x,1,quantile, c(0.025, 0.975),na.rm=T)[2,]


cdf.y.med = apply(p.y,1,median,na.rm=T)
cdf.y.Lower = apply(p.y,1,quantile, c(0.025, 0.975),na.rm=T)[1,]
cdf.y.Upper = apply(p.y,1,quantile, c(0.025, 0.975),na.rm=T)[2,]




est.x  = cdf.x.med
est.y  = cdf.y.med


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









