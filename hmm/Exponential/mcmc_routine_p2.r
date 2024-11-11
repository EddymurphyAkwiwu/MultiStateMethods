library(mvtnorm, quietly=T);
library(Rcpp, quietly=T)
library(RcppArmadillo, quietly = T)
library(RcppDist, quietly = T)
# library(deSolve)
library(foreach)
library(doParallel, quietly=T)
# sourceCpp("likelihood.cpp")



# -----------------------------------------------------------------------------
# The mcmc routine for sampling the parameters
# -----------------------------------------------------------------------------



mcmc_routine = function( y, x, t, id, init_par, prior_par, par_index, steps, burnin, disc=TRUE, disc_t) {
  
  
  disc=TRUE
  pars = init_par
  n = length(y)
  n_par = length(pars)
  chain = matrix( 0, steps, n_par)
  
  # group = list(c(par_index$beta))
  group = list(c(par_index$beta[-c(4 ,8)]))  #  with covaraite; exp distribution 

  n_group = length(group)
  
  pcov = list();	for(j in 1:n_group)  pcov[[j]] = diag(length(group[[j]]))
  pscale = rep( .0001, n_group)
  
  accept = rep( 0, n_group)
  eids = unique(id)
  
  # Evaluate the log posterior of the initial parameters
  
  log_post_prev = eddy::fn_log_post( pars, prior_par, par_index, x, y, t, id, eids, disc_t)
  
  if(!is.finite(log_post_prev)){
    print("Infinite log-posterior; choose better initial parameters")
    break
  }
  
  # Begin the MCMC algorithm --------------------------------------------------
  chain[1,] = pars
  for(ttt in 2:steps){
    for(j in 1:n_group){
      
      # Propose an update
      ind_j = group[[j]]
      proposal = pars
      proposal[ind_j] = rmvnorm( n=1, mean=pars[ind_j],sigma=pcov[[j]]*pscale[j])
      
      # Compute the log density for the proposal
      
      log_post = eddy::fn_log_post(proposal, prior_par, par_index, x, y, t, id, eids, disc_t)   
      
      
      # Only propose valid parameters during the burnin period
      if(ttt < burnin){
        while(!is.finite(log_post)){
          print('bad proposal')
          proposal = pars
          proposal[ind_j] = rmvnorm( n=1, mean=pars[ind_j],
                                     sigma=pcov[[j]]*pscale[j])
          
          log_post = eddy::fn_log_post(proposal, prior_par, par_index, x, y, t, id, eids, disc_t)   
          
          
        }
      }
      
      # Evaluate the Metropolis-Hastings ratio
      if( log_post - log_post_prev > log(runif(1,0,1)) ){
        log_post_prev = log_post
        pars[ind_j] = proposal[ind_j]
        accept[j] = accept[j] +1
      }
      chain[ttt,ind_j] = pars[ind_j]
      
      # Proposal tuning scheme ------------------------------------------------
      if(ttt < burnin){
        # During the burnin period, update the proposal covariance in each step
        # to capture the relationships within the parameters vectors for each
        # transition.  This helps with mixing.
        if(ttt == 100)  pscale[j] = 1
        
        if(100 <= ttt & ttt <= 2000){
          temp_chain = chain[1:ttt,ind_j]
          pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
          
        } else if(2000 < ttt){
          temp_chain = chain[(ttt-2000):ttt,ind_j]
          pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
        }
        if( sum( is.na(pcov[[j]]) ) > 0)  pcov[[j]] = diag( length(ind_j) )
        
        # Tune the proposal covariance for each transition to achieve
        # reasonable acceptance ratios.
        if(ttt %% 30 == 0){
          if(ttt %% 480 == 0){
            accept[j] = 0
            
          } else if( accept[j] / (ttt %% 480) < .4 ){ 
            pscale[j] = (.75^2)*pscale[j]
            
          } else if( accept[j] / (ttt %% 480) > .5 ){ 
            pscale[j] = (1.25^2)*pscale[j]
          }
        }
      }
      # -----------------------------------------------------------------------
    }
    # Restart the acceptance ratio at burnin.
    if(ttt == burnin)  accept = rep( 0, n_group)
    
    # if(ttt%%1==0)  cat('--->',ttt,'\n')
  }
  # ---------------------------------------------------------------------------
  
  # print(accept/(steps-burnin))
  return(list( chain=chain[burnin:steps,], accept=accept/(steps-burnin),
               pscale=pscale, pcov = pcov))
}
# -----------------------------------------------------------------------------

