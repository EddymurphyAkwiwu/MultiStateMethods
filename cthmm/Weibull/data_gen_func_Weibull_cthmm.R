



##################  With 2 covariates: One normally distributed & one binary covariate ########## 


n = 1000  # Sample size

# n = 2000  # Sample size

p = 1 # number of covariates
p.discrete = 1
r=.1
s = 1
beta.X  = c(0.5,0.5)
beta.S  = c(0.5,0.5)

sigma.X = 0.2          # True scale parameter
mu.X    = 3           # True intercept parameter
sigma.S = 0.3        # True scale parameter
mu.S    = 1.2           # True intercept parameters
cor.X.S = 0

dist.X='weibull'
dist.S='weibull'

######  Strong censoring ###########
# # 
v.min   = 1           # Minimum time between screening moments
v.max   = 5.5           # Maximum time between screening moments
Tmax    = 2e2        # Maximum number of screening times (this is set to high value; but too high values increase computation time)
mean.rc = 25.9
do.b.XS = F; b.XS = NULL

#####  Medium censoring ########### 

# v.min   = 1           # Minimum time between screening moments
# v.max   = 7.3           # Maximum time between screening moments
# Tmax    = 2e2        # Maximum number of screening times (this is set to high value; but too high values increase computation time)
# mean.rc = 57.5
# do.b.XS = F; b.XS = NULL




if(p ==0){
  beta.X = as.matrix(mu.X)
  beta.S = as.matrix(mu.S)
}   else{
  beta.X = as.matrix(c(mu.X, beta.X))
  beta.S = as.matrix(c(mu.S, beta.S))
}

# Sim Z
R = matrix(r, p, p)
diag(R) =  1
S = rep(s, p)
Sigma = cor2cov(R, S)
if(p>0) {Z = mvrnorm(n, mu = rep(0,p), Sigma)} else Z = NULL
# Sim discrete Z
if(p.discrete == 1){
  Z.discrete = rbinom(n, 1, 0.5)
  Z = cbind( Z, Z.discrete)
  colnames(Z) = paste(1:ncol(Z))
}  
Z1 = cbind( as.matrix(rep(1,n)),Z)

# Sim.X
if( dist.X=='weibull' | dist.X=='loglog' | dist.X=='lognormal' ){
  if(dist.X=='weibull')    e.X = r.ev(n)
  if(dist.X=='loglog')     e.X = rlog(n)
  if(dist.X=='lognormal')  e.X = rnorm(n)
  X = Z1 %*% beta.X + sigma.X * e.X # log surv times
  X = exp(X) # surv times
  X  = as.numeric(X)
}
# Sim.S
if( dist.S=='weibull' | dist.S=='loglog' | dist.S=='lognormal' ){
  if(dist.S=='weibull')   e.S = r.ev(n)
  if(dist.S=='loglog')    e.S = rlog(n)
  if(dist.S=='lognormal') e.S = rnorm(n)
  S = Z1 %*% beta.S + sigma.S * e.S # log surv times
  if(do.b.XS){
    #logstdX = scale(X)
    S = Z1 %*% beta.S + b.XS * log(X) + sigma.S * e.S # log surv times
  } 
  S = exp(S) # surv times
  S = as.numeric(S)
}

if(dist.X=='bv-lognormal'){
  Sigma.e  = matrix( c(1, cor.X.S, cor.X.S, 1), 2, 2)
  e = mvrnorm(n, c(0,0), Sigma.e)
  X = Z1 %*% rbind(mu.X, beta.X) + sigma.X * e[,1] # log surv times
  S = Z1 %*% rbind(mu.S, beta.S) + sigma.S * e[,2] # log surv times
  X = as.numeric(exp(X))
  S = as.numeric(exp(S))
}

# Create total time
XS = as.numeric(X+S)  # Total time from baseline to event 2

# Generate screening sequences
t.rc = rexp(n, 1/mean.rc ) #runif(n, start.rc, vmax*(visitdist+leniancy))
V = as.matrix(runif(n, v.min, v.max ))
for(i in 2:Tmax){
  V = cbind(V, runif(n, V[,(i-1)]+v.min, V[,(i-1)]+v.max ))
}
rc = (V >= t.rc) * matrix( rep(t.rc, Tmax), ncol = Tmax, nrow = n)
V  = (V < t.rc) * V + rc 
V  = cbind(0, V)
V


# Find indices of events
ind.h  = rowSums(V < X) # Last healthy time index
ind.x  = ind.h+1        # First stage 1 
ind.xs = rowSums(V < XS) + 1 # First stage 2
ind.x[ind.x>ncol(V)] = ncol(V)

# Define delta (event)
d = rep(-99,n)
d[ind.x < ind.xs ]  = 2
d[ind.x == ind.xs  & ind.xs != (ncol(V)+1)] = 3
d[ind.h == ind.x]  = 1

# Define left and right interval bound
L = V[ cbind(1:nrow(V), ind.h) ]
R = V[ cbind(1:nrow(V), ind.x) ]
R[ ind.h == ind.x ] = Inf

if(p>0) dat = data.frame(L=L, R=R, d=d, X=X, S=S, Z=Z)
if(p==0) dat = data.frame(L=L, R=R, d=d, X=X, S=S)


dat2 = list()  

for(k in 1:nrow(V)){
  
  
  dd = dat$d[k]
  if(dd==2 ){
    obs.times = V[ k, 1: ind.x[k] ]
    obs.data = rep(1, length(obs.times)) 
    obs.data[length(obs.times)] = 2
  }else if(dd==3){
    obs.times = V[ k, 1: ind.x[k] ]
    obs.data = rep(1, length(obs.times)) 
    obs.data[length(obs.times)] = 3
  }else{
    obs.times = unique(V[ k, 1: ind.x[k] ])
    obs.data = rep(1, length(obs.times))
  }
  
  
  dt =  list(obs.times=obs.times, obs.data = obs.data)
  
  dat2[[k]]   = dt
  
}  

if(p==0) cov.data = data.frame(id = 1:n)
if(p>0) cov.data = data.frame(id = 1:n, Z=Z)

