
# rm(list = ls())       
library(msm)         
library(expm)         

# Load the data generated from an exponentially distribution progression times in a 3-state model
### Note that data is without covariates. Models with covariates follow similar approach below
load('data.RData')    

# Display the first few rows of data to inspect
head(data)



### Model Setup ###

# Define a 3-state Markov model with the transition intensity matrix Q
# Here, each entry specifies the transition rates between states
# Note: 0.01 is used as a starting value for non-zero transitions; can be adjusted 
Q = rbind(
  c(0, 0.01, 0),
  c(0, 0, 0.01),
  c(0, 0, 0)
)

# Fit the model to the data
msm.mod = msm(state ~ time, subject = id, data = data, qmatrix = Q)


### CDF Estimation  using three different approaches/methods ###

# Example time point for CDF estimation
t.time = 3


## #####################  Method 1: #####################################################
#      Using the  passage probability function `ppass.msm` 
#############################################################################################

#  Refer to help(ppass.msm) in msm package for more details
#  # See section 3.7 of  https://chjackson.github.io/msm/msmcourse/outputs.html#passage-probabilities 

CDF.method_one = ppass.msm(msm.mod, tot=t.time)

# Extract the CDF for the transition from state 1 to state 2
CDF.method_one[1, 2]

# Extract the CDF for the transition from state 2 to state 3
CDF.method_one[2, 3]



## #####################  Method 2: #####################################################
#      Using the Estimated Transition Intensity Matrix
#############################################################################################
 # # See help(qmatrix.msm) of the msm package and 
 ## section 7.1.1 of the cthmn package manual by Lange et (2013)

# The sojourn time in state A has the same distribution as time to absorption in a process that
# starts in state A and treats state B as an absorbing state. Thus, we need to construct a new transition intensity matrix
# from the original transition intensity matrix Q.est by setting the rows of the original matrix to zero for the state
# corresponding to the  state B.


# Step 1: Extract the estimated transition intensity matrix from the fitted model
Q.est = qmatrix.msm(msm.mod)[[1]]  # [[1]] extracts only the estimated rates, excluding confidence intervals

# Step 2: Define and modify transition intensity matrix Q.est for state 1 to state 2 transition
# For state 1 to state 2, set all transitions out of state 2 to zero, treating it as an absorbing state

Q.est.new = Q.est        # Copy the estimated Q matrix
Q.est.new[2, ] = 0       # Set second row to zero to make state 2 an absorbing state

# Step 3: Compute transition probability matrix at t = 3 using matrix exponential
CDF.method_two = expm(Q.est.new * t.time)

# Extract CDF for the transition from state 1 to state 2
CDF.method_two[1, 2]

# Repeat process to estimate time distribution from state 2 to state 3
Q.est.new = Q.est         # Reset to original Q matrix
Q.est.new[3, ] = 0        # Set third row to zero to make state 3 an absorbing state

# Compute the CDF at t = 3 for this modified matrix
CDF.method_two = expm(Q.est.new * t.time)

# Extract CDF for the transition from state 2 to state 3
CDF.method_two[2, 3]



## #####################  Method 3: #####################################################
#                 Using the Estimated Transition Probability Matrix
#############################################################################################
### See help(pmatrix.msm) of the msm package 

# Compute transition probability matrix directly at t = 3 using pmatrix.msm
CDF.method_third = pmatrix.msm(msm.mod, t=t.time)

# CDF from state 1 to state 2 is the sum of probabilities reaching either state 2 or state 3 from state 1
CDF.method_third[1, 2] + CDF.method_third[1, 3]

# CDF from state 2 to state 3
CDF.method_third[2, 3]


######################### Plotting Estimated CDFs vs True Empirical CDFs #########################

# Load assumed true CDF data for comparison in a large population setting
load('ExpTrueCDF_SCen_p_0.RData')    # Assumed true CDF in a  large sample setting

# Define time points for plotting
t.x = seq(0, 100, length.out=100)    # Time for state 1 to state 2 transition
t.y = seq(0, 20, length.out=100)     # Time for state 2 to state 3 transition

# Define empirical cumulative distribution functions from the data
F1 = ecdf(dat$X)                     # Empirical CDF for state 1 to 2 transition
F2 = ecdf(dat$S)                     # Empirical CDF for state 2 to 3 transition

# Generate true CDF values at specified time points
TrueX = F1(t.x)
TrueY = F2(t.y)


# Estimate CDFs using pmatrix.msm for a range of time points
cdf.x = cdf.y = numeric(length(t.x)) # Initialize vectors for storing estimated CDFs

for (i in 1:length(t.x)) {
  prob.x = pmatrix.msm(msm.mod, t=t.x[i])
  prob.y = pmatrix.msm(msm.mod, t=t.y[i])
  
  cdf.x[i] = prob.x[1, 2] + prob.x[1, 3]  # CDF from state 1 to 2
  cdf.y[i] = prob.y[2, 3]                # CDF from state 2 to 3
}


# Plot estimated and true CDFs for state 1 to 2 transition
plot(t.x, cdf.x, type="l", col=1, lty=1, lwd=2, ylim=c(0, 1),
     xlab="time since entry in state 1", ylab="CDF",
     main="Transition from state 1 to 2")
lines(t.x, TrueX, type="l", lwd=2, col='blue', lty=2)   # True empirical CDF (data)

# Add legend
legend("bottomright", legend=c("Model Estimate", "True"),
       col=c("black", "blue"), lty=c(1, 2), lwd=2)


# Plot estimated and true CDFs for state 2 to 3 transition
plot(t.y, cdf.y, type="l", col=1, lty=1, lwd=2, ylim=c(0, 1),
     xlab="time since entry in state 2", ylab="CDF",
     main="Transition from state 2 to 3")
lines(t.y, TrueY, type="l", lwd=2, col='blue', lty=2)   # True empirical CDF (data)

# Add legend
legend("bottomright", legend=c("Model Estimate", "True"),
       col=c("black", "blue"), lty=c(1, 2), lwd=2)
