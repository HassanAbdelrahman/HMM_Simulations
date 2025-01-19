library(dplyr)
library(tidyverse)
library(lubridate)
library(markovchain)
library(matrixStats)

set.seed(12345) # For reproducibility

# Number of individuals and time points
n_individuals <- 1000 
N <- 100

# Covariates
x1 <- sample(1:50, n_individuals, replace = TRUE)
x2 <- sample(c("A", "B", "C"), n_individuals, replace = TRUE)
x3 <- sample(c("X", "Y"), n_individuals, replace = TRUE)

# Encode categorical variables as indicators
x2_B <- as.numeric(x2 == "B")
x2_C <- as.numeric(x2 == "C")
x3_Y <- as.numeric(x3 == "Y")

# Covariates
covariates <- data.frame(ID = paste0(x1,x2,x3), x1, x2, x3)
covariates <- unique(covariates)

# Hidden Markov process parameters
pi <- c(1, 0) # Initial distribution
gamma <- matrix(c(0.7, 0.3,
                  0.4, 0.6),
                nrow = 2, byrow = TRUE) # Transition matrix

# Log-intensity functions
log_lambda <- function(state, x1, x2_B, x2_C, x3_Y) {
  if (state == 1) {
    return(0.20 + 0.0010 * x1 + 0.05 * x2_B + 0.10 * x2_C - 0.05 * x3_Y)
  } else if (state == 2) {
    return(0.40 + 0.0020 * x1 + 0.10 * x2_B + 0.15 * x2_C - 0.05 * x3_Y)
  } 
}

# Simulate hidden states for N time points
simulate_hidden_states <- function(N, pi, gamma) {
  states <- numeric(N)
  states[1] <- sample(1:2, 1, prob = pi)
  for (t in 2:N) {
    states[t] <- sample(1:2, 1, prob = gamma[states[t - 1], ])
  }
  return(states)
}

# Simulate Poisson counts
simulate_counts <- function(hidden_states, x1, x2_B, x2_C, x3_Y) {
  counts <- rep(0,length(hidden_states))
  for (t in 1:length(counts)) {
    lambda_log <- log_lambda(hidden_states[t], x1, x2_B, x2_C, x3_Y)
    lambda <- exp(lambda_log)
    counts[t] <- rpois(1, lambda)
  }
  return(counts)
}

# Simulate data for all individuals
# Simulate wide-format data
simulated_data <- data.frame(x1 = x1, x2 = x2, x3 = x3)
for (t in 1:N) {
  simulated_data[[paste0("count_t", t)]] <- NA
}

# Fill in counts
hidden_states <- simulate_hidden_states(N, pi, gamma)

for (i in 1:n_individuals) {
  xi1 <- x1[i]
  xi2_B <- as.numeric(x2[i] == "B")
  xi2_C <- as.numeric(x2[i] == "C")
  xi3_Y <- as.numeric(x3[i] == "Y")
  
  counts <- simulate_counts(hidden_states, xi1, xi2_B, xi2_C, xi3_Y)
  
  simulated_data[i, 4:(N + 3)] <- counts
}

# Convert to long format using pivot_longer
long_data <- simulated_data %>%
  pivot_longer(cols = starts_with("count_t"), 
               names_to = "time", 
               names_prefix = "count_t", 
               values_to = "count") %>%
  mutate(time = as.integer(time)) # Ensure time is an integer

# Combine individuals with the same covariates
long_data <- long_data %>% mutate(ID = paste0(x1,x2,x3)) %>%
  group_by(ID, time) %>% summarise(count = sum(count), exposure = n())
long_data <- long_data %>% left_join(covariates, by = "ID")

rm(list=setdiff(ls(), c("long_data", "N")))

# EM Algorithm
# This function calculate sum(log.prob) for each period and for all lambda^{j}. 
lprobmatrix <- function(data, lambda){
  
  # lambda's dimension is nrow(data) x g 
  time = data$time
  lp <- 1:length(unique(data$time))
  for(j in 1:ncol(lambda)){
    temp <- data
    temp$lambda <- lambda[,j]
    temp <- temp %>% mutate(lprob = dpois(count, exposure * lambda,
                                          log = TRUE)) %>%
      group_by(time) %>% summarise(sum.lp = sum(lprob)) 
    lp <- cbind(lp, temp$sum.lp)
  }
  
  return(lp[,-1]) # T x g matrix
}

# Function for the forward and backward probabilities
pois.HMM.lalphabeta <-function(t, log.prob.matrix, g, gamma, pi=NULL)
{
  # log.prob.matrix is T x g matrix .. gamma is gxg matrix
  if(is.null(pi)) pi <-solve(t(diag(g)-gamma +1),rep(1,g))
  
  lalpha <- lbeta <- matrix(NA ,g,t)
  
  lalpha[,1] <- log(pi)+log.prob.matrix[1,]
  
  for (n in 2:t){
    log_m = max(lalpha[,n-1])
    sum_exp=rep(0,g)
    for(j in 1:g){
      sum_exp = sum_exp + exp(lalpha[j,n-1]-log_m+log(gamma[j,]))
    }
    lalpha[,n] <- log_m + log.prob.matrix[n,] + log(sum_exp)
  }
  
  llik = logSumExp(lalpha[,t])
  
  # c <- max(lalpha[,n])
  # llik <- c+log(sum(exp(lalpha[,n]-c)))
  
  lbeta[,t] <- rep(0,g)
  
  for (n in (t-1) :1){
    log_m = max(lbeta[,n+1]+log.prob.matrix[n+1,])
    sum_exp=rep(0,g)
    for(j in 1:g){
      sum_exp = sum_exp + exp(lbeta[j,n+1]-log_m+log(gamma[,j])+ log.prob.matrix[n+1,j])
    }
    lbeta[,n] <- log_m + log(sum_exp)
  }
  
  list(la=lalpha ,lb=lbeta, loglik = llik)
}

# EM estimation
pois.HMM.EM <- function(n, g, data, log.probs, theta,  
                        gamma, pi, maxiter =1000, tol=10^(-6) ,...)
{
  ncoef <- nrow(theta)
  theta.next <- theta
  gamma.next <- gamma
  pi.next <- pi
  start = Sys.time()
  for (iter in 1: maxiter){
    
    # obtain forward and backward probabilities and the likelihood
    fb <- pois.HMM.lalphabeta(t = n, log.prob.matrix =log.probs, 
                              g = g, gamma = gamma, pi=pi)
    lalpha <- fb$la # g x T matrix
    lbeta <- fb$lb # g x T matrix
    llik <- fb$loglik
    
    # next gamma
    for (j in 1:g){
      for (k in 1:g){
        gamma.next[j,k] <- gamma[j,k]*sum(exp(lalpha[j ,1:(n-1)]+log.probs[2:n,k]
                                              + lbeta[k,2:n] -llik))
      }
    }
    gamma.next <- gamma.next/apply(gamma.next ,1,sum)
    
    # next pi
    pi.next <- exp(lalpha[,1]+lbeta[,1]-llik)
    pi.next <- pi.next/sum(pi.next)
    
    # next thetas 
    sum_u = colSums(exp(lalpha+lbeta-llik))
    lambda <- as.data.frame(matrix(0, nrow = nrow(data), ncol = g))
    for(j in 1:g){
      u = exp(lalpha[j,]+lbeta[j,]-llik)/sum_u
      w = u[data$time]
      mod = glm(count ~ x1 + x2 + x3, offset = log(exposure),
                     data = data, family = poisson("log"), weights = w)
      theta.next[,j] = as.numeric(mod$coefficients)
      lambda[,j] = exp(predict(mod, newdata = data, type = "link"))/data$exposure
    }
    
    # checking convergence 
    crit <- sum(abs(theta -theta.next)) + sum(abs(gamma -gamma.next)) +
      sum(abs(pi -pi.next))
    # next itteration
    theta <- theta.next
    gamma <- gamma.next
    pi <- pi.next
    
    # calc prob matrix
    log.probs <- lprobmatrix(data, lambda)
    
    if(crit <tol){
      end = Sys.time()
      np <- g*g+g*ncoef-1 
      BIC <- -2*llik+np*log(n)
      return(list(BIC=BIC,theta=theta,gamma=gamma,pi=pi,run.time = (end-start)/iter, iter=iter))
    }
  }
  print(paste ("No convergence after ",maxiter ," iterations "))
  NA
}

############################################################################################

## Initialization
g = 2 # this can be changed to any number of states

# pi 
pi.init <- rep(1,g)/g #Initial distribution

# theta 
p  = 5 #no. of regression coefficients
theta.init <- matrix(0, nrow = p, ncol = g)
avg_lambda <- long_data %>% group_by(time) %>%
  summarise(lambda = sum(count)/sum(exposure))
cluster_id <- kmeans(avg_lambda$lambda, centers = sort(kmeans(avg_lambda$lambda, centers = g)$centers))$cluster
avg_lambda$cluster_id <- cluster_id
avg_lambda <- avg_lambda %>% group_by(cluster_id) %>%
  summarise(lambda = mean(lambda))
theta.init[1,] <- log(avg_lambda$lambda) #initial regression coefficients

# gamma
gamma.init <- markovchainFit(cluster_id) #Transition matrix
gamma.init <- gamma.init$estimate[1:g]

# lambda
lambda.init <- matrix(0, nrow = 1, ncol = g)
lambda.init[1,] <- avg_lambda$lambda
lambda.init <- do.call(rbind, replicate(nrow(long_data), lambda.init, simplify = F)) # initial lambda for each ID & g
rm(cluster_id, avg_lambda)

# log probibility matrix P
prob.init <- lprobmatrix(long_data, lambda.init)

# running the model
run.model <- pois.HMM.EM(n = N, g=g, data = long_data, 
                         log.probs = prob.init, theta = theta.init,
                         gamma = gamma.init, pi = pi.init)


# print results
print(run.model$BIC)
print(run.model$theta)
print(run.model$gamma)
print(run.model$pi)
print(run.model$iter)
print(run.model$run.time)
