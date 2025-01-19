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

rm(list=setdiff(ls(), c("long_data", "N", "pi","gamma","hidden_states")))

# Now simulate from the Multinomial distribution
D = 4 # Claims occurring in t can be reported in t+d (d=0,1,2,3,4)

# Define targeted p_t values for x2
p_A <- rep(1,D+1)/sum(rep(1,D+1))
p_C <- (D+1):1/sum((D+1):1)
p_B <- (p_A + p_C)/2

# Define q_t values for x2 
q_target <- function(D, p){
  q = rep(0,D+1)
  q[1] = 1
  q[D+1] = p[D+1]
  for(t in 2:D){
    q[t] = p[t]/sum(p[1:t])
  }
  return(q)
}
q_A = q_target(D,p_A)
q_B = q_target(D,p_B)
q_C = q_target(D,p_C)

# log(q_t) = delta[1] + delta[2] x 1{x2 = B} + delta[3] x 1{x2 = C} regression coefficients
delta <- matrix(0, nrow = 3, ncol = D)
delta[1,] <- log(q_A[2:(D+1)])
delta[2,] <- log(q_B[2:(D+1)]) - delta[1,]
delta[3,] <- log(q_C[2:(D+1)]) - delta[1,]

# Simulate Z_{i,t,d}
Z <- matrix(0, nrow = nrow(long_data), ncol = D + 1)
for (i in 1:nrow(long_data)) {
  prob = if(long_data$x2[i] == "A"){
    prob = p_A
  } else{
    if(long_data$x2[i] == "B"){
      prob = p_B
    } else {
      prob = p_C
    }
  }
  Z[i,] <- rmultinom(1, long_data$count[i], prob = prob)
}
long_data <- cbind(long_data,Z)
colnames(long_data) <- c("ID","time","count","exposure","x1","x2","x3",c(paste0("d",0:D)))

# Censor data at time = N
censored_data <- long_data
for(i in 1:nrow(censored_data)){
  if(censored_data$time[i] > N-D){
    censored_data[i,"count"] <- censored_data[i,"count"] - sum(censored_data[i,c(paste("d",(N-censored_data$time[i]+1):D, sep=""))])
    censored_data[i,c(paste("d",(N-censored_data$time[i]+1):D, sep=""))] <- NA
  }
}

sim.parameters <- list(pi = pi, gamma = gamma, delta=delta,hidden_states = hidden_states)
rm(list=setdiff(ls(), c("censored_data", "N","D","sim.parameters")))

###################################################################################################
## EM Algorithm
# This function calculate sum(log.prob) for each period and for all lambda^{j}. 
lprobmatrix <- function(long.data, lambda){
  
  # lambda is nrow(long.data) x g 
  time = long.data$time
  lp <- 1:length(unique(long.data$time))
  
  for(j in 1:ncol(lambda)){
    temp <- long.data
    temp$lambda <- lambda[,j]
    temp <- temp %>% mutate(lprob = dpois(count, ofst * lambda,
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

# Prob vector estimate
prob.vec <- function(long.data, longer.data, D){
  
  # for i = 1
  temp.data <- longer.data %>% filter(delay.period == "d1")
  temp.data <- temp.data %>% filter(cum.sum != 0)
  reg <- glm(cbind(freq.counts, cum.sum - freq.counts) ~x2,
             data = temp.data, family = quasibinomial("log"))
  pred <- predict(reg, newdata = long.data, type = "response")
  delta <- reg$coefficients
  
  # for i = 2 to D
  for(i in 2:D){
    temp.data <- longer.data %>% filter(delay.period == paste("d",i,sep=""))
    temp.data <- temp.data %>% filter(cum.sum != 0)
    reg <- glm(cbind(freq.counts, cum.sum - freq.counts) ~ x2, 
               data = temp.data, family = quasibinomial("log"))
    pred <- cbind(pred, predict(reg, newdata = long.data, type = "response"))
    delta <- cbind(delta, reg$coefficients)
  }
  delta <- as.data.frame(delta)
  colnames(delta) <- c(paste("E",1:D,sep=""))
  
  prob.est <- as.data.frame(matrix(0, nrow = nrow(pred), ncol= D+1))
  prob.est[,D+1] <- pred[,D]
  prob.est[,D] <- pred[,D-1]*(1-prob.est[,D+1])
  for(i in (D-1):2){
    prob.est[,i] <- pred[,i-1]*(1-rowSums(prob.est[,(i+1):(D+1)]))
  }
  prob.est[,1] <- 1-rowSums(prob.est[,2:(D+1)])
  prob.est <- prob.est/rowSums(prob.est)
  return(list(prob.init = prob.est, delta = delta))
}
# EM estimation
pois.HMM.EM <- function(n, g, long.data, log.probs, theta,  
                        gamma, pi, delta, lambda, delay, D, maxiter =1000, tol=10^(-6) ,...)
{
  ncoef <- nrow(theta) 
  theta.next <- theta
  gamma.next <- gamma
  pi.next <- pi
  delta.next <- delta
  
  # Arrange data for binomial reg
  longer.data <- gather(long.data,delay.period,freq.counts, d0:paste("d",D,sep=""), factor_key=TRUE)
  z_r <- longer.data$freq.counts
  start = Sys.time()
  ## itteration
  for (iter in 1: maxiter){
    
    # obtain forward and backward probabilities and the likelihood
    fb <- pois.HMM.lalphabeta(t = n, log.prob.matrix =log.probs, 
                              g = g, gamma = gamma, pi=pi)
    lalpha <- fb$la 
    lbeta <- fb$lb 
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
    
    ## E-step 
    
    # next thetas 
    ## find full N_{i,t}
    temp <- data.frame(time = long.data$time, count = long.data$count, 
                       rep.prob = long.data$rep.prob,
                       exposure = long.data$exposure)
    temp <- cbind(temp, lambda)
    for (i in 1:g) {
      temp[[paste0("freq_", i)]] <- temp$count + temp[,paste0("L", i)]*temp$exposure*(1-temp$rep.prob)
    }
    
    sum_u = colSums(exp(lalpha+lbeta-llik))
    
    lambda.temp <- lambda
    
    for(j in 1:g){
      long.data$freq.adj <- temp[,paste("freq_",j,sep="")]
      u = exp(lalpha[j,]+lbeta[j,]-llik)/sum_u
      w = u[long.data$time]
      mod = glm(freq.adj ~ x1+x2+x3, offset = log(exposure),
                data = long.data, family = quasipoisson("log"), weights = w)
      theta.next[,j] = as.numeric(mod$coefficients)
      lambda[,j] = exp(predict(mod, newdata = long.data, type = "link"))/long.data$exposure
    }
    rm(temp,mod)
    
    ## next delta's
    u = exp(lalpha+lbeta-llik)/colSums(exp(lalpha+lbeta-llik))
    w = u[,long.data$time]
    weighted.lambda <- data.frame(L1 = w[1,]*lambda.temp[,1],
                                  L2 = w[2,]*lambda.temp[,2]) 
    weighted.lambda <- rowSums(weighted.lambda)
    weighted.lambda <- weighted.lambda * long.data$exposure
    z <- delay * weighted.lambda
    colnames(z) <- paste("d",0:D,sep="")
    z <- gather(z, delay.period, freq.counts, d0:paste("d",D,sep=""), factor_key=TRUE)
    z <- z$freq.counts
    rm(u,w, weighted.lambda, lambda.temp)
    z[which(is.na(z_r) != T)] <- z_r[which(is.na(z_r) != T)]
    longer.data$freq.counts <- z
    longer.data <- longer.data %>% group_by(ID, time) %>%
      mutate(cum.sum = cumsum(freq.counts))
    rm(z)
    
    x <- prob.vec(long.data = long.data, longer.data = longer.data, D)
    delta.next <- x$delta
    delay <- x$prob.init
    rm(x)
    
    # checking convergence 
    crit <- sum(abs(theta -theta.next)) + sum(abs(gamma -gamma.next)) +
      sum(abs(pi -pi.next)) + sum(abs(delta -delta.next))
    
    # next itteration
    theta <- theta.next
    gamma <- gamma.next
    pi <- pi.next
    delta <- delta.next 
    
    # calc prob matrix
    cumsum.delay <- t(apply(delay, 1, cumsum))
    time <- long.data$time
    rep.prob <- rep(0, length(time))
    for(i in 1:length(rep.prob)){
      rep.prob[i] <- ifelse(time[i] <= n - D, 1, cumsum.delay[i, n - time[i] + 1])
    }
    long.data$rep.prob <- rep.prob
    long.data$ofst <- long.data$exposure * long.data$rep.prob
    rm(time, rep.prob)
    log.probs <- lprobmatrix(long.data, lambda)
    
    if(crit <tol){
      end = Sys.time()
      np <- g*g+g*ncoef-1 
      BIC <- -2*llik+np*log(n)
      return(list(BIC = BIC, theta=theta, pi = pi, gamma=gamma, delta=delta, lambda=lambda, delay=delay,
                  run.time = (end-start)/iter, iter=iter))
    }
    
  }
  print(paste ("No convergence after ",maxiter ," iterations "))
  NA
}
###################################################################################################

truncated_data <- censored_data %>% dplyr::filter(time <= N-D) #consider periods where we know everything
## Initialization 
g= 2 #number of states

# pi (initial distribution)
pi.init <- rep(1,g)/g #Initial distribution

# delta (regression coefficients for q_t's) .. initialize delta[2] = delta[3] = 0 
# delta[1] is initialized based on empirical data
p  = 3 #no. of regression coefficients
delta.init <- matrix(0, nrow = p, ncol = D)
sum_d <- colSums(truncated_data[,c(paste("d",0:D,sep = ""))])
cumsum_d <- cumsum(sum_d)
for(t in 1:D){
  delta.init[1,t] <- log(sum_d[t+1]/cumsum_d[t+1])
}

# Calculate probability of reporting based on initialization of delta (i.e) \sum_{d=0}^{min{D,N}} p_t(d;x_i)
rep.prob <- cumsum(sum_d/sum(sum_d))
censored_data$rep.prob <- 1
for(i in 1:nrow(censored_data)){
  if(censored_data$time[i] > N-D){
    censored_data$rep.prob[i] <- rep.prob[N-censored_data$time[i]+1]
  }
}

# initial p_t(d;x_i)
p_t <- sum_d/sum(sum_d)
delay.init <- do.call(rbind, replicate(nrow(censored_data), p_t, simplify = F)) # initial lambda for each ID & g
delay.init <- data.frame(delay.init)
rm(rep.prob,sum_d,cumsum_d,p_t)

# theta (regression coefficients for lambda)
# Assume all coefficients are zero except for the intercept (initial value obtained from empirical data)
p  = 5 #no. of regression coefficients
theta.init <- matrix(0, nrow = p, ncol = g)
avg_lambda <- truncated_data %>% group_by(time) %>%
  summarise(lambda = sum(count)/sum(exposure))
cluster_id <- kmeans(avg_lambda$lambda, centers = sort(kmeans(avg_lambda$lambda, centers = g)$centers))$cluster
avg_lambda$cluster_id <- cluster_id
avg_lambda <- avg_lambda %>% group_by(cluster_id) %>%
  summarise(lambda = mean(lambda))
theta.init[1,] <- log(avg_lambda$lambda) #initial regression coefficients

# gamma (transition matrix)
gamma.init <- markovchainFit(cluster_id) #Transition matrix
gamma.init <- gamma.init$estimate[1:g]

# lambda (initial estimates of lambdas based on theta.init)
lambda.init <- matrix(0, nrow = 1, ncol = g)
lambda.init[1,] <- avg_lambda$lambda
lambda.init <- do.call(rbind, replicate(nrow(censored_data), lambda.init, simplify = F)) # initial lambda for each ID & g
lambda.init <- as.data.frame(lambda.init)
colnames(lambda.init) <- c(paste("L",1:g,sep=""))
rm(cluster_id, avg_lambda,p)

#######################################################################################################
## Running the model
censored_data <- censored_data %>% mutate(ofst = exposure*rep.prob)
log.probs <- lprobmatrix(censored_data, lambda.init)
run.model <- pois.HMM.EM(n = N,g=g,long.data = censored_data, 
                         log.probs = log.probs, theta = theta.init,
                         gamma = gamma.init, delta = delta.init,
                         pi = pi.init, lambda = lambda.init, delay = delay.init, 
                         D = D, tol = 0.01)
print(run.model$run.time)


####
print(run.model$pi)
print(run.model$gamma)
print(run.model$delta)
print(run.model$theta)
print(run.model$run.time)
print(run.model$iter)

q_A <- exp(run.model$delta[1,])
q_B <- exp(run.model$delta[1,]+run.model$delta[2,])
q_C <- exp(run.model$delta[1,]+run.model$delta[3,])

p_fn <- function(q){
  q <- as.numeric(q)
  p <- rep(0,5)
  p[5] = q[4]
  s = sum(p)
  for(i in 4:2){
    p[i] <- (1-s)*q[i-1]
    s = sum(p)
  }
  p[1] <- 1-s
  return(p)
}
print(p_fn(q_A))
print(p_fn(q_B))
print(p_fn(q_C))
