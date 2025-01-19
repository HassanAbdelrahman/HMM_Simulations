library(dplyr)
library(tidyverse)
library(lubridate)
library(markovchain)
library(matrixStats)
library(DirichletReg)
library(foreach)
library(doParallel)

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

# Delta's (Dirichlet regression parameters)
eta_0 = 10000
eta_A <- eta_0 * p_A 
eta_B <- eta_0 * p_B
eta_C <- eta_0 * p_C
delta <- rbind(log(eta_A),log(eta_B)-log(eta_A),log(eta_C)-log(eta_A))

# Simulate Z_{i,t,d}
Z <- matrix(0, nrow = nrow(long_data), ncol = D + 1)
for (i in 1:nrow(long_data)) {
  prob = if(long_data$x2[i] == "A"){
    eta = eta_A
  } else{
    if(long_data$x2[i] == "B"){
      eta = eta_B
    } else {
      eta = eta_C
    }
  }
  if(long_data$count[i] > 1){
    prob <- rdirichlet(long_data$count[i], eta)
    for(t in 1:nrow(prob)){
      Z[i,] <- Z[i,] + rmultinom(1, 1, prob = prob[t,])
    }
  }
  
}

long_data <- cbind(long_data, data.frame(Z))
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
lprobmatrix <- function(long.data, lambda, alpha, N, D){
  
  # lambda is nrow(long.data) x g 
  time = long.data$time
  lp <- 1:length(unique(long.data$time))
  
  index <- which((time > (N-D)) & (time <N))
  index_n <- which(time ==N)
  
  for(j in 1:ncol(lambda)){
    temp <- long.data
    temp$lambda <- lambda[,j]
    temp <- temp %>% mutate(lprob = dpois(count, exposure * lambda,
                                          log = TRUE))
    for(k in index){
      u <- rdirichlet(1000, alpha = alpha[k,])
      u <- rowSums(u[,1:(N-time[k]+1)])
      temp$lprob[k] <- mean(dpois(temp$count[k], temp$exposure[k] * temp$lambda[k] * u,
                                  log = T))
      
    }
    
    for(k in index_n){
      u <- rdirichlet(1000, alpha = alpha[k,])
      u <- u[,1]
      temp$lprob[k] <- mean(dpois(temp$count[k], temp$exposure[k] * temp$lambda[k] * u,
                                  log = T))
      
    }
    temp <- temp %>% group_by(time) %>% summarise(sum.lp = sum(lprob)) 
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

# Acceptance Rejection
acc_rej <- function(p_r, pi, int, n_rep, p_min=0) {
  
  # Ensure inputs are numeric vectors
  int <- as.numeric(int[1,])
  pi <- as.numeric(pi)
  
  g <- sapply(1:length(pi), function(i) {
    log(pi[i]) + n_rep * log(int[i]) - int[i] * p_r - lfactorial(n_rep)
  })
  g <- exp(g)
  g <- rowSums(g)

  sup_g <- log(pi) + n_rep * log(int) - int * p_min - lfactorial(n_rep)
  sup_g <- exp(sup_g)
  sup_g <- sum(sup_g)
  
  return(g/sup_g)
}

# EM estimation
pois.HMM.EM <- function(n, g, long.data, log.probs, theta,  
                        gamma, pi, eta, mu,delta, lambda, run_off, D, maxiter =1000, tol=10^(-6) ,...)
{
  ncoef <- nrow(theta) 
  theta.next <- theta
  gamma.next <- gamma
  pi.next <- pi
  delta.next <- delta
  
  # reported claims
  n_r = rowSums(run_off)
  start = Sys.time()
  ## itteration
  for (iter in 1: maxiter){
    # obtain forward and backward probabilities and the likelihood
    fb <- pois.HMM.lalphabeta(t = n, log.prob.matrix =log.probs, 
                              g = g, gamma = gamma, pi=pi)
    lalpha <- fb$la 
    lbeta <- fb$lb 
    llik <- fb$loglik
    
    # Posterior distribution for p_t
    lambda.temp <- lambda*long.data$exposure
    
    p_hat <- matrix(0, nrow = nrow(long.data), ncol = D+1)

    pi_t <- exp(lalpha+lbeta-llik)
    pi_t <- pi_t/colSums(pi_t)
    
    time <- long.data$time
    index <- which(long.data$time <= (n-D))
    for(t in index){
      p_hat[t,] <- (eta$alpha[t,]+run_off[t,])/sum(eta$alpha[t,]+run_off[t,])
    }
    
    index <- which(long.data$time == n)
    for(t in index){
      temp <- data.frame(x =integer())
      #p_min <- min(rdirichlet(100000, alpha = eta$alpha[t,]+run_off[t,])[,1])
      while(nrow(temp)<=1){
        u <- rdirichlet(1000, alpha = eta$alpha[t,]+run_off[t,])
        w <- runif(1000, min = 0, max = 1)
        p_r <- u[,1]
        f <- acc_rej(p_r = p_r,pi = pi_t[,n], int = lambda.temp[t,], n_rep = n_r[t])
        temp <- as.data.frame(u[which(w<f),])
      }
      p_hat[t,] <- colMeans(temp)
    }
    
    index <- which((long.data$time > (n-D)) & (long.data$time < n))
    for(t in index){
      temp <- data.frame(x =integer())
      #p_min <- min(rowSums(rdirichlet(100000, alpha = eta$alpha[t,]+run_off[t,])[,1:(n-time[t]+1)]))
      while(nrow(temp)<=1){
        u <- rdirichlet(1000, alpha = eta$alpha[t,]+run_off[t,])
        w <- runif(1000, min = 0, max = 1)
        p_r <- rowSums(u[,1:(n-time[t]+1)])
        f <- acc_rej(p_r = p_r, pi = pi_t[,time[t]], int = lambda.temp[t,], n_rep = n_r[t])
        temp <- as.data.frame(u[which(w<f),])
      }
      p_hat[t,] <- colMeans(temp)
    }

    rm(lambda.temp, temp, pi_t, u, w, f, index,p_r,t)
    
    # estimate expected value for reporting probability
    index= which(time > (n-D))
    cumsum.delay <- t(apply(p_hat, 1, cumsum))
    rep.prob <- rep(1, nrow(long.data))
    for(t in index){
      rep.prob[t] <- cumsum.delay[t, n - time[t] + 1]
    }
    long.data$rep.prob <- rep.prob 
    long.data$ofst <- long.data$exposure * long.data$rep.prob
    rm(time,index,rep.prob)
    
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
    
    ## next eta's
    delay.vec <- DR_data(p_hat)
    dir_reg <- DirichReg(delay.vec ~ x2, data = long.data)
    delta.next <- matrix(as.numeric(dir_reg$coefficients),ncol=(D+1))
    eta <- predict.DirichletRegModel(dir_reg, newdata = long.data, alpha = TRUE)
    mu.next <- unique(eta$mu) 
    
    # checking convergence 
    crit <- sum(abs(theta -theta.next)) + sum(abs(gamma -gamma.next)) +
      sum(abs(pi -pi.next)) + sum(abs(mu -mu.next))
    print(crit)
    print(mu.next)
    # next itteration
    theta <- theta.next
    gamma <- gamma.next
    pi <- pi.next
    delta <- delta.next
    mu <- mu.next
    
    #log_probs
    log.probs <- lprobmatrix(long.data, lambda, alpha = eta$alpha, N =n, D = D)
    
    if(crit <tol){
      end = Sys.time()
      np <- g*g+g*ncoef-1 
      BIC <- -2*llik+np*log(n)
      return(list(BIC = BIC, theta=theta, pi = pi, gamma=gamma, mu=mu, lambda=lambda,
                  run.time = (end-start)/iter,iter=iter))
    }
    
  }
  print(paste ("No convergence after ",maxiter ," iterations "))
  NA
}
###################################################################################################
truncated_data <- censored_data %>% dplyr::filter(time <= N-D) #consider periods where we know everything

# delta (regression coefficients for eta's) 
p  = 3 #no. of regression coefficients
delta.init <- matrix(0, nrow = p, ncol = D+1)
eta.init <- 0.1*colSums(truncated_data[,paste0("d",0:D)])/sum(colSums(truncated_data[,paste0("d",0:D)]))
delta.init[1,] <- log(eta.init)
eta.init <- do.call(rbind, replicate(nrow(censored_data), eta.init, simplify = F))
eta.init <- list(alpha = eta.init)
mu.init <- colSums(truncated_data[,paste0("d",0:D)])/sum(colSums(truncated_data[,paste0("d",0:D)]))
mu.init <- matrix(rep(mu.init,3),nrow=3,ncol=D+1, byrow = T)

# Create Run off triangle
run_off <- censored_data %>% group_by(ID,time) %>%
  summarise(across(d0:paste0("d",D), ~sum(.x, na.rm = T)))
run_off <- as.matrix(run_off[,-c(1,2)])


#Initialization
g=2 # number of states

pi.init <- rep(1,g)/g #Initial distribution

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
log.probs <- lprobmatrix(censored_data, lambda.init,alpha = eta.init$alpha , N = N, D = D)
run.model <- pois.HMM.EM(n = N,g=g,long.data = censored_data, 
                         log.probs = log.probs, theta = theta.init,
                         gamma = gamma.init, eta = eta.init, mu = mu.init, delta = delta.init,
                         pi = pi.init, lambda = lambda.init, run_off = run_off, 
                         D = D, tol = 0.01)

#Print output
print(run.model$pi)
print(run.model$gamma)
print(run.model$theta)
print(run.model$mu)
print(run.model$run.time)
print(run.model$iter)
