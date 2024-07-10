rm(list = ls())

# load packages
library(rstan)
library(dplyr)
library(ggplot2)
library(ggpubr)

############### prepare data ############
behaviour_path <- "/Volumes/LaCie/nature_neuro_submission/data/behaviour/fmri_experiment/fmri_behaviour.csv"
d <- read.csv(behaviour_path)

############### define functions ##################

# Gernerate folds
sample_into_folds <- function(sessions, size) {
  # Initialize list to store the sampled subsets
  folds <- list()
  
  # Perform sampling until there are no elements left
  while (length(sessions) >= size) {
    set.seed(123)  # Set a seed for reproducibility
    
    # Determine the sample size (size or the remainder if fewer than size elements are left)
    sample_size <- ifelse(length(sessions) > size, size, length(sessions))
    
    # Randomly sample elements without replacement
    sampled_subset <- sample(sessions, size = sample_size, replace = FALSE)
    
    # Store the sampled subset
    folds <- c(folds, list(sampled_subset))
    
    # Remove the sampled elements from the vector
    sessions <- setdiff(sessions, sampled_subset)
  }
  
  # If there's a remainder, add it to the last subset
  if (length(sessions) > 0) {
    last_subset <- unlist(folds[length(folds)])
    last_subset <- c(last_subset, sessions)
    folds[length(folds)] <- list(last_subset)
  }
  
  return(folds)
}

log_sum_exp <- function(x) {
  max_x <- max(x)
  sum_exp <- sum(exp(x - max_x))
  log_sum_exp <- max_x + log(sum_exp)
  return(log_sum_exp)
}

fwd_algorithm <- function(y, ev, richness, stochasticity, n_state, pi, A, alpha, beta_ev, beta_rich, beta_stoch){
  
  session_length = length(y)
  logalpha <- matrix(0, nrow = session_length, ncol = n_state)
  # Forward algorithm implementation
  for (j in 1:n_state) {
    logalpha[1, j] <- log(pi[j]) +
      dbinom(y[1], size = 1,
             prob = plogis(alpha[j] + beta_ev * ev[1] +
                             beta_rich * richness[1] +
                             beta_stoch * stochasticity[1]),
             log = TRUE)
  }
  
  for (t in 2:session_length) {
    for (j in 1:n_state) {
      accumulator <- numeric(n_state)
      for (i in 1:n_state) {
        accumulator[i] <- logalpha[t - 1, i] +
          log(A[i, j]) +
          dbinom(y[t], size = 1,
                 prob = plogis(alpha[j] + beta_ev * ev[t] +
                                 beta_rich * richness[t] +
                                 beta_stoch * stochasticity[t]),
                 log = TRUE)
      }
      logalpha[t, j] <- log_sum_exp(accumulator)
    }
  }
  log_lik = log_sum_exp(logalpha[session_length,])
  return(log_lik)
}

############### specify glm-hmm ################

glm_hmm <- '

data {
    int<lower=1> n_states; // number of states
    int<lower=1> S; // count of sessions
    int<lower=1> T; // count of trials
    int<lower=1> session_index[S, 2]; // session indices
    int<lower=1> session_length[S]; // length of each session
    int y[T]; // observations
    real ev[T]; // expected value
    real richness[T]; // environment richness
    real stochasticity[T]; //environment stochasticity
}

parameters {

    // state-transition model
    simplex[n_states] pi1; // initial state distribution
    simplex[n_states] A[n_states]; // transition probability matrix

    // state-dependent observation model
    ordered[n_states] alpha; // state-dependent intercept
    real beta_ev; // state-dependent parameters
    real beta_rich; 
    real beta_stoch; 
    
}

model {

    // priors
    for(j in 1:n_states) {
        alpha[j] ~ normal(0, 3);
    }
    beta_ev ~ normal(0, 2);
    beta_rich ~ normal(0, 2); 
    beta_stoch ~ normal(0, 2); 
    
    //forward

    for (s in 1:S){

    vector[n_states] logalpha[session_length[s]];
    vector[n_states] observation_likelihood[session_length[s]];

    { // Observation likelihood
    for(t in 1:session_length[s]){
        for(j in 1:n_states) {
            observation_likelihood[t, j] = bernoulli_logit_lpmf(y[session_index[s,1]:session_index[s,2]][t] | alpha[j] + beta_ev * ev[session_index[s,1]:session_index[s,2]][t] + beta_rich * richness[session_index[s,1]:session_index[s,2]][t] + beta_stoch * stochasticity[session_index[s,1]:session_index[s,2]][t]);
                }
            }
     }

    { // forward algorithm implementation log p(z_t = j|x_{1:t})
    
    real accumulator[n_states];
    
    // special case for first obs
    for(j in 1:n_states)
        logalpha[1, j] = log(pi1[j]) + observation_likelihood[1, j];

    for (t in 2:session_length[s]){
        for (j in 1:n_states){
            for (i in 1:n_states){
                accumulator[i] = logalpha[t-1, i] + log(A[i, j]) + observation_likelihood[t, j];
                }
                logalpha[t, j] = log_sum_exp(accumulator);
            }
         }
      } 
    target += log_sum_exp(logalpha[session_length[s]]); 
  }
}
'

############### perform 5-fold cross validation ###################
write_dir <- "/Volumes/LaCie/nature_neuro_submission/data/behaviour/glm_hmm_crossval"

K <- 5 # how many folds?
n_states <- 5 # how many GLM-HMM states to test?
m <- 'Vampire' # which animals to perform cross-validation on

# divide data into folds
d_monkey <- d %>% filter(monkey == m)
sessions = as.vector(unique(d_monkey$session))
fold_size <- 3
folds <- sample_into_folds(sessions, fold_size)

# get sessions indices
session_idx = matrix(nrow = length(sessions), ncol = 2)
session_length <- vector()
for (s in 1:length(sessions)) {
  session_idx[s, 1] <- min(which(d_monkey$session == sessions[s]))
  session_idx[s, 2] <- max(which(d_monkey$session == sessions[s]))
  session_length[s] <- nrow(d_monkey[d_monkey$session ==sessions[s],])
}

for(s in 1:n_states){
  for(k in 1:length(folds)){
    
    stan_data <- list(
      session_index = as.array(session_idx[-folds[[k]],]),
      session_length = as.array(session_length[-folds[[k]]]),
      n_states = s,
      S = length(sessions[-folds[[k]]]),
      T = nrow(d_monkey), 
      y = d_monkey$response_mod, 
      ev = d_monkey$ev, 
      richness = d_monkey$richness_bin, 
      volatility = d_monkey$stochasticity_bin
    )
    
    print(paste('fitting', m, 'fold', k))
    
    hmm_fit <-
      stan(
        model_code = glm_hmm,
        data = stan_data,
        iter = 1e4,
        cores = 4,
        chains = 4
      )
    
    filename = paste(m, '_', s, '_states_fold_', k, '.rds', sep = '')
    hmm_fit@stanmodel@dso <- new("cxxdso")
    saveRDS(hmm_fit, file = filename)
    rm(filename, hmm_fit, stan_data)
  }
}



############### run fwd algorithm on held-out sessions ##################

n_folds <- length(folds) # how many folds are there? 

# pre-allocate array of log-likelihoods for left out sessions
log_likelihoods <- vector("list", length = n_states)
for (model in 1:n_states) {
  log_likelihoods[[model]] <- vector("list", length = n_folds)
  for (fold in 1:n_folds) {
    n_sessions_in_fold <- length(folds[[fold]])
    log_likelihoods[[model]][[fold]] <- vector("list", length = n_sessions_in_fold)
    for (session in 1:n_sessions_in_fold) {
      log_likelihoods[[model]][[fold]][[session]] <- numeric()
    }
  }
}

for(n in 1:n_states){
  for(k in 1:n_folds){
    # get hmm fit + parameters
    hmm_path <- paste(write_dir,'/', m, '_', n, 'states_fold_', k, '.rds', sep = '')
    print(paste(m, n, 'states fold', k, sep = ' '))
    hmm <- readRDS(paste(hmm_path))
    A <- matrix(nrow = n, ncol = n)
    for(i in 1:dim(A)[1]){
      for(j in 1:dim(A)[2]){
        A[i,j] <- get_posterior_mean(hmm, pars = paste('A[',i,',',j,']', sep = ''))[5]
      }
    }
    pi <- numeric(length = n)
    for(i in 1:length(pi)){
      pi[i] <- get_posterior_mean(hmm, pars = paste('pi1[',i,']', sep = ''))[5]
    }
    alpha <- numeric(length = n)
    for(i in 1:length(alpha)){
      alpha[i] <- get_posterior_mean(hmm, pars = paste('alpha[',i,']', sep = ''))[5]
    }
    beta_ev <- get_posterior_mean(hmm, pars = 'beta_ev')[5]
    beta_richness <- get_posterior_mean(hmm, pars = 'beta_rich')[5]
    beta_stochasticity <- get_posterior_mean(hmm, pars = 'beta_vol')[5]
    n_sessions_in_fold <- length(folds[[k]])
    for(f in 1:n_sessions_in_fold){
      log_likelihoods[[n]][[k]][[f]] <- fwd_algorithm(y = d_monkey$response_mod[d_monkey$session == folds[[k]][f]],
                                                      ev = d_monkey$ev[d_monkey$session == folds[[k]][f]],
                                                      richness = d_monkey$richness_bin[d_monkey$session == folds[[k]][f]],
                                                      stochasticity = d_monkey$stochasticity_bin[d_monkey$session == folds[[k]][f]],
                                                      n_state = n,
                                                      pi = pi,
                                                      A = A,
                                                      alpha = alpha,
                                                      beta_ev = beta_ev,
                                                      beta_rich = beta_richness,
                                                      beta_stoch = beta_stochasticity
      )
    }
    rm(hmm, A, pi, alpha, beta_ev, beta_richness, beta_stochasticity)
  }
}


cross_val_log_lik <- data.frame(
  model = character(), 
  session = integer(), 
  fold = integer(), 
  log_likelihood = numeric()
)

for (model in 1:n_states) {
  for (fold in 1:n_folds) {
    n_sessions_in_fold <- length(folds[[fold]])
    for (session in 1:n_sessions_in_fold) {
      log_lik <- log_likelihoods[[model]][[fold]][[session]]
      cross_val_log_lik <- rbind(cross_val_log_lik, data.frame(model, session, fold, log_lik))
    }
  }
}

cross_val_log_lik$session_id <- NA
for(fold in unique(cross_val_log_lik$fold)){
  for(session in unique(cross_val_log_lik$session[cross_val_log_lik$fold==fold])){
    cross_val_log_lik$session_id[cross_val_log_lik$fold==fold & cross_val_log_lik$session==session] <- folds[[fold]][session]
  }
}

ggplot(data = cross_val_log_lik, aes(x = as.factor(model), y = log_lik)) + 
  geom_point(aes(x = as.factor(model), y = log_lik, group = session_id), color = 'pink', alpha = 0.75) + 
  geom_line(aes(x = as.factor(model), y = log_lik, group = session_id), color = 'pink', alpha = 0.75) + 
  stat_summary(aes(group = 1), fun = mean, geom = 'line', color = 'black') + 
  stat_summary(aes(group = 1), fun = mean, geom = 'point',shape = 22, color = 'black', fill = 'dark red', size = 3) + 
  scale_x_discrete(name = 'HMM-states [N]') + 
  scale_color_gradient(name = 'Session') + 
  theme_pubr() + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    strip.text = element_text(size = 14),
    strip.background = element_rect(colour='black', fill='grey99'),
    aspect.ratio = 1
  )
