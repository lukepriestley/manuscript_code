# Code for fiting a GLM-HMM for behavioural data recorded during fMRI sessions 

rm(list = ls())

# load packages
library(lme4)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggdist)
library(gghalves)
library(rstan)

############### prepare data ############

behaviour_path <- "/Volumes/LaCie/nature_neuro_submission/data/behaviour/fmri_experiment/fmri_behaviour.csv"
d <- read.csv(behaviour_path)

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

############### fit hmm-glm ###################
n_states <- 3:5 # how many GLM-HMM states to test?
write_dir <- "/Volumes/LaCie/nature_neuro_submission/data/behaviour/fmri_experiment/" # where to write models to?
setwd(write_dir)
monkeys <- unique(d$monkey)

for(m in unique(monkeys)){ # loop through monkeys
  
  d_monkey <- d %>% filter(monkey == m)
  sessions = as.vector(unique(d_monkey$session))
  # get sessions indices
  session_idx = matrix(nrow = length(sessions), ncol = 2)
  session_length <- vector()
  for (s in 1:length(sessions)) {
    session_idx[s, 1] <- min(which(d_monkey$session == sessions[s]))
    session_idx[s, 2] <- max(which(d_monkey$session == sessions[s]))
    session_length[s] <- nrow(d_monkey[d_monkey$session ==sessions[s],])
  }
  
  for(s in n_states){
    
    stan_data <- list(
      session_index = as.array(session_idx),
      session_length = as.array(session_length),
      n_states = s,
      S = length(sessions),
      T = nrow(d_monkey), 
      y = d_monkey$response_mod, 
      ev = d_monkey$ev_z, 
      richness = d_monkey$richness_bin, 
      stochasticity = d_monkey$stochasticity_bin
    )
    
    print(paste('fitting', m, 'with', s, 'states'))
    
    hmm_fit <-
      stan(
        model_code = glm_hmm,
        data = stan_data,
        iter = 1e4,
        cores = 4,
        chains = 4
      )
    
    filename = paste(m, '_', s, '_states_fmri.rds', sep = '')
    hmm_fit@stanmodel@dso <- new("cxxdso")
    saveRDS(hmm_fit, file = filename)
    rm(filename, hmm_fit, stan_data)
  }
}
