rm(list = ls())

# load packages
library(rstan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggdist)
library(gghalves)
library(scales)
library(latex2exp)

create_sim <- F # should new simulations be run? N.B option to use existing sims
fit_hmm <- F # should new GLM-HMMs be fitted to simulations? N.B. option to use existing fits

################### Read data ################
behaviour_path <- "/Volumes/LaCie/neuron_submission/data/behaviour/fmri_experiment/fmri_behaviour.csv"
d <- read.csv(behaviour_path)

################ specify functions ############

glm_hmm_pred <- function(ev, richness_bin, stochasticity_bin, n_states, pi, A, alpha, beta_ev, beta_rich, beta_stoch){
  
  n_trials <- length(ev) 
  y <- integer(n_trials)
  hmm_state <- integer(n_trials)
  
  # simulate evolution of hmm states
  hmm_state[1] <- sample(1:n_states, size = 1, prob = pi)
  for(t in 2:n_trials){
    hmm_state[t] <- sample(1:n_states, size =1, prob = A[hmm_state[t - 1], ])
  }
  
  # simulate state-dependent observations
  for (t in 1:n_trials) {
    state <- hmm_state[t]
    logit_prob <- alpha[state] + beta_ev * ev[t] + beta_rich * richness_bin[t] + beta_stoch * stochasticity_bin[t]
    prob <- plogis(logit_prob)
    y[t] <- rbinom(1, size = 1, prob = prob)
  }
  
  return(data.frame(y, hmm_state))
  
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

viterbi_algorithm <- function(y, ev, richness, stochasticity, K, A, alpha, beta_ev, beta_richness, beta_stochasticity) {
  
  zstar <- vector("integer", length(y))  # Store most likely state sequence
  session_length <- length(y)
  bpointer <- matrix(0, nrow = session_length, ncol = K)  # Backpointer matrix
  delta <- matrix(-Inf, nrow = session_length, ncol = K)   # Delta matrix
  
  for (j in 1:K) {
    delta[1, j] <- dbinom(y[1], size = 1,
                          prob = plogis(alpha[j] + beta_ev * ev[1] +
                                          beta_richness * richness[1] +
                                          beta_stochasticity * stochasticity[1]),
                          log = TRUE)
  }
  
  for (t in 2:length(y)) {
    for (j in 1:K) {
      delta[t, j] <- -Inf
      for (i in 1:K) {
        logp <- delta[t - 1, i] + log(A[i, j]) +
          dbinom(y[t], size = 1,
                 prob = plogis(alpha[j] + beta_ev * ev[t] +
                                 beta_richness * richness[t] +
                                 beta_stochasticity * stochasticity[t]),
                 log = TRUE)
        if (logp > delta[t, j]) {
          bpointer[t, j] = i
          delta[t, j] = logp
        }
      }
    }
  }
  
  logp_zstar <- max(delta[session_length, ])
  
  for(j in 1:K){
    if(delta[session_length, j] == logp_zstar){zstar[session_length] = j}
  }
  
  for (t in (1:(session_length - 1))) {
    zstar[session_length-t] <- bpointer[session_length - t + 1, zstar[session_length - t + 1]]
  }
  return(zstar)
}

################ run simulations ##################

if(create_sim){
  # select subset to simulate on
  sim_monkey <- 'Ultra'
  d_sim <- d %>% filter(monkey == sim_monkey) %>% select(session, trial, ev, mag, prob, richness_bin, stochasticity_bin)
  
  # select monkey to fit with 
  fit_monkey <- 'Winky'
  
  # Create an empty data frame to store simulation results
  sim_results <- data.frame(ev = numeric(),
                            richness = integer(),
                            volatility = integer(),
                            y = numeric(),
                            hmm_state = numeric(),
                            gen_hmm_states = numeric(),
                            session_id = numeric(),
                            simulation_number = numeric(),
                            stringsAsFactors = FALSE)
  
  n_states <- 2 # max number of HMM-states to simulate with
  n_sim <- 10 # max number of simulations
  hmm_path <- "/Volumes/LaCie/neuron_submission/data/behaviour/glm_hmm_fits/" # path to fitted HMM objects
  
  for (n in 1:n_states) {
    
    # get hmm fit + parameters
    model_path <- paste(hmm_path,fit_monkey, '_', n, '_states_fmri.rds', sep = '')
    hmm <- readRDS(model_path)
    
    A <- matrix(nrow = n, ncol = n)
    for (i in 1:dim(A)[1]) {
      for (j in 1:dim(A)[2]) {
        A[i, j] <- get_posterior_mean(hmm, pars = paste('A[', i, ',', j, ']', sep = ''))[5]
      }
    }
    
    pi <- numeric(length = n)
    for (i in 1:length(pi)) {
      pi[i] <- get_posterior_mean(hmm, pars = paste('pi1[', i, ']', sep = ''))[5]
    }
    
    alpha <- numeric(length = n)
    for (i in 1:length(alpha)) {
      alpha[i] <- get_posterior_mean(hmm, pars = paste('alpha[', i, ']', sep = ''))[5]
    }
    
    beta_ev <- get_posterior_mean(hmm, pars = 'beta_ev')[5]
    beta_richness <- get_posterior_mean(hmm, pars = 'beta_rich')[5]
    beta_stochasticity <- get_posterior_mean(hmm, pars = 'beta_stoch')[5]
    
    n_sessions <- length(unique(d_sim$session))
    
    for (s in 1:n_sessions) {
      for (simulation_number in 1:n_sim) {
        print(paste('session', s, 'simulation', simulation_number))
        tmp_d <- glm_hmm_pred(ev = d_sim$ev[d_sim$session == s],
                              richness_bin = d_sim$richness_bin[d_sim$session == s],
                              stochasticity_bin = d_sim$stochasticity_bin[d_sim$session == s],
                              n_states = n,
                              pi = pi,
                              A = A,
                              alpha = alpha,
                              beta_ev = beta_ev,
                              beta_rich = beta_richness,
                              beta_stoch = beta_stochasticity)
        
        # Append the simulation results to the data frame
        sim_results <- rbind(sim_results, data.frame(ev = d_sim$ev[d_sim$session == s], 
                                                     richness = d_sim$richness_bin[d_sim$session == s], 
                                                     stochasticity = d_sim$stochasticity_bin[d_sim$session == s],
                                                     y = tmp_d$y,
                                                     hmm_state = tmp_d$hmm_state,
                                                     gen_hmm_states = n,
                                                     session_id = s,
                                                     simulation_number = simulation_number))
      }
    }
  }
  
  sim_results %>% group_by(gen_hmm_states) %>% summarise(m = mean(y, na.rm = T))
  rm(hmm)
  write.csv(sim_results, file = 'param_recovery_sim.csv')
}

################ specify hmm-glm ################################
#################################################################

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

############### fit hmm-glm to simulated data ###################

if(create_sim){
  n_states <- 2 # how many GLM-HMM states to test?
  n_sim <- 10
  write_dir <- "/Volumes/LaCie/neuron_submission/data/behaviour/glm_hmm_param_recovery/"
  
  for(gen in 1:n_states){ # loop over generative numbers of HMM-states
    
    for(sim in 1:n_sim){ # loop over simulations
      d_fit <- sim_results %>% filter(gen_hmm_states == gen & simulation_number == sim)
      sessions = unique(d_fit$session_id)
      # get sessions indices
      session_idx = matrix(nrow = length(sessions), ncol = 2)
      session_length <- vector()
      for (sess in 1:length(sessions)) {
        session_idx[sess, 1] <- min(which(d_fit$session == sess))
        session_idx[sess, 2] <- max(which(d_fit$session == sess))
        session_length[sess] <- nrow(d_fit[d_fit$session == sess,])
      }
      
      for (fit in 1:n_states){
        
        stan_data <- list(
          session_index = as.array(session_idx),
          session_length = as.array(session_length),
          n_states = fit,
          S = length(sessions),
          T = nrow(d_fit), 
          y = d_fit$y, 
          ev = d_fit$ev, 
          richness = d_fit$richness, 
          stochasticity = d_fit$stochasticity
        )
        
        print(paste('fitting simulated data from simulation number', sim, ' with', gen, 'generative states with', fit, 'state glm-hmm'))
        
        hmm_fit <-
          stan(
            model_code = glm_hmm,
            data = stan_data,
            iter = 1e4,
            cores = 4,
            chains = 4
          )
        
        filename = paste(write_dir, gen, '_gen_', sim, '_sim_', fit, '_states.rds', sep = '')
        hmm_fit@stanmodel@dso <- new("cxxdso")
        saveRDS(hmm_fit, file = filename)
        rm(filename, hmm_fit, stan_data)
      }
    }
  }
}

#################### test log-likelihood #####################
# Can we recover the true number of HMM-states that generated the data using
# the log-likelihood of fitted models? 

simulation_path <- "/Volumes/LaCie/neuron_submission/data/behaviour/glm_hmm_param_recovery/" # where simulated behaviour + model fits to simulated behaviour are located
sim_results_path <- paste(simulation_path, 'param_recovery_sim.csv', sep = '')
sim_results <- read.csv(sim_results_path)

n_states <- 2
n_sim <- 10
n_sessions <- max(sim_results$session_id) # set to however many sessions in each simulation (i.e. max session_id)

# Pre-allocate a list that stores session-wise log-likelihood of the simulated data under
# 1-state and 2-state GLM-HMM

log_likelihoods <- vector("list", length = n_states)
for (gen in 1:n_states) { # loop through generative states (i.e. how many HMM-states generated the data?)
  for(fit in 1:n_states){ # loop through fitted states (i.e. how many HMM-states were fit to the data?)
    log_likelihoods[[gen]][[fit]] <- vector("list", length = n_sim)
    for(sim in 1:n_sim){ # loop through simulations
      log_likelihoods[[gen]][[fit]][[sim]] <- vector("list", length = n_sessions)
      for(session in 1:n_sessions){
        log_likelihoods[[gen]][[fit]][[sim]][[session]] <- numeric()
      }
    }
  }
}

for(gen in 1:n_states){
  for(fit in 1:n_states){
    for(sim in 1:n_sim){
      # get hmm fit + parameters
      print(paste(gen, ' generative states & ', fit, 'fitted states on simulation', sim, sep = ' '))
      hmm_path <- paste(simulation_path, gen, '_gen_', sim, '_sim_', fit, '_states.rds', sep = '')
      hmm <- readRDS(hmm_path)
      A <- matrix(nrow = fit, ncol = fit)
      for(i in 1:dim(A)[1]){
        for(j in 1:dim(A)[2]){
          A[i,j] <- get_posterior_mean(hmm, pars = paste('A[',i,',',j,']', sep = ''))[5]
        }
      }
      pi <- numeric(length = fit)
      for(i in 1:length(pi)){
        pi[i] <- get_posterior_mean(hmm, pars = paste('pi1[',i,']', sep = ''))[5]
      }
      alpha <- numeric(length = fit)
      for(i in 1:length(alpha)){
        alpha[i] <- get_posterior_mean(hmm, pars = paste('alpha[',i,']', sep = ''))[5]
      }
      beta_ev <- get_posterior_mean(hmm, pars = 'beta_ev')[5]
      beta_richness <- get_posterior_mean(hmm, pars = 'beta_rich')[5]
      beta_stochasticity <- get_posterior_mean(hmm, pars = 'beta_vol')[5]
      for(session in 1:n_sessions){
        log_likelihoods[[gen]][[fit]][[sim]][[session]] <- fwd_algorithm(y = sim_results$y[sim_results$gen_hmm_states==gen & sim_results$sim==sim & sim_results$session_id==session],
                                                                         ev = sim_results$ev[sim_results$gen_hmm_states==gen & sim_results$sim==sim & sim_results$session_id==session],
                                                                         richness = sim_results$richness[sim_results$gen_hmm_states==gen & sim_results$sim==sim & sim_results$session_id==session],
                                                                         stochasticity = sim_results$stochasticity[sim_results$gen_hmm_states==gen & sim_results$sim==sim & sim_results$session_id==session],
                                                                         n_state = fit,
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
}

sim_log_lik <- data.frame(
  gen = character(), 
  fit = character(),
  simulation = integer(),
  session = integer(), 
  log_likelihood = numeric()
)

for (gen in 1:n_states) {
  for (fit in 1:n_states) {
    for(sim in 1:n_sim){
      for (session in 1:n_sessions) {
        log_lik <- log_likelihoods[[gen]][[fit]][[sim]][[session]]
        sim_log_lik <- rbind(sim_log_lik, data.frame(gen, fit, sim, session, log_lik))
      }
    }
  }
}

############## FIGURE S3A #################

sim_log_lik$session_in_sim <- paste('sim_', sim_log_lik$sim, '_session_', sim_log_lik$session, sep = '')
sim_log_lik$gen <- factor(sim_log_lik$gen, labels = c('Gen. states = 1', 'Gen. states = 2'))
fig.s3a <- ggplot(data = sim_log_lik, aes(x = factor(fit), y = log_lik, group = session_in_sim)) + 
  geom_point(color = 'pink', alpha = 0.25) + 
  geom_line(color = 'pink', alpha = 0.25) + 
  stat_summary(aes(group = 1), fun = mean, geom = 'line', color = 'black') + 
  stat_summary(aes(group = 1), fun = mean, geom = 'point',shape = 22, color = 'black', fill = 'dark red', size = 3) + 
  scale_x_discrete(name = 'HMM-states [N]') + 
  scale_y_continuous(name = 'Log. lik.') + 
  theme_pubr() + 
  facet_grid(~gen) + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    plot.background = element_blank(),
    aspect.ratio = 1,
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines")
  )
fig.s3a 

##################### compare recovered vs true params (FIGURE S3B-C) #################
# Within 2-state GLM-HMMs, can we recover the parameters which generated the data?

recovered_params <- data.frame(
  sim_number = integer(),
  hmm_state = integer(),
  pi = numeric(),
  p_stay = numeric(),
  p_transition = numeric(),
  alpha = numeric(),
  beta_rich = numeric(),
  beta_stoch = integer(),
  stringsAsFactors = FALSE
)

d_alpha <- data.frame(simulation = factor(),
                      hmm_state = factor(),
                      alpha_value = numeric()
)

d_tm <- data.frame(
  simulation = factor(),
  transition = factor(),
  sample = numeric()
)

n_states <- 2

for(gen in n_states){
  for(fit in n_states){
    for (sim in 1:n_sim){
      
      hmm_path <- paste(simulation_path, gen, '_gen_', sim, '_sim_', fit, '_states.rds', sep = '')
      hmm <- readRDS(hmm_path)
      
      # get posterior over alpha values
      alpha <- rstan::extract(hmm, pars = 'alpha') %>% data.frame() %>% pivot_longer(cols = everything(), names_to = 'hmm_state', values_to = 'alpha_value')
      d_alpha <- rbind(d_alpha, 
                       data.frame(simulation = rep(sim, each = nrow(alpha)),
                                  hmm_state = alpha$hmm_state,
                                  alpha_value = alpha$alpha_value)
      )
      
      tm_samples <- rstan::extract(hmm, pars = 'A')
      sample_matrix <- matrix(unlist(tm_samples), ncol = 4, byrow = FALSE)
      sample_matrix <- reshape2::melt(sample_matrix) %>% select(-Var1)
      d_tm <- rbind(d_tm, 
                    data.frame(simulation = rep(sim, each = nrow(sample_matrix)),
                               transition = sample_matrix$Var2,
                               sample = sample_matrix$value)
      )
      
      # get specific param values
      A <- matrix(nrow = fit, ncol = fit)
      for(i in 1:dim(A)[1]){
        for(j in 1:dim(A)[2]){
          A[i,j] <- get_posterior_mean(hmm, pars = paste('A[',i,',',j,']', sep = ''))[5]
        }
      }
      pi <- numeric(length = fit)
      for(i in 1:length(pi)){
        pi[i] <- get_posterior_mean(hmm, pars = paste('pi1[',i,']', sep = ''))[5]
      }
      alpha <- numeric(length = fit)
      for(i in 1:length(alpha)){
        alpha[i] <- get_posterior_mean(hmm, pars = paste('alpha[',i,']', sep = ''))[5]
      }
      beta_ev <- get_posterior_mean(hmm, pars = 'beta_ev')[5]
      beta_richness <- get_posterior_mean(hmm, pars = 'beta_rich')[5]
      beta_stochasticity <- get_posterior_mean(hmm, pars = 'beta_vol')[5]
      
      for(state in 1:n_states){
        recovered_params <- rbind(recovered_params, data.frame(
          sim_number = sim,
          hmm_state = state,
          pi = pi[state],
          p_stay = A[state, state],
          p_transition = 1-A[state, state],
          alpha = alpha[state],
          beta_rich = beta_richness,
          beta_stoch = beta_stochasticity
        ))
      }
    }
  }
}

# Graph recovered vs true transition probabilites
d_tm <- d_tm %>%
  mutate(transition_facet = case_when(
    transition == 1 ~ 'Low-to-low',
    transition == 2 ~ 'High-to-low',
    transition == 3 ~ 'Low-to-high',
    transition == 4 ~ 'High-to-high',
    TRUE ~ NA_character_
  ))

d_tm %>% group_by(transition_facet) %>% summarise(m = mean(sample)) # check mean recovered transition probabilities

d_hline_sim <- data.frame(transition_facet = rep(c('Low-to-low','High-to-low','Low-to-high','High-to-high')),
                          transition = c(1, 2, 3, 4),
                          value = c(0.922, 0.0471, 0.0781, 0.953)) # mean values per above

d_hline_true <- data.frame(transition_facet = rep(c('Low-to-low','High-to-low','Low-to-high','High-to-high')),
                           transition = c(1, 2, 3, 4),
                           value = c(0.926, 0.0455, 0.0738, 0.9544)) # true values for 2-state GLM-HMM of fit monkey

fig.s3b <- ggplot(d_tm %>% filter(sample > 0.87|sample < 0.14), aes(y = sample, group = simulation, fill = transition)) + 
  stat_halfeye(adjust = 1, justification = 0, .width = 0, alpha = 0.20, outline_bars = T, point_alpha = 0.6, point_color = 'black', shape = 21) + 
  labs(
    y = TeX("p(Transition)")
  ) +
  scale_y_continuous(breaks = seq(0, 1, 0.05)) +
  theme_pubr() + 
  scale_fill_gradient(low = "pink", high = "darkred") +
  scale_color_gradient(low = "pink", high = "darkred") +
  border(color = 'black') + 
  geom_hline(data = d_hline_sim, aes(yintercept = value, color = transition), linetype = 2) + 
  geom_hline(data = d_hline_true, aes(yintercept = value), linetype = 4, color = 'black') + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.position = 'none',
    legend.background = element_blank(),
    aspect.ratio = 1,
  ) + 
  facet_wrap(~transition_facet, nrow = 2, scales = 'free_y')
fig.s3b


# graph recovered vs true state-specific bias parameters
fig.s3c <- ggplot(d_alpha, aes(y = alpha_value, x = hmm_state, fill = hmm_state, color = hmm_state, group = simulation)) + 
  stat_halfeye(adjust = 1, justification = 0, .width = 0, alpha = 0.25, outline_bars = T, point_alpha = 0.6, point_color = 'black', shape = 21) + 
  labs(
    x = 'HMM-state', 
    y = TeX("$\\alpha$-coeff. [a.u.]")
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), breaks = seq(-1.5, 1.5, 1.5)) +
  scale_fill_manual(name = 'State', values = c('pink', 'darkred'), labels = c('Low', 'High')) + 
  scale_color_manual(name = 'State', values = c('pink', 'darkred'), labels = c('Low', 'High')) + 
  geom_hline(yintercept = 1.50, color = 'darkred', linetype = 2) + 
  geom_hline(yintercept = -1.30, color = 'pink', linetype = 2) + 
  geom_hline(yintercept = 1.46, color = 'black', linetype = 4) + 
  geom_hline(yintercept = -1.25, color = 'black', linetype = 4) + 
  theme_pubr() + 
  border(color = 'black') + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.position = 'none',
    legend.background = element_blank(),
    aspect.ratio = 1.5,
  )
fig.s3c
