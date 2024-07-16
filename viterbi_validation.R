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

################### Read data ################
simulation_path <- "/Volumes/LaCie/neuron_submission/data/behaviour/glm_hmm_param_recovery/" # where simulated behaviour + model fits to simulated behaviour are located
sim_results_path <- paste(simulation_path, 'param_recovery_sim.csv', sep = '')
sim_results <- read.csv(sim_results_path)

################# Specify functions ###################
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

########################## viterbi decoding ############################
n_sim <- max(sim_results$simulation_number) # how many simulations were performed?
n_sessions <- max(sim_results$session_id) # how many sessions per simulation?
n_states <- 2 # how many HMM-states are we performing Viterbi decoding for? 

# pre-allocate array storing true-vs-decoded state on each trial in each simulation
sim_viterbi <- data.frame(
  simulation = integer(),
  session = integer(), 
  true_state = integer(),
  decoded_state = integer()
)

simulation_dir <- "/Volumes/LaCie/neuron_submission/data/behaviour/glm_hmm_param_recovery/"

for(sim in 1:n_sim){
  print(paste('simulation ', sim, sep = ''))
  hmm_path <- paste(simulation_dir, '2_gen_', sim, '_sim_', n_states,'_states.rds', sep = '')
  hmm <- readRDS(hmm_path)
  A <- matrix(nrow = n_states, ncol = n_states)
  for(i in 1:dim(A)[1]){
    for(j in 1:dim(A)[2]){
      A[i,j] <- get_posterior_mean(hmm, pars = paste('A[',i,',',j,']', sep = ''))[5]
    }
  }
  pi <- numeric(length = n_states)
  for(i in 1:length(pi)){
    pi[i] <- get_posterior_mean(hmm, pars = paste('pi1[',i,']', sep = ''))[5]
  }
  alpha <- numeric(length = n_states)
  for(i in 1:length(alpha)){
    alpha[i] <- get_posterior_mean(hmm, pars = paste('alpha[',i,']', sep = ''))[5]
  }
  beta_ev <- get_posterior_mean(hmm, pars = 'beta_ev')[5]
  beta_richness <- get_posterior_mean(hmm, pars = 'beta_rich')[5]
  beta_stochasticity <- get_posterior_mean(hmm, pars = 'beta_vol')[5]
  for(session in 1:n_sessions){
    true_state <- sim_results$hmm_state[sim_results$gen_hmm_states==n_states & sim_results$simulation_number==sim & sim_results$session_id==session]
    decoded_state <- viterbi_algorithm(y = sim_results$y[sim_results$gen_hmm_states==n_states & sim_results$simulation_number==sim & sim_results$session_id==session],
                                       ev = sim_results$ev[sim_results$gen_hmm_states==n_states & sim_results$simulation_number==sim & sim_results$session_id==session],
                                       richness = sim_results$richness[sim_results$gen_hmm_states==n_states & sim_results$simulation_number==sim & sim_results$session_id==session],
                                       stochasticity = sim_results$stochasticity[sim_results$gen_hmm_states==n_states & sim_results$simulation_number==sim & sim_results$session_id==session],
                                       K = n_states,
                                       A = A, 
                                       alpha = alpha, 
                                       beta_ev = beta_ev, 
                                       beta_richness = beta_richness, 
                                       beta_stochasticity = beta_stochasticity)
    sim_viterbi <- rbind(sim_viterbi, data.frame(simulation = rep(sim, length(true_state)), 
                                                 session = rep(session, length(true_state)), 
                                                 true_state = true_state, 
                                                 decoded_state = decoded_state))
  }
}

########################## FIGURE 3 ############################

# Graph mean decoding accuracy
decoding_accuracy <- sim_viterbi %>%
  group_by(simulation, session) %>%
  summarise(
    mean_decoding_accuracy = mean(true_state==decoded_state, na.rm = TRUE)
  )

se_decoding_accuracy = sd(decoding_accuracy$mean_decoding_accuracy, na.rm = TRUE) / sqrt(sum(!is.na(decoding_accuracy$mean_decoding_accuracy)))

fig3a <- ggplot(decoding_accuracy, aes(x = mean_decoding_accuracy*100)) + 
  geom_histogram(aes(y = after_stat(count / sum(count))*100), color = 'black', bins = 20, fill = 'skyblue1', alpha = 0.15) + 
  geom_vline(xintercept = 91, color = 'blue', linetype = 2) +
  scale_y_continuous(name = '% of sessions', breaks = seq(0, 15, 5), limits = c(0, 15)) + 
  scale_x_continuous(name = 'Decoding accuracy [%]', limits = c(75, 100), breaks = seq(80, 100, 10)) + 
  theme_pubr() + 
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )
fig3a

# Graph Viterbi-decoded vs true HMM-state time series for example session
example_data <- sim_viterbi %>% filter(simulation==1, session==4)
example_data$time <- 1:nrow(example_data)
example_data <- example_data %>% pivot_longer(cols = ends_with('_state'), values_to = 'state_value', names_to = 'source')
example_data$state_value[example_data$source=='decoded_state'] <- example_data$state_value[example_data$source=='decoded_state'] + 0.005

fig3b <- ggplot(example_data, aes(x = time, y = state_value-1, color = source, linetype = source)) + 
  geom_line(position = position_dodge(width = 1)) + 
  scale_y_continuous(name = 'HMM-state', limits = c(0, 1.1), breaks = c(0, 1.1), labels = c('Low', 'High')) + 
  scale_x_continuous(name = 'Trial [N]', limits = c(0, 191), breaks = seq(0, 200, 40)) + 
  theme_pubr() + 
  scale_color_manual(name = 'HMM-state', values = c('skyblue1', 'darkblue'), labels = c('Decoded', 'True')) + 
  scale_linetype_manual(name = 'HMM-state', values = c(1, 2), labels = c('Decoded', 'True')) + 
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    axis.text.y = element_text(size = 14, angle = 90),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.position = 'none',
    aspect.ratio = 1,
  )
fig3b

# Identify state transitions in true and decoded HMM-state time series
sim_viterbi$true_state_difference <- NA
sim_viterbi$decoded_state_difference <- NA

for (simulation in 1:n_sim) { # loop through simulations
  for (session in 1:n_sessions) { # loop through sessions
    subset_data <- sim_viterbi[sim_viterbi$simulation == simulation & sim_viterbi$session == session, ]
    if (nrow(subset_data) > 1) {
      subset_data$true_state_difference <- c(NA, diff(subset_data$true_state))
      subset_data$decoded_state_difference <- c(NA, diff(subset_data$decoded_state))
      sim_viterbi[sim_viterbi$simulation == simulation & sim_viterbi$session == session, ] <- subset_data
    }
    rm(subset_data)
  }
}

# create binarised indicator of state-transition
sim_viterbi$true_state_transition <- ifelse(sim_viterbi$true_state_difference!=0, 1, 0)
sim_viterbi$decoded_state_transition <- ifelse(sim_viterbi$decoded_state_difference!=0, 1, 0)

# calculate distance between state transitions 
d_transition <- data.frame(simulation = integer(),
                           session = integer(),
                           transition_type = integer(),
                           decoded_index = integer(),
                           true_index = integer(),
                           distance = integer())

n_simulations <- max(sim_viterbi$simulation)
n_sessions <- max(sim_viterbi$session)

for(sim in 1:n_simulations){
  for(sess in 1:n_sessions){
    subset_data <- sim_viterbi %>% filter(simulation==sim & session==sess)
    
    decoded_indices <- which(subset_data$decoded_state_difference != 0)
    
    for(i in decoded_indices){
      # Find the indices with same state difference
      decoded_difference = subset_data$decoded_state_difference[i]
      true_indices <- which(subset_data$true_state_difference == decoded_difference)
      
      # Calculate the absolute differences in row indices
      row_diff <- abs(i - true_indices)
      
      # Find the minimum absolute difference
      min_diff <- min(row_diff)
      
      # Find the index of the nearest row with 'true_state_transition' == 1
      nearest_index <- true_indices[which(row_diff == min_diff)[1]]
      
      d_transition <- rbind(d_transition, 
                            data.frame(simulation = sim,
                                       session = sess,
                                       transition_type = decoded_difference,
                                       decoded_index = i,
                                       true_index = nearest_index,
                                       distance = min_diff))
    }
  }
}

d_transition$signed_difference <- d_transition$decoded_index - d_transition$true_index

d_transition %>% group_by(transition_type) %>% summarise(m = mean(signed_difference),
                                                         sd = sd(signed_difference))
fig3c <- ggplot(d_transition, aes(x = signed_difference)) + 
  geom_histogram(aes(y = after_stat(count / sum(count))*100), color = 'black', bins = 20, fill = 'skyblue1', alpha = 0.25) + 
  geom_vline(xintercept = -0.026, color = 'darkblue', linetype = 2) +
  scale_y_continuous(name = '% of transitions') + 
  scale_x_continuous(limits = c(-10, 10), breaks = seq(-10, 10, 5)) + 
  labs(
    y = 'Count', 
    x = "[Decoded - True]"
  ) + 
  theme_pubr() + 
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    aspect.ratio = 1.2
  )
fig3c
