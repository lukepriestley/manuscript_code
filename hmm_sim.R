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
behaviour_path <- "/Volumes/LaCie/nature_neuro_submission/data/behaviour/fmri_experiment/fmri_behaviour.csv"
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
# Function to calculate run lengths within a vector
calculate_run_lengths <- function(vec) {
  runs <- rle(vec)
  run_lengths <- runs$lengths
  return(run_lengths)
}

############### simulate data ##################
hmm_path <- '/Volumes/Lacie/nature_neuro_submission/data/behaviour/glm_hmm_fits/fmri/'

# select monkey to simulate from
monkey <- 'Winky'

# select subset to simulate on
d_sim <- d %>% filter(monkey == 'Winky') %>% select(session, trial, mag, prob, ev, richness_bin, stochasticity_bin, response_mod)
sessions <- unique(d_sim$session)

# set N HMM states to simulate from + N simulations
n_states <- c(1,2)
n_simulations <- 10

# Create an empty data frame to store simulation results
sim_results <- data.frame(y = numeric(),
                          hmm_state = numeric(),
                          n_hmm_states = numeric(),
                          simulation_number = integer(),
                          session_id = character(),
                          stringsAsFactors = FALSE)

set.seed(1) # set seed for reproduction

for(n in n_states){
  # get hmm fit + parameters
  filename <- paste(hmm_path, monkey, '_', n, '_states_fmri.rds', sep = '')
  hmm <- readRDS(filename)
  
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
  for(s in sessions){
    for (simulation_number in 1:n_simulations) {
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
      sim_results <- rbind(sim_results, data.frame(y = tmp_d$y,
                                                   hmm_state = tmp_d$hmm_state,
                                                   n_hmm_states = n,
                                                   simulation_number = simulation_number,
                                                   session_id = s))
    }
  }
}

############## compare ACF ##############

# Initialize empty data frame to store ACF results
acf_results_df <- data.frame()
max_lag <- 10

for (n_states in unique(sim_results$n_hmm_states)) {
  
  for (s in unique(sim_results$session_id)) {
    
    y <- sim_results$y[sim_results$session == s & sim_results$n_hmm_states == n_states]
    
    # Calculate ACF with a lag of 10 for each series
    acf_y <- acf(y, lag.max = max_lag, plot = FALSE)$acf
    
    # Create a data frame for the current time series and ACF results
    acf_data <- data.frame(
      session = rep(s, length(acf_y)),
      n_hmm_states = rep(n_states, length(acf_y)),
      Lag = 0:max_lag,
      ACF = acf_y
    )
    
    # Append the data frame to the results data frame
    acf_results_df <- dplyr::bind_rows(acf_results_df, acf_data)
  }
}

sim_acf <- acf_results_df %>%
  group_by(n_hmm_states, Lag) %>%
  summarize(
    Mean_ACF = mean(ACF),
    SEM_ACF = sd(ACF) / sqrt(n())
  ) %>%
  ungroup()

# Calculate ACF for observed data
observed_acf <- d_sim %>%
  group_by(session) %>%
  summarize(acf_values = list(acf(response_mod, lag.max = 10, plot = FALSE)$acf[1:11]))

# Convert ACF values to a matrix for easy calculation
acf_matrix <- do.call(rbind, observed_acf$acf_values)

# Calculate mean and standard error for each lag
mean_acf <- colMeans(acf_matrix, na.rm = TRUE)
sem_acf <- apply(acf_matrix, 2, function(col) sd(col, na.rm = TRUE) / sqrt(sum(!is.na(col))))

# Create a data frame for the mean ACF values and standard errors
mean_acf_df <- data.frame(
  Lag = 0:10,
  Mean_ACF = mean_acf,
  SEM_ACF = sem_acf
)

# Combine simulation and observed ACF results
combined_acf_results <- bind_rows(mean_acf_df %>%
                                    mutate(n_hmm_states = 0), sim_acf)

# graph observed vs simulated ACFs
fig2f <- ggplot(combined_acf_results %>% filter(Lag != '0' & n_hmm_states < 3), 
       aes(x = Lag, y = Mean_ACF, ymin = Mean_ACF-1.96*(SEM_ACF), ymax = Mean_ACF+1.96*(SEM_ACF), color = factor(n_hmm_states), fill = factor(n_hmm_states), group = factor(n_hmm_states))) +
  geom_line(color = 'black') + 
  geom_ribbon(alpha = 0.50, color = NA) + 
  #geom_pointrange(shape=21, size = 0.5, color = 'black', position = position_dodge(width = 0.15)) + 
  scale_y_continuous(name = "ACF", limits = c(-0.01, 0.42), breaks = seq(from = 0.00, to = 0.40, by = 0.10), labels = scales::number_format(accuracy = 0.1)) + 
  scale_x_continuous(name = "Lag", breaks = seq(0, 10, 2), limits = c(1, 10)) + 
  scale_color_manual(values = c('grey90', 'slategray1', 'darkblue'), name = NULL, labels = c('Observed data', 'Binomial GLM', 'GLM-HMM')) + 
  scale_fill_manual(values = c('grey90', 'slategray1', 'darkblue'), name = NULL, labels = c('Observed data', 'Binomial GLM', 'GLM-HMM')) + 
  theme_pubr() + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.position = c(0.6, 0.85),
    axis.text.y = element_text(size = 14),
    legend.background = element_blank(),
    aspect.ratio = 1)
fig2f
############## compare run lengths ##############

# Initialize an empty list to store results
sim_run_lengths <- list()

# Loop through unique combinations of session, model, and simulation
for (session_id in unique(sim_results$session_id)) {
  for (n_hmm_states in unique(sim_results$n_hmm_states)) {
    for (simulation_number in unique(sim_results$simulation_number)) {
      
      # Subset the data for the current combination
      subset_data <- sim_results[
        sim_results$session_id == session_id &
          sim_results$n_hmm_states == n_hmm_states &
          sim_results$simulation_number == simulation_number, 
      ]
      
      # Calculate run lengths using rle
      rle_result <- rle(subset_data$y)
      
      # Extract run lengths and values
      run_lengths <- rle_result$lengths
      values <- rle_result$values
      
      # Store the results in the list
      sim_run_lengths[[length(sim_run_lengths) + 1]] <- data.frame(
        session_id = session_id,
        n_hmm_states = n_hmm_states,
        simulation_number = simulation_number,
        run_length = run_lengths,
        value = values
      )
    }
  }
}

# Combine the list of results into a data frame
sim_run_lengths <- do.call(rbind, sim_run_lengths)

# convert to proportions
sim_proportions <- sim_run_lengths %>%
  group_by(n_hmm_states, run_length) %>%
  summarize(total_run_length = n()) %>%
  group_by(n_hmm_states) %>%
  mutate(total_model_run_length = sum(total_run_length)) %>%
  mutate(proportion = total_run_length / total_model_run_length) %>%
  select(n_hmm_states, run_length, proportion)

# get run-lengths in observed data
# Initialize an empty list to store results
observed_run_lengths <- list()

# Loop through unique sessions
for (s in unique(d_sim$session)) {
  # Subset data for the current session
  subset_data <- d_sim[d_sim$session == s, ]
  
  # Calculate run lengths using rle
  rle_result <- rle(subset_data$response_mod)
  
  # Extract run lengths
  run_lengths <- rle_result$lengths
  
  # Store the results in the list
  observed_run_lengths[[length(observed_run_lengths) + 1]] <- data.frame(
    session = s,
    run_length = run_lengths
  )
}

observed_run_lengths <- do.call(rbind, observed_run_lengths)

observed_proportions <-  observed_run_lengths %>% 
  group_by(run_length) %>%
  summarise(n = n()) %>%
  mutate(proportion = n / sum(n)) %>% 
  select(-n)

# combine simulated and observed proportions
overall_proportions <- observed_proportions %>% 
  mutate(n_hmm_states = 0) %>% 
  rbind(sim_proportions)

fig2g <- ggplot(overall_proportions %>% 
         filter(run_length >= 5) %>% 
         mutate(proportion = round(proportion, 3)), aes(x = run_length, y = proportion, group = factor(n_hmm_states), fill = factor(n_hmm_states), color = factor(n_hmm_states))) +
  geom_area(alpha = 0.25, position = 'identity') + 
  #geom_bar(stat = 'identity', position = 'dodge2', color = 'black') +
  scale_x_continuous(breaks = seq(5, 20, 5), limits = c(5, 20)) + 
  scale_y_continuous(breaks = seq(0, 0.05, 0.025)) + 
  labs(x = "Run Length", y = "Proportion", fill = "Model") +
  theme_pubr() +
  scale_color_manual(values = c('grey15', 'skyblue2', 'darkblue'), name = NULL, labels = c('Observed data', 'Binomial GLM', 'GLM-HMM')) + 
  scale_fill_manual(values = c('grey90', 'skyblue1', 'darkblue'), name = NULL, labels = c('Observed data', 'Binomial GLM', 'GLM-HMM')) + 
  theme(legend.position = "top") + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.position = c(0.58, 0.85),
    axis.text.y = element_text(size = 14),
    legend.background = element_blank(),
    aspect.ratio = 1
  )
fig2g
################## Richness and Prev. response effects ###################

# model-free glm over one-state
simulations = unique(sim_results$simulation_number)
sessions = unique(sim_results$session_id)
n_states = unique(sim_results$n_hmm_states)

sim_results$prev_response <- NA
sim_results$mag <- NA
sim_results$prob <- NA
sim_results$reward <- NA
sim_results$reward_received <- NA
sim_results$ave_rw <- NA

set.seed(1) #set seed for reprodiction

for(sim in simulations){
  sim_data <- sim_results %>% filter(simulation_number==sim)
  sim_data$mag <- d_sim$mag
  sim_data$prob <- d_sim$prob
  sim_data$reward <- ifelse(sim_data$y==0,0,
                            ifelse(sim_data$y==1, rbinom(n = length(sim_data$y[sim_data$y==1]), size = 1, prob = sim_data$prob), NA))
  sim_data$reward_received <- sim_data$reward * sim_data$mag
  for(sess in sessions){
    sess_data <- sim_data %>% filter(session_id==sess)
    for(t in 2:nrow(sess_data)){
      sess_data$prev_response[t] <- sess_data$y[t-1]
    }
    for(t in 6:nrow(sess_data)){
      sess_data$ave_rw[t] <- mean(sess_data$reward_received[(t-5):(t-1)])
    }
    sim_data[sim_data$session_id==sess,] <- sess_data
  }
  sim_results[sim_results$simulation_number==sim,] <- sim_data
}

sim_results$mag_z <- NA
sim_results$prob_z <- NA
sim_results$ave_rw_z <- NA
for(n in n_states){
  for(s in simulations){
    tmp_data <- sim_results %>% filter(simulation_number==s & n_hmm_states==n)
    tmp_data$mag_z <- scale(tmp_data$mag)
    tmp_data$prob_z <- scale(tmp_data$prob)
    tmp_data$ave_rw_z <- scale(tmp_data$ave_rw)
    sim_results[sim_results$simulation_number==s & sim_results$n_hmm_states==n,] <- tmp_data
    rm(tmp_data)
  }
}

sim_results$n_hmm_states <- ifelse(sim_results$n_hmm_states==1, '1-state', '2-states')

true_behav <- d %>% filter(monkey == 'Winky') %>% select(session, trial, ev, mag, prob, richness_bin, stochasticity_bin, ave_rw, prev_response)

fig2h <-  ggplot(sim_results %>%
                      drop_na(ave_rw_z) %>%
                      mutate(ave_rw_ind = cut(ave_rw_z, breaks = unique(quantile(ave_rw_z, probs = seq(0, 1, length.out = 11))), labels = FALSE, include.lowest = TRUE)) %>%
                      group_by(ave_rw_ind) %>%
                      mutate(ave_rw_graph = mean(ave_rw_z, na.rm = T)) %>%
                      ungroup() %>%
                      group_by(n_hmm_states, ave_rw_graph) %>% summarise(m = mean(y, na.rm = T)*100), 
                    aes(x = ave_rw_graph, y = m, fill = n_hmm_states, color = n_hmm_states)) + 
  stat_summary(fun.data = 'mean_se', geom = 'pointrange', size = 0.75, color = 'black', shape = 21) + 
  geom_smooth(method = 'lm', linewidth = 0.65, color = 'black') + 
  scale_y_continuous(name = 'Pursue-rate [%]', limits = c(0, 100), breaks = c(20, 60, 100)) + 
  scale_x_continuous(name = 'Env. richness [Z]', breaks = seq(-2, 3,1), labels = scales::number_format(accuracy = 0.1)) + 
  scale_fill_manual(name = '', values = c('skyblue1', 'darkblue'), labels = c('Binomial GLM', 'GLM-HMM [2-states]')) + 
  theme_pubr() + 
  theme(
    text = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 11),
    legend.position = c(0.45, 0.95),
    aspect.ratio = 1,
    legend.background = element_blank()
  )   + 
  coord_cartesian(ylim=c(50, 100)) + 
  border('black')
fig2h

fig2i <- ggplot(sim_results %>%
                          mutate(n_hmm_states = ifelse(n_hmm_states=="1-state", 'Binomial-GLM', 'GLM-HMM')) %>%
                          drop_na(prev_response) %>%
                          group_by(n_hmm_states, simulation_number, session_id, prev_response) %>% 
                          summarise(m = mean(y, na.rm = T)*100),
                        aes(x = factor(prev_response), y = m, fill = n_hmm_states)) +  
  stat_summary(fun = 'mean', geom = 'bar', color = 'black', width = 0.60, alpha = 0.75) + 
  geom_point(shape = 21, color = 'black', size = 2, alpha = 0.35) + 
  scale_y_continuous(name = 'Pursue-rate[%]', breaks = c(20, 60, 100)) + 
  scale_x_discrete(name = 'Behavioural hist.', labels = c('Rej.', 'Purs.')) + 
  scale_fill_manual(name = '', values = c('skyblue1', 'darkblue'), labels = c('Binomial GLM', 'GLM-HMM')) + 
  theme_pubr() + 
  theme(
    text = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 12),
    aspect.ratio = 1.5,
    legend.background = element_blank(), 
    legend.position = 'none') + 
  border('black') + 
  facet_wrap(~n_hmm_states) + 
  coord_cartesian(ylim = c(20, 100))
fig2i
