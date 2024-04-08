# Code for all analysis involving GLM2, GLM-HMM for behavioural data recorded 
# during fMRI sessions ("A Hidden Markov Model identifies motivation-states in 
# behaviour"),fig.2A-E and fig.3D-I, and Supplementary Figure S2B-D. 

rm(list = ls())

# load packages
library(lme4)
library(rstan)
library(car)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gghalves)
library(ggdist)
library(ggpubr)
library(scales)
library(latex2exp)

################### Read data ################
behaviour_path <- "/Volumes/LaCie/cell_submission/data/behaviour/fmri_experiment/fmri_behaviour.csv"
d <- read.csv(behaviour_path)

############### specify functions ############

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

calculate_AIC <- function(log_likelihood, num_parameters) {
  AIC_score <- -2 * log_likelihood + 2 * num_parameters
  return(AIC_score)
}

############### implement Viterbi decoding ############
glm_hmm_dir <- "/Volumes/LaCie/cell_submission/data/behaviour/fmri_experiment/glm_hmm/"
d$viterbi_state <- NA
n_states <- 2
monkeys <- unique(d$monkey)
for(m in monkeys){
  
  # get hmm-glm
  hmm <- readRDS(paste(glm_hmm_dir,m, '_', n_states, '_states_fmri.rds', sep = ''))
  
  A <- matrix(nrow = n_states, ncol = n_states)
  for (i in 1:dim(A)[1]) {
    for (j in 1:dim(A)[2]) {
      A[i, j] <- get_posterior_mean(hmm, pars = paste('A[', i, ',', j, ']', sep = ''))[5]
    }
  }
  
  pi <- numeric(length = n_states)
  for (i in 1:length(pi)) {
    pi[i] <- get_posterior_mean(hmm, pars = paste('pi1[', i, ']', sep = ''))[5]
  }
  
  alpha <- numeric(length = n_states)
  for (i in 1:length(alpha)) {
    alpha[i] <- get_posterior_mean(hmm, pars = paste('alpha[', i, ']', sep = ''))[5]
  }
  
  beta_ev <- get_posterior_mean(hmm, pars = 'beta_ev')[5]
  beta_richness <- get_posterior_mean(hmm, pars = 'beta_rich')[5]
  beta_stochasticity <- get_posterior_mean(hmm, pars = 'beta_stoch')[5]
  
  sessions = unique(d$session[d$monkey==m])
  
  for(s in sessions){
    
    y <- d$response_mod[d$monkey==m & d$session==s]
    ev_z <- d$ev_z[d$monkey==m & d$session==s]
    richness <- d$richness_bin[d$monkey==m & d$session==s]
    stochasticity <- d$stochasticity_bin[d$monkey==m & d$session==s]
    
    zstar <- viterbi_algorithm(
      y = y, 
      ev = ev_z,
      richness = richness,
      stochasticity = stochasticity, 
      K = n_states,
      A = A, 
      alpha = alpha, 
      beta_ev = beta_ev, 
      beta_richness = beta_richness,
      beta_stochasticity = beta_stochasticity)
    d$viterbi_state[d$session == s & d$monkey==m] <- zstar
    rm(y, ev_z, richness, stochasticity, zstar)
  }
  
  rm(hmm, A, pi, alpha, beta_ev, beta_richness, beta_stochasticity, sessions)
  
}

# identify state changes
d$state_difference <- NA # Initialize the state_difference column with NA
monkeys <- unique(d$monkey)
# Loop through unique combinations of 'monkey' and 'session_id'
for (m in monkeys) {
  sessions = unique(d$session[d$monkey==m])
  for (s in sessions) {
    subset_data <- d[d$monkey == m & d$session == s, ]
    
    if (nrow(subset_data) > 1) {
      subset_data$state_difference <- c(NA, diff(subset_data$viterbi_state))
      d[d$monkey == m & d$session == s, ] <- subset_data
    }
    rm(subset_data)
  }
  rm(sessions)
}

d$state_transition <- ifelse(d$state_difference!=0, 1, 0)
d$state_increase <- ifelse(d$state_difference==1, 1, 0)
d$state_decrease <- ifelse(d$state_difference==-1, 1, 0)

#write.csv(d, 'fmri_behaviour_with_glmhmm_beta.csv') # write data for fMRI analysis

################### GLMs over Viterbi states ######################
d$viterbi_state_bin <- ifelse(d$viterbi_state==1, 0, 
                              ifelse(d$viterbi_state==2, 1, NA))
for(m in unique(d$monkey)){
  d$RT_z[d$monkey==m] <- scale(log(d$RT[d$monkey==m]))
  d$viterbi_state_z[d$monkey==m] <- scale(d$viterbi_state[d$monkey==m])
}

# Define a function to reset count when viterbi_state changes
reset_count <- function(states) {
  n <- length(states)
  counts <- integer(n)
  counts[1] <- 1
  count <- 1
  
  for (i in 2:n) {
    if (states[i] == states[i - 1]) {
      count <- count + 1
    } else {
      count <- 1
    }
    counts[i] <- count
  }
  
  return(counts)
}

# Initialize columns in the data frame 'd'
d$time_in_state <- NA
d$state_id <- NA
d$transition_countdown <- NA
d$trial_per_block <- NA

# Get unique monkeys from the data frame
monkeys <- unique(d$monkey)

# Process data for each monkey
for(m in monkeys) {
  monkey_data <- d %>% filter(monkey == m)
  sessions <- unique(monkey_data$session)
  
  # Process data for each session
  for(s in sessions) {
    session_data <- monkey_data %>% filter(session == s)
    
    # Calculate time in state, state_id, and transition countdown
    session_data$time_in_state <- reset_count(session_data$viterbi_state)
    session_data$state_id <- cumsum(session_data$time_in_state == 1)
    session_data$transition_countdown <- -unlist(tapply(session_data$state_id, session_data$state_id, function(x) rev(seq_along(x) - 1)))
    
    # Process data for each block type
    blocks <- unique(session_data$blockType)
    for(b in blocks) {
      indices <- session_data$blockType == b
      session_data$trial_per_block[indices] <- 1:length(session_data$trial_per_block[indices])
    }
    
    # Update monkey_data with processed session_data
    monkey_data[monkey_data$session == s, ] <- session_data
  }
  
  # Update the main data frame with processed monkey_data
  d[d$monkey == m, ] <- monkey_data
}

# create indicator for time-windows around state increases and state decreases
d$increase_bin <- ifelse(d$viterbi_state==1&d$transition_countdown>=-5, -1,
                         ifelse(d$viterbi_state==2&d$time_in_state<=5, 1, NA))

d$decrease_bin <- ifelse(d$viterbi_state==2&d$transition_countdown>=-5, -1,
                         ifelse(d$viterbi_state==1&d$time_in_state<=5, 1, NA))

GLM2.1 <- glmer(response ~ mag_z * prob_z + viterbi_state_bin + trial_z + (viterbi_state_bin|monkey), data = d, family = 'binomial')
summary(GLM2.1)

# GLM2.2a: is behaviour different before-vs-after increases in motivation-state?
GLM2.2a <- glmer(response ~ mag_z*prob_z + increase_bin + trial_z + (increase_bin|monkey), data = d, family = 'binomial')
summary(GLM2.2a)
# GLM2.2a: is behaviour different before-vs-after decreases in motivation-state?
GLM2.2b <- glmer(response ~ mag_z*prob_z + decrease_bin + trial_z + (0 + decrease_bin|monkey), data = d, family = 'binomial')
summary(GLM2.2b)

GLM2.3 <- lmer(pupil_size_decision_z ~ mag_z * prob_z + viterbi_state + trial_z + (0 + viterbi_state|monkey), data = d)  
summary(GLM2.3)
Anova(GLM2.3)

GLM2.4 <- lmer(RT_z ~ mag_z * prob_z + viterbi_state + trial_z + (0 + viterbi_state|monkey), data = d)  
summary(GLM2.4)
Anova(GLM2.4)

GLM2.5 <- glmer(viterbi_state_bin ~ ev_z + ave_ev_z + trial_z + (ave_ev_z|monkey), data = d, family = 'binomial')
summary(GLM2.5)

################## FIGURES #####################
################################################

monkeys <- c('Ulrich', 'Ultra', 'Winky', 'Vampire')

# Initialize data frames
n_states <- 1:2
d_aic <- data.frame(monkey = factor(),
                    session = integer(),
                    n_hmm_states = integer(),
                    log_lik = numeric(),
                    aic = numeric()
)

d_tm <- data.frame(monkey = factor(),
                   s_t = factor(),
                   s_t_plus_1 = factor(),
                   prob = numeric()

)

d_alpha <- data.frame(monkey = factor(),
                      hmm_state = factor(),
                      alpha_value = numeric()
)

for (m in monkeys) {
  
  monkey_data <- d %>% filter(monkey == m)
  sessions <- unique(monkey_data$session)
  
  for(state in n_states){
    
    filename <- paste(glm_hmm_dir,m, '_', state, '_states_fmri.rds', sep = '')
    hmm <- readRDS(filename)
    
    # extract parameters
    A <- matrix(nrow = state, ncol = state)
    for (i in 1:dim(A)[1]) {
      for (j in 1:dim(A)[2]) {
        A[i, j] <- get_posterior_mean(hmm, pars = paste('A[', i, ',', j, ']', sep = ''))[5]
      }
    }
    
    pi <- numeric(length = state)
    for (i in 1:length(pi)) {
      pi[i] <- get_posterior_mean(hmm, pars = paste('pi1[', i, ']', sep = ''))[5]
    }
    
    alpha <- numeric(length = state)
    for (i in 1:length(alpha)) {
      alpha[i] <- get_posterior_mean(hmm, pars = paste('alpha[', i, ']', sep = ''))[5]
    }
    
    beta_ev <- get_posterior_mean(hmm, pars = 'beta_ev')[5]
    beta_richness <- get_posterior_mean(hmm, pars = 'beta_rich')[5]
    beta_stochasticity <- get_posterior_mean(hmm, pars = 'beta_stoch')[5]
    
    # calculate number of parameters
    t_matrix <- (state^2); init_state_distribution <- state; state_intercepts <- state; glm_params <- 3
    n_parameters <- t_matrix + init_state_distribution + state_intercepts + glm_params
    
    # calculate session-wise log-likelihood and AIC scores
    for(s in sessions){
      subset_monkey_data <- monkey_data %>% filter(session == s)
      log_lik <- fwd_algorithm(y = subset_monkey_data$response_mod,
                               ev = subset_monkey_data$ev_z,
                               richness = subset_monkey_data$richness_bin,
                               stochasticity = subset_monkey_data$stochasticity_bin, 
                               n_state = state,
                               A = A,
                               pi = pi,
                               alpha = alpha,
                               beta_ev = beta_ev,
                               beta_rich = beta_richness,
                               beta_stoch = beta_stochasticity
                               )
      aic <- calculate_AIC(log_lik, n_parameters)
      d_aic <- rbind(d_aic, data.frame(monkey = m,
                                       session = s,
                                       n_hmm_states = state,
                                       log_lik = log_lik,
                                       aic = aic))
      rm(subset_monkey_data, log_lik, aic)
    }
    
    if(state == 2){
      
      # extract transition matrix params for graphs
      dimnames(A) <- list(c("state_1", "state_2"), c("state_1", "state_2"))
      tmp_data <- data.frame(
        monkey = rep(m, times = state^2),
        prob = as.vector(A),
        s_t = rep(colnames(A), each = nrow(A)),
        s_t_plus_1 = rep(rownames(A), times = ncol(A))
      )
      # add transition matrix params to d_tm
      d_tm <- rbind(d_tm, tmp_data)
      
      # get posterior density of alpha coefficients
      alpha <- rstan::extract(hmm, pars = 'alpha') %>% data.frame() %>% pivot_longer(cols = everything(), names_to = 'hmm_state', values_to = 'alpha_value')
      d_alpha <- rbind(d_alpha, 
                       data.frame(monkey = rep(m, each = nrow(alpha)),
                       hmm_state = alpha$hmm_state,
                       alpha_value = alpha$alpha_value)
      )
      #rm(tmp_data, alpha)
    }
    #rm(filename, hmm, A, pi, alpha, beta_ev, beta_richness, beta_stochasticity, t_matrix, init_state_distribution, state_intercepts, glm_params, n_parameters)
  }
  #rm(sessions, monkey_data)
}


# figure 2a: example session
figure_data <- d %>% filter(monkey == 'Ulrich' & session == 15) %>% mutate(ev_norm = (ev)/3)

# Function to identify segments for a given variable
identify_segments <- function(data, var_name) {
  data %>%
    arrange(trial) %>%
    mutate(change = ifelse(lag(!!sym(var_name), default = first(!!sym(var_name))) != !!sym(var_name), 1, 0),
           segment = cumsum(change)) %>%
    group_by(segment) %>%
    summarize(start_trial = min(trial), end_trial = max(trial), state = first(!!sym(var_name))) %>%
    ungroup()
}

# Identifying segments for viterbi_state_bin
viterbi_segments <- identify_segments(figure_data, "viterbi_state_bin")

# Identifying segments for richness_bin
richness_segments <- identify_segments(figure_data, "richness_bin")

# Creating the plot
fig2a <- ggplot() +
  geom_line(data = figure_data, aes(x = trial, y = ev_norm), color = 'black', linewidth = 0.50) + 
  geom_point(data = figure_data, aes(x = trial, y = response_mod, shape = factor(response_mod)), size = 2.5, alpha = 0.50, fill = 'darkgrey') + 
  geom_rect(data = viterbi_segments, aes(xmin = start_trial, xmax = end_trial, ymin = 0.05, ymax = 0.48, fill = as.factor(state)), alpha = 0.50) +
  geom_rect(data = richness_segments, aes(xmin = start_trial, xmax = end_trial, ymin = 0.52, ymax = 0.95, fill = as.factor(state+10)), alpha = 0.50) +
  scale_fill_manual(values = c("0" = "pink", "1" = "darkred", "10" = "slateblue1", "11" = "purple4")) + 
  scale_shape_manual(values=c(1, 1)) + 
  scale_x_continuous(name = 'Trial', breaks = seq(0, 200, 50)) + 
  scale_y_continuous(name = 'EV [norm.]', breaks = seq(0, 1, 0.25), limits = c(-0.1, 1.1)) +
  theme_pubr() + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = 'none',
    legend.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    aspect.ratio = 0.35,
  )
fig2a

# figure 2c: model selection with AIC
fig2c <- ggplot(d_aic) + 
  stat_summary(aes(x = as.numeric(n_hmm_states), y = aic, color = monkey, group = monkey), fun = 'mean', geom = 'line', linewidth = 0.5, alpha = 1, linetype = 2, color = 'black', position = position_dodge(width = 0.05)) + 
  stat_summary(aes(x = as.numeric(n_hmm_states), y = aic, fill = monkey, group = monkey), fun.data = 'mean_se', geom = 'point', shape = 21, color = 'black', size = 4, position = position_dodge(width = 0.05)) + 
  theme_pubr() + 
  scale_x_continuous(name = 'HMM-states [N]') + 
  scale_y_continuous(name = 'AIC-score') + 
  scale_color_brewer(name = 'Monkey', palette = 'Reds', labels = c('M1', 'M2', 'M3', 'M4')) + 
  scale_fill_brewer(name = 'Monkey', palette = 'Reds',labels = c('M1', 'M2', 'M3', 'M4')) + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = c(0.22, 0.725),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    aspect.ratio = 1.2,
  ) + 
  border(color = 'black') + 
  coord_cartesian(xlim = c(0.8, 5.2))
fig2c

# figure 2d: state-specific bias coefficients
d_alpha$monkey <- factor(d_alpha$monkey, levels = c('Ultra', 'Ulrich', 'Vampire', 'Winky'),
                       labels = c('M1', 'M2', 'M3', 'M4'))
fig2d <- ggplot(d_alpha %>% 
                  filter(alpha_value > -6.5), aes(y = alpha_value, x = hmm_state, fill = hmm_state, color = hmm_state)) + 
  stat_halfeye(adjust = 0.75, justification = 0, .width = 0, point_colour = NA, alpha = 0.60) + 
  geom_boxplot(width = 0.1, alpha = 0.60, outlier.color = NA, color = 'black') + 
  geom_half_point(data = d_alpha %>% filter(alpha_value > -6.5) %>% sample_n(size = 1e4), aes(y = alpha_value, x = hmm_state, fill = hmm_state, color = hmm_state),
                  side = "l", range_scale = .4, size = 0.05, alpha = 0.2) + 
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.75, color = 'red') + 
  labs(
    x = 'Motivation-state', 
    y = TeX("bias-coeff. [a.u.]")
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  scale_x_discrete(labels = c('Low', 'High')) + 
  scale_fill_manual(name = 'Motivation-state', values = c('pink', 'darkred'), labels = c('Low', 'High')) + 
  scale_color_manual(name = 'Motivation-state', values = c('pink', 'darkred'), labels = c('Low', 'High')) + 
  facet_wrap(~monkey, nrow = 2, scales = c("free_y")) + 
  theme_pubr() + 
  border(color = 'black') + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.text = element_text(size = 14),
    #axis.title.x = element_text(size = 18),
    #axis.text.x = element_text(size = 14, angle = 35, vjust = 0.75),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 14),
    strip.background = element_rect(colour='black', fill='grey99'),
    legend.position = 'none',
    legend.background = element_blank(),
    aspect.ratio = 0.8,
    panel.spacing.x = unit(0.2, "lines")
  )
fig2d

# fig2e: transition matrices
d_tm$monkey <- factor(d_tm$monkey, levels = c('Ultra', 'Ulrich', 'Vampire', 'Winky'),
                         labels = c('M1', 'M2', 'M3', 'M4'))
fig2e <- ggplot(d_tm, aes(x = s_t, y = s_t_plus_1, fill = prob)) + 
  geom_tile(alpha = 0.75, color = 'black', linewidth = 0.25) + 
  geom_text(aes(label = round(prob, 2)), size = 4) +
  scale_fill_continuous(low = 'pink', high = 'darkred', limits = c(0.00, 1.00), breaks = c(seq(0, 1, 0.25)), labels = c('0.00', '0.25', '0.50', '0.75', '1.00')) + 
  coord_fixed() + 
  scale_x_discrete(labels = c('Low', 'High')) + 
  scale_y_discrete(labels = c('Low', 'High')) + 
  labs(x = 'State [t]', y = ' State [t+1]', fill = "P(state[t+1] | state[t])") + 
  facet_wrap(~monkey, nrow = 2) + 
  theme_pubr() + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text = element_text(size = 14),
    #legend.text = element_text(size = 14),
    axis.text.y = element_text(size = 14, angle = 90),
    strip.text = element_text(size = 14),
    strip.background = element_rect(colour='black', fill='grey99'),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.position = 'right',
    aspect.ratio = 1,
    panel.spacing.x = unit(0.2, "lines"),
    legend.key.height = unit(1, 'cm'),
    legend.key.width = unit(0.5, 'cm')
  )
fig2e

# figure 2f: pursue-rates in decoded HMM states
fig2f <- ggplot(d %>%
                  mutate(ev_ind = cut(ev_z, breaks = unique(quantile(ev_z, probs = seq(0, 1, length.out = 20))), labels = FALSE, include.lowest = TRUE)) %>%
                  group_by(ev_ind) %>%
                  mutate(ev_graph = mean(ev_z, na.rm = T)) %>%
                  ungroup() %>%
                  group_by(monkey, ev_graph, viterbi_state_bin) %>% summarise(m = mean(response_mod, na.rm = T)*100), 
                aes(x = ev_graph, y = m, group = factor(viterbi_state_bin), fill = factor(viterbi_state_bin))) + 
  geom_point(shape = 21, alpha = 0.80, color = 'black', position = position_jitter(width = 0.05, height = 0.05)) + 
  geom_smooth(method = 'lm', linewidth = 0.65, color = 'black') + 
  scale_y_continuous(name = 'Pursue-rate [%]', limits = c(10, 100), breaks = c(20, 40, 60, 80, 100)) + 
  scale_x_continuous(name = 'EV [Z]', breaks = c(-1.0, 0, 1.0), labels = number_format(accuracy = 1)) + 
  scale_fill_manual(name = 'Motivation-state', labels = c('Low', 'High'), values = c('pink', 'darkred')) + 
  theme_pubr() + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.background = element_blank(),
    legend.position = 'none',
    aspect.ratio = 1
  ) + 
  coord_cartesian(xlim=c(-1.4, 1.4))
fig2f

# figure 2g: RTs in decoded states
fig2h <- ggplot(d %>% drop_na(RT_z) %>%
                  group_by(monkey, session, viterbi_state) %>%
                  summarise(m = mean(RT_z, na.rm = T)),
                aes(x = factor(viterbi_state), y = m, fill = factor(viterbi_state))) + 
  geom_hline(yintercept = 0, color = 'red', linetype = 2, linewidth = 0.65, alpha = 0.50) + 
  stat_halfeye(adjust = 1, justification = 0, .width = 0, point_colour = NA, alpha = 0.75) + 
  geom_boxplot(width = 0.1, alpha = 0.60, outlier.color = NA, fill = NA, color = 'black') + 
  geom_half_point(side = "l", range_scale = .4, size = 1.5, alpha = 0.50, shape = 21, color = 'black') + 
  scale_y_continuous(name = 'log(RT) [Z]', limits = c(-0.60, 1.5), breaks = seq(-.5, 2, .5)) + 
  scale_x_discrete(name = 'Motivation-state', labels = c('Low', 'High')) + 
  scale_fill_manual(name = 'State', labels = c('Low', 'High'), values = c('pink', 'darkred')) + 
  scale_color_manual(name = 'State', labels = c('Low', 'High'), values = c('pink', 'darkred')) + 
  theme_pubr() + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.background = element_blank(),
    legend.position = 'none',
    aspect.ratio = 1
  )
fig2h

# figure 2h: pupil-size in decoded states
ps_resid <- d %>%
  drop_na(pupil_size_decision_z, prev_response_mod) %>%
  mutate(ps_resid = resid(lm(pupil_size_decision_z ~ mag_z + prob_z + trial_z + prev_response_mod, data = d)))

# Create the plot
fig2i <- ggplot(ps_resid %>%
         group_by(monkey, session, viterbi_state) %>%
         summarise(m = mean(ps_resid, na.rm = T)), aes(x = factor(viterbi_state), y = m, fill = factor(viterbi_state))) +
  geom_hline(yintercept = 0, color = 'red', linetype = 2, linewidth = 0.65, alpha = 0.50) + 
  stat_halfeye(adjust = 1, justification = 0, .width = 0, point_colour = NA, alpha = 0.75) + 
  geom_boxplot(width = 0.1, alpha = 0.60, outlier.color = NA, fill = NA, color = 'black') + 
  geom_half_point(side = "l", range_scale = .4, size = 1.5, alpha = 0.50, shape = 21, color = 'black') + 
  scale_y_continuous(name = 'Pupil-size [Z]') +
  scale_x_discrete(name = 'Motivation-state', labels = c('Low', 'High')) + 
  scale_fill_manual(name = 'State', labels = c('Low', 'High'), values = c('pink', 'darkred')) + 
  theme_pubr() +
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.background = element_blank(),
    legend.position = 'none',
    aspect.ratio = 1
  )
fig2i

# fig 2j: state occupancy as a function of availability of rewards (i.e. average expected value)
fig2j <- ggplot(d %>%
                  mutate(ave_ev_graph = cut(ave_ev_z, breaks = unique(quantile(ave_ev_z, probs = seq(0, 1, length.out = 20), na.rm = T)), labels = FALSE, include.lowest = TRUE)) %>%
                  group_by(ave_ev_graph) %>%
                  mutate(ave_ev_graph = mean(ave_ev_z)) %>%
                  group_by(monkey, ave_ev_graph) %>% 
                  summarise(m = mean(viterbi_state_bin, na.rm = T)), aes(x = ave_ev_graph, y = m, fill = m)) + 
  geom_point(shape = 21, position = position_jitter(width = 0.10), color = 'black', alpha = 0.65, fill = 'pink') + 
  geom_smooth(method = 'lm', color = 'black', linewidth = 0.65, alpha = 0.50, fill = 'pink') + 
  geom_hline(yintercept = 0.5, linetype = 2, color = 'red') + 
  scale_y_continuous(name = 'P(high-motivation)', breaks = seq(0, 1, 0.25), limits = c(0.25, 1.0)) + 
  scale_x_continuous(name = 'Ave.EV [Z]', breaks = seq(-1, 1, 1)) + 
  theme_pubr() + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.background = element_blank(),
    legend.position = 'none',
    aspect.ratio = 1
  ) + 
  coord_cartesian(ylim = c(0.20, 0.90))
fig2j

# fig 2k: average expected value as a function of state-transitions
# Create a vector of NAs with the same length as 'd'
increase_time <- rep(NA, nrow(d))
decrease_time <- rep(NA, nrow(d))

# Find the row indices where state_increase == 1 and state_decrease == 1
increase_indices <- which(d$state_increase == 1)
decrease_indices <- which(d$state_decrease == 1)

# Create 'transition_time' centered around each 'state_transition' == 1
for (index in increase_indices) {
  increase_time[index + (-5:5)] <- -5:5
}

for (index in decrease_indices) {
  decrease_time[index + (-5:5)] <- -5:5
}

# Add 'transition_time' as new columns in the data frame 'd'
d$increase_time <- increase_time
d$decrease_time <- decrease_time

# Combine the results for increases and decreases, adding a variable 'transition_type'
mean_se_data <- rbind(
  d %>%
    group_by(transition_time = increase_time) %>%
    summarize(
      mean_ave_ev_z = mean(ave_ev_z, na.rm = TRUE),
      se_ave_ev_z = sd(ave_ev_z, na.rm = TRUE) / sqrt(sum(!is.na(ave_ev_z))),
      transition_type = "increase"
    ),
  d %>%
    group_by(transition_time = decrease_time) %>%
    summarize(
      mean_ave_ev_z = mean(ave_ev_z, na.rm = TRUE),
      se_ave_ev_z = sd(ave_ev_z, na.rm = TRUE) / sqrt(sum(!is.na(ave_ev_z))),
      transition_type = "decrease"
    )
)

fig2k <- ggplot(data = mean_se_data, aes(x = transition_time, y = mean_ave_ev_z, fill = transition_type)) + 
  geom_vline(xintercept = 0, linetype = 2, color = 'red') + 
  geom_line(color = 'black', linewidth = 0.65) +
  geom_ribbon(aes(x = transition_time, ymin = mean_ave_ev_z - se_ave_ev_z, ymax = mean_ave_ev_z + se_ave_ev_z), alpha = 0.5) +
  scale_fill_manual(name = 'Transition-type', values = c('pink', 'firebrick4'), labels = c('High-to-low', 'Low-to-high')) + 
  scale_x_continuous(name = 'Transition-time [T]', breaks = seq(-4, 4, 2)) + 
  scale_y_continuous(name = 'Ave.EV [Z]', labels = scales::number_format(accuracy = 0.01)) + 
  theme_pubr() + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(size = 14),
    legend.background = element_blank(),
    legend.position = 'none',
    aspect.ratio = 1
  ) + 
  coord_cartesian(ylim = c(-0.5, 0.5))
fig2k

# fig 2l: average expected value as a function of state-transitions
# Create a vector of NAs with the same length as 'd'
increase_time <- rep(NA, nrow(d))
decrease_time <- rep(NA, nrow(d))

# Find the row indices where state_increase == 1 and state_decrease == 1
increase_indices <- which(d$state_increase == 1)
decrease_indices <- which(d$state_decrease == 1)

# Create 'transition_time' centered around each 'state_transition' == 1
for (index in increase_indices) {
  increase_time[index + (-5:5)] <- -5:5
}

for (index in decrease_indices) {
  decrease_time[index + (-5:5)] <- -5:5
}

# Add 'transition_time' as new columns in the data frame 'd'
d$increase_time <- increase_time
d$decrease_time <- decrease_time

# Combine the results for increases and decreases, adding a variable 'transition_type'
mean_se_data <- rbind(
  d %>%
    group_by(transition_time = increase_time) %>%
    summarize(
      mean_response = mean(response, na.rm = TRUE),
      se_response = sd(response, na.rm = TRUE) / sqrt(sum(!is.na(response))),
      transition_type = "increase"
    ),
  d %>%
    group_by(transition_time = decrease_time) %>%
    summarize(
      mean_response = mean(response, na.rm = TRUE),
      se_response = sd(response, na.rm = TRUE) / sqrt(sum(!is.na(response))),
      transition_type = "decrease"
    )
)

fig2l <- ggplot(data = mean_se_data, aes(x = transition_time, y = mean_response*100, fill = transition_type, color = transition_type)) + 
  geom_vline(xintercept = 0, linetype = 2, color = 'red') + 
  geom_line(color = 'black', linewidth = 0.65) +
  geom_ribbon(aes(x = transition_time, ymin = (mean_response - se_response)*100, ymax = (mean_response + se_response)*100), alpha = 0.5) +
  scale_fill_manual(name = 'Transition-type', values = c('pink', 'firebrick4'), labels = c('High-to-low', 'Low-to-high')) + 
  scale_color_manual(name = 'Transition-type', values = c('pink', 'firebrick4'), labels = c('High-to-low', 'Low-to-high')) + 
  scale_x_continuous(name = 'Transition-time [T]', breaks = seq(-4, 4, 2)) + 
  scale_y_continuous(name = 'Pursue-rate [%]', breaks = seq(0, 100, 25)) + 
  theme_pubr() + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(size = 14),
    legend.background = element_blank(),
    legend.position = c(0.725, 0.5),
    legend.title = element_blank(),
    aspect.ratio = 0.75
  ) + 
  coord_cartesian(ylim = c(-20, 120))
fig2l

######################## SUPPLEMENTARY FIGURES #################### 

# fig S2C: graph mean state occupancy as a function of time
bin_width <- 10
figs2b <- ggplot(d %>% 
         mutate(bin = findInterval(trial, seq(min(trial), max(trial) + bin_width, by = bin_width))) %>% 
         group_by(bin, session) %>% 
         summarize(mean_viterbi = mean(viterbi_state_bin, na.rm = TRUE)), 
       aes(x = bin, y = mean_viterbi)) + 
  geom_jitter(width = 0.5, height = 0.05, size = 1, shape = 21, fill = 'pink', color = 'pink', alpha = 0.75) + 
  geom_hline(yintercept = 0.56, linetype = 2, color = 'red') + 
  stat_summary(fun.data = 'mean_se', geom = 'ribbon', fill = 'darkred', alpha = 0.60) + 
  stat_summary(fun.data = 'mean_se', geom = 'line', color = 'black') + 
  scale_x_continuous(name = 'Trial [N]', breaks = seq(0, 20, 5), limits = c(0, 20), labels = seq(from = 0, to = 200, by = 50)) + 
  scale_y_continuous(name = 'p(HMM-state==High)', breaks = seq(0, 1, 0.25), limits = c(0, 1)) + 
  theme_pubr() + 
    theme(
      axis.title.y = element_text(size = 18),
      axis.title.x = element_text(size = 18),
      axis.text = element_text(size = 14),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      axis.text.y = element_text(size = 14),
      legend.background = element_blank(),
      legend.position = c(0.25, 0.975),
      aspect.ratio = 0.85
    )
figs2b

# fig S2C: graph likelihood of state transition as a function of time
figs2c <- ggplot(d %>% 
         mutate(bin = findInterval(trial, seq(min(trial), max(trial) + bin_width, by = bin_width))) %>% 
         group_by(bin, session) %>% 
         summarize(p_change = mean(state_transition, na.rm = TRUE)), 
       aes(x = bin, y = p_change)) + 
  #geom_jitter(width = 0.5, height = 0.05, size = 1, shape = 21, fill = 'pink', color = 'pink', alpha = 0.75) + 
  geom_hline(yintercept = 0.016, linetype = 2, color = 'red') + 
  stat_summary(fun.data = 'mean_se', geom = 'ribbon', fill = 'grey', alpha = 0.60) + 
  stat_summary(fun.data = 'mean_se', geom = 'line', color = 'black') + 
  scale_x_continuous(name = 'Trial [N]', breaks = seq(0, 20, 5), limits = c(0, 20), labels = seq(from = 0, to = 200, by = 50)) + 
  scale_y_continuous(name = 'p(State-transition)') + 
  theme_pubr() + 
    theme(
      axis.title.y = element_text(size = 18),
      axis.title.x = element_text(size = 18),
      axis.text = element_text(size = 14),
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      axis.text.y = element_text(size = 14),
      legend.background = element_blank(),
      legend.position = c(0.25, 0.975),
      aspect.ratio = 0.7
    )
figs2c

# fig S2D: count of state transitions per session
figs2d <- ggplot(d %>% 
         group_by(monkey, session) %>% 
         summarise(n = sum(state_transition, na.rm = TRUE)), 
       aes(x = n)) + 
  geom_histogram(binwidth = 1, fill = 'pink', alpha = 0.5, color = 'black') + 
  scale_y_continuous(name = 'Count', limits = c(0, 23), breaks = seq(0, 20, 5)) + 
  scale_x_continuous(name = 'State-transitions [N]', breaks = seq(0, 10, 2)) + 
  theme_pubr() + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(size = 14),
    legend.background = element_blank(),
    aspect.ratio = 1.5)
figs2d
