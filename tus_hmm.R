# Code for all analysis involving GLM-HMM for behavioural data recorded during 
# TUS sessions ("Non-invasive disruption of DRN perturbs motivation-state 
# transitions "), Figure 4F-H ("Non-invasive perturbation of dorsal raphe 
# nucleus shows that it is causally involved in the influence of a monkey’s 
# environment on its behaviour"). 

rm(list = ls())

# load packages
library(tidyr)
library(lme4)
library(rstan)
library(dplyr)
library(ggplot2)
library(gghalves)
library(ggdist)
library(ggpubr)

################### Read data ################
behaviour_path <- "/Volumes/LaCie/nature_neuro_submission/data/behaviour/tus_experiment/tus_behaviour.csv"
d <- read.csv(behaviour_path)
d$stim_cond <- factor(d$stim_cond, levels = c('sham', 'sts', 'drn', 'vta'), labels = c('Sham', 'STS', 'DRN', 'VTA')) # code TUS-condition such that sham is the reference group

################### specify functions ################

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

################## perform viterbi decoding ################
glm_hmm_dir <- "/Volumes/LaCie/nature_neuro_submission/data/behaviour/tus_experiment/glm_hmm/"
d$viterbi_state <- NA
n_states <- 2

for(cond in unique(d$stim_cond)){
  
  # get hmm-glm
  hmm <- readRDS(paste(glm_hmm_dir,cond, '_', n_states, '_state_tus.rds', sep = ''))
  
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
  
  for(s in unique(d$session_id[d$stim_cond==cond])){
    
    y <- d$response_mod[d$stim_cond==cond & d$session_id==s]
    ev_z <- d$ev_z[d$stim_cond==cond & d$session_id==s]
    richness <- d$richness_bin[d$stim_cond==cond & d$session_id==s]
    stochasticity <- d$stochasticity_bin[d$stim_cond==cond & d$session_id==s]
    
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
    d$viterbi_state[d$stim_cond==cond & d$session_id==s] <- zstar
    rm(y, ev_z, richness, stochasticity, zstar)
  }
  
  rm(hmm, A, pi, alpha, beta_ev, beta_richness, beta_stochasticity)
  
}

# identify state changes
d$state_difference <- NA # Initialize the state_difference column with NA

# Loop through unique combinations of 'monkey' and 'session_id'
unique_monkeys <- unique(d$monkey)
unique_sessions <- unique(d$session_id)

for (m in unique_monkeys) {
  for (s in unique_sessions) {
    subset_data <- d[d$monkey == m & d$session_id == s, ]
    
    if (nrow(subset_data) > 1) {
      subset_data$state_difference <- c(NA, diff(subset_data$viterbi_state))
      d[d$monkey == m & d$session_id == s, ] <- subset_data
    }
    rm(subset_data)
  }
}

d$state_transition <- ifelse(d$state_difference!=0, 1, 0)
d$state_increase <- ifelse(d$state_difference==1, 1, 0)
d$state_decrease <- ifelse(d$state_difference==-1, 1, 0)

################### GLMs over Viterbi states ######################
d$viterbi_state_bin <- ifelse(d$viterbi_state==1, 0, 
                                           ifelse(d$viterbi_state==2, 1, NA))

GLM4.3 <- glmer(state_transition ~ mag_z*prob_z + stim_cond + trial_z + (1|monkey), data = d, family = 'binomial')
summary(GLM4.3)

# GLM4.3: DRN-vs-VTA
GLM4.3a <- glmer(state_transition ~ mag_z*prob_z + stim_cond + trial_z + (1|monkey), data = d %>% filter(stim_cond=='DRN'|stim_cond=='VTA'), family = 'binomial')
summary(GLM4.3a)

# GLM4.3: DRN-vs-STS
GLM4.3b <- glmer(state_transition ~ mag_z*prob_z + stim_cond + trial_z + (1|monkey), data = d %>% filter(stim_cond=='DRN'|stim_cond=='STS'), family = 'binomial')
summary(GLM4.3b)

GLM4.4 <- glmer(state_decrease ~ mag_z*prob_z + stim_cond + trial_z + (1|monkey), data = d, family = 'binomial')
summary(GLM4.4)

GLM4.5 <- glmer(viterbi_state_bin ~ ev_z + ave_ev_z*stim_cond + trial_z + (1|monkey), data = d, family = 'binomial')
summary(GLM4.5)

# GLM4.5: DRN-vs-STS
GLM4.5a <- glmer(viterbi_state_bin ~ ev_z + ave_ev_z*stim_cond + trial_z + (1|monkey), data = d %>% filter(stim_cond == 'DRN'|stim_cond == 'STS'), family = 'binomial')
summary(GLM4.5a)

# GLM4.5: DRN-vs-VTA
GLM4.5b <- glmer(viterbi_state_bin ~ ev_z + ave_ev_z*stim_cond + trial_z + (1|monkey), data = d %>% filter(stim_cond == 'DRN'|stim_cond == 'VTA'), family = 'binomial')
summary(GLM4.5b)

# GLM4.5: low-EV environments (i.e. z-scored average-EV < 0)
GLM4.5c <- glmer(viterbi_state_bin ~ ev_z + ave_ev_z*stim_cond + trial_z + (1|monkey), data = d %>% filter(ave_ev_z < 0), family = 'binomial')
summary(GLM4.5c)

# GLM4.5: high-EV environments (i.e. z-scored average-EV >= 0)
GLM4.5d <- glmer(viterbi_state_bin ~ ev_z + ave_ev_z*stim_cond + trial_z + (1|monkey), data = d %>% filter(ave_ev_z >= 0), family = 'binomial')
summary(GLM4.5d)

# GLM4.6
GLM4.6 <- glmer(viterbi_state_bin ~ ev_z + richness_bin*stim_cond + stochasticity_bin + trial_z + (1|monkey), data = d, family = 'binomial')
summary(GLM4.6)

GLM4.6 <- glmer(viterbi_state_bin ~ ev_z + richness_bin*stim_cond + stochasticity_bin + trial_z + (1|monkey), data = d %>% filter(stim_cond %in% c('DRN', 'VTA')), family = 'binomial')
summary(GLM4.6)

################ FIGURES #####################

# Figure 4F: posterior distribution of transition probabilities for different TUS conditions
tus_conditions <- unique(d$stim_cond)
n_states <- 2
n_repeats <- 2e4

# Create a data frame to store results
d_tm <- data.frame(
  cond = rep(tus_conditions, each = n_repeats * (n_states^2)),
  transition = NA,
  sample = NA
)

for (cond in tus_conditions) {
  model_name <- paste(glm_hmm_dir,cond, '_', n_states, '_state_tus.rds', sep = '')
  hmm <- readRDS(model_name)
  tm_samples <- rstan::extract(hmm, pars = 'A')
  sample_matrix <- matrix(unlist(tm_samples), ncol = 4, byrow = FALSE)
  sample_matrix <- reshape2::melt(sample_matrix) %>% select(-Var1)
  d_tm[d_tm$cond == cond, c("transition", "sample")] <- sample_matrix
  rm(tm_samples, sample_matrix, hmm)
}

d_tm <- d_tm %>%
  mutate(transition = case_when(
    transition == 1 ~ 'Low-to-low',
    transition == 2 ~ 'High-to-low',
    transition == 3 ~ 'Low-to-high',
    transition == 4 ~ 'High-to-high',
    TRUE ~ NA_character_
  ))

pink_to_red_palette <- colorRampPalette(c('white', 'darkorchid4', "darkred"))(4)
fig4f <- ggplot(d_tm %>% filter(transition=='Low-to-high'|transition=='High-to-low') %>%
         mutate(cond = factor(cond, levels = c('Sham', 'STS', 'VTA', 'DRN'))), aes(x = sample, group = cond, fill = cond)) + 
  geom_density(adjust = 3, color = 'black', alpha = 0.75) + 
  scale_x_continuous(name = 'p(Transition)', breaks = seq(0, 0.06, 0.02), limits = c(0, 0.052)) +
  scale_y_continuous(name = 'Posterior density [a.u.]', limits = c(0, 115), breaks = seq(0, 90, 30)) + 
  scale_fill_manual(name = 'TUS condition', 
                    values = pink_to_red_palette) + 
  #guides(fill = 'none') + 
  theme_pubr() + 
  border(color = 'black') + 
  theme(text = element_text(size = 16), axis.text = element_text(size = 14), panel.spacing = unit(1, "lines"), aspect.ratio = 0.7) + 
  facet_grid(~transition)
fig4f

# Figure 4G: count of decoded motivation-state transitions
pink_to_red_palette <- colorRampPalette(c('white', 'darkorchid4', "darkred"))(4)
fig4g <- ggplot(d %>% 
                  mutate(stim_cond = factor(stim_cond, levels = c('Sham', 'STS', 'VTA', 'DRN'))) %>% 
                  group_by(stim_cond) %>%
                  summarise(n_transition = sum(state_transition, na.rm = T)), 
                aes(x = stim_cond, fill = stim_cond, y = n_transition)) + 
  geom_bar(stat = 'identity', position = 'dodge', color = 'black', alpha = 0.75, width = 0.75) + 
  theme_pubr() + 
  scale_x_discrete(name = 'TUS condition') + 
  scale_y_continuous(name = 'Transitions (N)',  breaks = seq(0, 40, 20), limits = c(0, 52)) + 
  scale_fill_manual(values = pink_to_red_palette) + 
  theme(legend.position = 'none') + 
  theme(text = element_text(size = 16), axis.text = element_text(size = 14), aspect.ratio = 1)
fig4g

# Figure 4H: relationship between motivation-state occupancy and availability of rewards (Ave. EV)
fig4h <- ggplot(d %>%
         filter(stim_cond == 'DRN'|stim_cond == 'Sham') %>%
         drop_na(ave_ev_z) %>%
         mutate(ave_ev_graph = cut(ave_ev_z, breaks = quantile(ave_ev_z, probs = seq(0, 1, length.out = 11), na.rm = T), labels = FALSE, include.lowest = TRUE)) %>%
         group_by(ave_ev_graph) %>%
         mutate(ave_ev_graph = mean(ave_ev_z, na.rm = T)) %>%
         group_by(stim_cond, ave_ev_graph) %>% 
         summarise(m = mean(viterbi_state_bin, na.rm = T)), aes(x = ave_ev_graph, y = m, fill = stim_cond, color = stim_cond)) + 
  stat_summary(fun.data = 'mean_se', geom = 'pointrange', color = 'black', size = 0.4, shape = 21) + 
  geom_smooth(method = 'lm', color = 'black', linewidth = 0.65, alpha = 0.50) + 
  scale_y_continuous(name = 'P(High-motivation)', limits = c(0, 1), breaks = seq(0, 1, 0.25)) + 
  scale_x_continuous(name = 'Ave.EV [Z]', breaks = seq(-2, 2, 2), limits = c(-2.2, 2.2)) + 
  scale_color_manual(name = 'TUS condition', values = c('thistle2', 'darkred'), labels = c('Sham', 'DRN')) + 
  scale_fill_manual(name = 'TUS condition', values = c('thistle2', 'darkred'), labels = c('Sham', 'DRN')) + 
  theme_pubr() + 
  coord_cartesian(ylim = c(0.15, 0.75)) + 
  theme(strip.text = element_blank(), text = element_text(size = 16), axis.text = element_text(size = 14), panel.spacing = unit(0.25, "lines"), aspect.ratio = 2, legend.position = 'top', legend.direction = 'horizontal', legend.title = element_blank()) + 
  facet_grid(~stim_cond) + 
  border(color = 'black')
fig4h

fig4h <- ggplot(d %>%
                  filter(stim_cond == 'DRN'|stim_cond == 'Sham') %>%
                  drop_na(ave_ev_z) %>%
                  mutate(ave_ev_graph = cut(ave_ev_z, breaks = quantile(ave_ev_z, probs = seq(0, 1, length.out = 11), na.rm = T), labels = FALSE, include.lowest = TRUE)) %>%
                  group_by(ave_ev_graph) %>%
                  mutate(ave_ev_graph = mean(ave_ev_z, na.rm = T)), aes(x = ave_ev_graph, y = viterbi_state_bin, fill = stim_cond, color = stim_cond)) + 
  stat_summary(fun.data = 'mean_se', geom = 'pointrange', color = 'black', size = 0.4, shape = 21) + 
  geom_smooth(method = 'lm', color = 'black', linewidth = 0.65, alpha = 0.50) + 
  scale_y_continuous(name = 'P(High-motivation)', limits = c(0, 1), breaks = seq(0, 1, 0.25)) + 
  scale_x_continuous(name = 'Ave.EV [Z]', breaks = seq(-1.5, 1.5, 1.5)) + 
  scale_color_manual(name = 'TUS condition', values = c('thistle2', 'darkred'), labels = c('Sham', 'DRN')) + 
  scale_fill_manual(name = 'TUS condition', values = c('thistle2', 'darkred'), labels = c('Sham', 'DRN')) + 
  theme_pubr() + 
  coord_cartesian(ylim = c(0.24, 0.50), xlim = c(-2, 2)) + 
  theme(strip.text = element_blank(), text = element_text(size = 16), axis.text = element_text(size = 14), panel.spacing = unit(0.25, "lines"), aspect.ratio = 1.5, legend.position = 'none', legend.direction = 'horizontal', legend.title = element_blank()) + 
  facet_grid(~stim_cond)
fig4h

