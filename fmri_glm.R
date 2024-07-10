# Code for all analysis in GLM1 ("Motivation for rewards is modulated by the 
# environment), Figure 1B-H ("Animals are more likely to pursue rewards in rich 
# environments), and Supplementary Figure 1 ("Additional analysis of behaviour)

rm(list = ls())

# load packages
library(lme4)
library(dplyr)
library(ggplot2)
library(gghalves)
library(ggdist)
library(ggpubr)
library(tidyverse)
library(scales)

################### Read data ################
behaviour_path <- "/Volumes/LaCie/nature_neuro_submission/data/behaviour/fmri_experiment/fmri_behaviour.csv"
d <- read.csv(behaviour_path)

################### TIME-HORIZON ANALYSIS ####################### 
monkeys <- unique(d$monkey)
d$rw_1back <- NA
d$rw_2back <- NA
d$rw_3back <- NA
d$rw_4back <- NA
d$rw_5back <- NA
d$rw_6back <- NA
d$rw_7back <- NA
d$rw_8back <- NA
d$rw_9back <- NA
d$rw_10back <- NA

model_coefficients <- data.frame(
  predictor=character(),
  coefficient=numeric(),
  se=numeric(),
  p_val=numeric()
)

for(m in monkeys){
  monkey_data <- d %>% filter(monkey==m)
  sessions <- unique(monkey_data$session)
  for(s in sessions){
    session_data <- monkey_data %>% filter(session==s)
    trials <- max(session_data$trial)
    for(t in 11:trials){
      session_data$rw_1back[t] <- session_data$rw_received[t-1]
      session_data$rw_2back[t] <- session_data$rw_received[t-2]
      session_data$rw_3back[t] <- session_data$rw_received[t-3]
      session_data$rw_4back[t] <- session_data$rw_received[t-4]
      session_data$rw_5back[t] <- session_data$rw_received[t-5]
      session_data$rw_6back[t] <- session_data$rw_received[t-6]
      session_data$rw_7back[t] <- session_data$rw_received[t-7]
      session_data$rw_8back[t] <- session_data$rw_received[t-8]
      session_data$rw_9back[t] <- session_data$rw_received[t-9]
      session_data$rw_10back[t] <- session_data$rw_received[t-10]
    }
    monkey_data[monkey_data$session==s,] <- session_data
  }
  d[d$monkey==m,] <- monkey_data
}

base_formula_str <- 'response ~ mag_z*prob_z + rw_1back + '

for(trial in 1:10){
  random_slope <- paste('rw_',trial,'back', sep  = '')
  if(trial==1){
    formula_str <- paste(base_formula_str, '(0 +', random_slope, '|monkey)')
  }else{
    base_formula_str <- paste(base_formula_str, random_slope, ' + ', sep = '')
    formula_str <- paste(base_formula_str, '(0 +', random_slope, '|monkey)')
  }
  formula <- as.formula(formula_str)
  model <- glmer(formula, data = d, family = 'binomial')
  coefficient <- fixef(model)[(3+trial)]
  se <- coef(summary(model))[(3+trial),2]
  p_val <- coef(summary(model))[(3+trial),4]
  model_coefficients <- rbind(model_coefficients, 
                              data.frame(
                                predictor = names(coefficient),
                                se=as.vector(se),
                                t_stat = as.vector(coefficient),
                                p_val = as.vector(p_val))
  )
} 

model_coefficients$predictor <- factor(model_coefficients$predictor,
                                       levels = c("rw_1back",
                                                  "rw_2back",
                                                  "rw_3back",
                                                  "rw_4back",
                                                  "rw_5back",
                                                  "rw_6back",
                                                  "rw_7back",
                                                  "rw_8back",
                                                  "rw_9back",
                                                  "rw_10back"),
                                       labels = c('t-1',
                                                  't-2',
                                                  't-3',
                                                  't-4',
                                                  't-5',
                                                  't-6',
                                                  't-7',
                                                  't-8',
                                                  't-9',
                                                  't-10'))

fig1d <- ggplot(model_coefficients, aes(x = predictor, y = t_stat, ymin = t_stat-(1.96*se), ymax = t_stat + (1.96*se))) + 
  geom_line(aes(group = 1), alpha = 1, color = 'black', linetype = 2) +
  geom_pointrange(shape = 21, size = 1, fill = 'lightblue', color = 'black') + 
  geom_hline(yintercept = 0, color = 'blue', linetype = 2) + 
  scale_y_continuous(name = 'Effect-size [a.u.]', breaks = seq(0, 0.6, 0.2)) + 
  scale_x_discrete(name = 'Outcome', breaks = c('t-1', 't-5', 't-10')) + 
  theme_pubr() + 
  theme(    text = element_text(size = 18),
            axis.text = element_text(size = 14),
            axis.text.x = element_text(size = 14, angle = 40, vjust = 0.5),
            legend.text = element_text(size = 14),
            legend.title = element_blank(),
            aspect.ratio = 0.75,
            legend.background = element_blank(), 
            legend.position = c(0.8, 0.75)) + 
  border('black') + 
  coord_cartesian(ylim = c(-0.1, 0.7), xlim = c(0, 11))
fig1d

################### GLMS ################
GLM1.2a <- glmer(response_mod ~ mag_z*prob_z + ave_rw_z + sd_ave_rw_z + trial_z + (ave_rw_z|monkey), data = d, family = 'binomial')
summary(GLM1.1a)

GLM1.2b <- glmer(response_mod ~ mag_z*prob_z + ave_rw_z + sd_ave_rw_z + trial_z + (sd_ave_rw_z|monkey), data = d, family = 'binomial')
summary(GLM1.1b)

GLM1.3 <- glmer(response_mod ~ mag_z*prob_z + ave_rw_z + sd_ave_rw_z + prev_response_mod + trial_z + (prev_response_mod|monkey), data = d, family = 'binomial')
summary(GLM1.3)

# GLM for supplementary figure S1A
S1A_data <- d %>% filter(stochasticity=='stochastic' & (prob >= 0.60) & (prob <=0.75))
GLM.S1A <- glmer(response_mod ~ mag + richness_bin + trial_z + (richness_bin|monkey), data = d, family = binomial)
summary(GLM.S1A)

# GLMs for supplementary figure S1Bi and S1Bii
GLM.S1Bi <- glmer(response_mod ~ mag_z * prob_z + prev_rw_z + trial_z + (prev_rw_z|monkey), data = d, family = binomial)
summary(GLM.S1Bi)

GLM.S1Bii <- glmer(response_mod ~ mag_z * prob_z + ave_rw_10back_z + trial_z + (prev_rw_z|monkey), data = d, family = binomial)
summary(GLM.S1Bii)

################### FIGURES ################

# Define constants for common plot settings
TEXT_SIZE_MAIN <- 18
TEXT_SIZE_AXIS <- 14
ASPECT_RATIO_STANDARD <- 0.65
ASPECT_RATIO_SQUARE <- 1.0
ALPHA_LEVEL <- 0.80
LINE_WIDTH <- 0.5

# Prepare data for fig. 1b
x <- seq(0, 1.01, length.out = 1e4)
rich_predictable <- dnorm(x, mean = 0.75, sd = 0.05)
rich_stochastic <- dunif(x, min = 0.60, max = 1.00)
poor_predictable <- dnorm(x, mean = 0.55, sd = 0.05)
poor_stochastic <- dunif(x, min = 0.35, max = 0.75)

# Combine data into a data frame
pdf_data <- data.frame(
  x = x,
  rich_pred = rich_predictable,
  rich_stoch = rich_stochastic,
  poor_pred = poor_predictable,
  poor_stoch = poor_stochastic
)

# Reshape data for plotting
pdf_data_long <- pdf_data %>%
  pivot_longer(cols = -x, names_to = 'env', values_to = 'pdf')

# Plot for fig. 1b: PDF over reward environments
fig1b <- ggplot(pdf_data_long, aes(x = x, y = pdf, fill = env)) + 
  geom_area(color = 'black', position = 'identity', alpha = ALPHA_LEVEL) + 
  theme_pubr() + 
  scale_x_continuous(labels = number_format(accuracy = 0.01), breaks = seq(0, 1, 0.20), limits = c(0.20, 1.01)) + 
  labs(x = 'p(Reward)', y = 'Prob. density') + 
  scale_fill_manual(
    values = c("#0b2370", "#B2D4FF", "#6b0f2e", "#e0a2b6"),
    labels = c('Poor-predictable', 'Poor-stochastic', 'Rich-predictable', 'Rich-stochastic'),
    name = 'Environment'
  ) + 
  theme(
    text = element_text(size = TEXT_SIZE_MAIN),
    axis.text = element_text(size = TEXT_SIZE_AXIS),
    legend.background = element_blank(),
    plot.background = element_blank(),
    legend.position = c(0.185, 0.70),
    legend.title = element_text(size = TEXT_SIZE_AXIS),
    legend.text = element_text(size = 11),
    aspect.ratio = 0.5,
  )
fig1b

# Preprocessing for fig. 1c
example_session <- d %>%
  filter(monkey == 'Winky' & session == 14) %>%
  select(trial, prob, richness, stochasticity, ave_rw, sd_ave_rw)

# Mapping environment conditions to labels
env_labels <- c('Poor-Predictable', 'Poor-Stochastic', 'Rich-Predictable', 'Rich-Stochastic')
example_session$env <- factor(
  interaction(example_session$richness, example_session$stochasticity, sep = '-'),
  levels = c('poor-predictable', 'poor-stochastic', 'rich-predictable', 'rich-stochastic'),
  labels = env_labels
)

# Plot for fig. 1c: Example session part 1
fig1c <- ggplot(example_session, aes(x = trial, y = prob)) +
  geom_vline(xintercept = c(47, 97, 148), color = 'red', linetype = 2) +
  geom_ribbon(aes(fill = env, ymax = 1, ymin = 0), linetype = 2, alpha = 0.15) +
  geom_line(linewidth = LINE_WIDTH) +
  theme_pubr() +
  labs(x = 'Trial', y = 'p(Reward)') +
  scale_y_continuous(labels = number_format(accuracy = 0.01)) + 
  scale_fill_manual(
    values = c("skyblue1", "skyblue1", "white", "white"),
    labels = env_labels,
    name = 'Environment'
  ) +
  guides(fill = 'none') +
  coord_cartesian(ylim = c(0.30, 1)) + 
  theme(
    text = element_text(size = TEXT_SIZE_MAIN),
    axis.text = element_text(size = TEXT_SIZE_AXIS),
    aspect.ratio = 0.5,
  ) 
fig1c

# Preprocessing for fig. 1e & fig. 1f: effect of richness and stochasticity
figure_data <- d %>%
  mutate(
    ave_rw_ind = ifelse(ave_rw_z < 0, -1, ifelse(ave_rw_z >= 0, 1, NA)),
    sd_ave_rw_ind = ifelse(sd_ave_rw_z < 0, -1, ifelse(ave_rw_z >= 0, 1, NA))
  )

# Plot for fig. 1e: effect of richness
fig1e <- ggplot(
  figure_data %>%
    drop_na(ave_rw_ind) %>%
    mutate(ev_ind = cut(ev_z, breaks = unique(quantile(ev_z, probs = seq(0, 1, length.out = 12))), labels = FALSE, include.lowest = TRUE)) %>%
    group_by(ev_ind) %>%
    mutate(ev_graph = mean(ev_z, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(monkey, ev_graph, ave_rw_ind) %>% 
    summarise(m = mean(response_mod, na.rm = TRUE) * 100), 
  aes(x = ev_graph, y = m, group = factor(ave_rw_ind), fill = factor(ave_rw_ind))
) +
  stat_summary(fun.data = 'mean_se', geom = 'pointrange', size = 0.75, shape = 21, alpha = ALPHA_LEVEL, color = 'black') + 
  stat_summary(fun.data = 'mean_se', geom = 'pointrange', size = 1, shape = 21, alpha = ALPHA_LEVEL, color = 'black') + 
  geom_smooth(method = 'lm', linewidth = LINE_WIDTH, color = 'black') + 
  scale_y_continuous(name = 'Pursue-rate [%]', limits = c(20, 100), breaks = c(20, 40, 60, 80, 100)) + 
  scale_x_continuous(name = 'EV [Z]', breaks = c(-1.0, 0, 1.0), labels = number_format(accuracy = 1)) + 
  scale_fill_manual(name = 'Richness', labels = c('Low', 'High'), values = c('skyblue1', 'blue4')) + 
  theme_pubr() + 
  theme(
    text = element_text(size = TEXT_SIZE_MAIN),
    axis.text = element_text(size = TEXT_SIZE_AXIS),
    legend.text = element_text(size = TEXT_SIZE_AXIS),
    aspect.ratio = 1,
    legend.background = element_blank(), 
    legend.position = c(0.30, 0.95),
    legend.direction = 'vertical'
  ) +
  coord_cartesian(xlim=c(-1.5, 1.5))
fig1e

# Plot for fig. 1f: effect of stochasticity
fig1f <- ggplot(
  figure_data %>%
    drop_na(sd_ave_rw_ind) %>%
    mutate(ev_ind = cut(ev_z, breaks = unique(quantile(ev_z, probs = seq(0, 1, length.out = 12))), labels = FALSE, include.lowest = TRUE)) %>%
    group_by(ev_ind) %>%
    mutate(ev_graph = mean(ev_z, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(monkey, ev_graph, sd_ave_rw_ind) %>% 
    summarise(m = mean(response_mod, na.rm = TRUE) * 100), 
  aes(x = ev_graph, y = m, group = factor(sd_ave_rw_ind), fill = factor(sd_ave_rw_ind))
) +
  stat_summary(fun.data = 'mean_se', geom = 'pointrange', size = 0.75, shape = 21, alpha = ALPHA_LEVEL, color = 'black') + 
  geom_smooth(method = 'lm', linewidth = LINE_WIDTH, color = 'black') + 
  scale_y_continuous(name = 'Pursue-rate [%]', limits = c(20, 100), breaks = c(20, 40, 60, 80, 100)) + 
  scale_x_continuous(name = 'EV [Z]', breaks = c(-1.0, 0, 1.0), labels = number_format(accuracy = 1)) + 
  scale_fill_manual(name = 'Stochasticity', labels = c('Low', 'High'), values = c('#6b0f2e', 'pink')) + 
  theme_pubr() + 
  theme(
    text = element_text(size = TEXT_SIZE_MAIN),
    axis.text = element_text(size = TEXT_SIZE_AXIS),
    legend.text = element_text(size = TEXT_SIZE_AXIS),
    aspect.ratio = ASPECT_RATIO_SQUARE,
    legend.background = element_blank(), 
    legend.position = c(0.35, 0.95),
    legend.direction = 'vertical'
  ) +
  coord_cartesian(xlim=c(-1.5, 1.5))
fig1f

# Plot for fig. 1g: Behavioural-history
fig1g <- ggplot(
  d %>%
    drop_na(prev_response_mod) %>%
    group_by(monkey, session, prev_response_mod) %>% 
    summarise(m = mean(response_mod, na.rm = TRUE) * 100),
  aes(x = factor(prev_response_mod), fill = factor(prev_response_mod), color = factor(prev_response_mod), y = m)
) +
  stat_summary(fun = 'mean', geom = 'bar', color = 'black', width = 0.65, alpha = 1) + 
  geom_point(size = 2.5, colour = 'grey50', shape = 21, alpha = 0.75) +
  scale_y_continuous(name = 'Pursue-rate[%]', breaks = seq(20, 100, 20)) + 
  scale_x_discrete(name = 'Behavioural hist.', labels = c('Reject', 'Pursue')) + 
  scale_fill_manual(name = 'Behavioural history', values = c('grey90', 'azure2'), labels = c('Reject', 'Pursue')) + 
  scale_color_manual(name = 'Behavioural history', values = c('grey90', 'azure2'), labels = c('Reject', 'Pursue')) + 
  theme_pubr() + 
  theme(
    text = element_text(size = TEXT_SIZE_MAIN),
    axis.text = element_text(size = TEXT_SIZE_AXIS),
    legend.title = element_blank(),
    legend.text = element_text(size = TEXT_SIZE_AXIS),
    aspect.ratio = ASPECT_RATIO_SQUARE,
    legend.background = element_blank(), 
    legend.position = c(0.3, 0.95)
  ) + 
  coord_cartesian(ylim = c(20, 100))
fig1g

######################## SUPPLEMENTARY FIGURES ####################

# supplementary figure S1A: effect of experimentally-defined environments
S1A_data <- d %>% filter(stochasticity=='stochastic' & (prob >= 0.60) & (prob <=0.75))
fig_s1a <- ggplot(S1A_data %>%
                  group_by(monkey, session, richness) %>% 
                  summarise(m = mean(response_mod, na.rm = T)*100),
                aes(x = factor(richness), fill = factor(richness), y = m)) + 
  stat_summary(fun = 'mean', geom = 'bar', color = 'black', width = 0.65, alpha = 0.75) + 
  geom_point(size = 2.5, colour = 'black', shape = 21, alpha = 0.35) +
  scale_y_continuous(name = 'Pursue-rate[%]', breaks = seq(20, 100, 20)) + 
  scale_x_discrete(name = 'Environment', labels = c('Poor', 'Rich')) + 
  scale_fill_manual(name = 'Environment', values = c('skyblue1', 'blue4'), labels = c('Poor', 'Rich')) + 
  theme_pubr() + 
  theme(
    text = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    aspect.ratio = 1,
    legend.background = element_blank(), 
    legend.position = c(0.40, 0.95),
    legend.direction = 'horizontal') + 
  guides(fill = 'none') + 
  coord_cartesian(ylim = c(20, 100))
fig_s1a

# supplementary figure S1B: testing different time-horizons for richness of env.

# Scale and pivot data for fig. s1bi
figure_data <- d %>%
  mutate(
    prev_rw_ind = ifelse(prev_rw_z < 0, -1, ifelse(prev_rw_z >= 0, 1, NA)),
  )
figs1bi <- ggplot(
  figure_data %>%
    drop_na(prev_rw_ind) %>%
    mutate(ev_ind = cut(ev_z, breaks = unique(quantile(ev_z, probs = seq(0, 1, length.out = 12))), labels = FALSE, include.lowest = TRUE)) %>%
    group_by(ev_ind) %>%
    mutate(ev_graph = mean(ev_z, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(monkey, ev_graph, prev_rw_ind) %>% 
    summarise(m = mean(response_mod, na.rm = TRUE) * 100), 
  aes(x = ev_graph, y = m, group = factor(prev_rw_ind), fill = factor(prev_rw_ind))
) +
  stat_summary(fun.data = 'mean_se', geom = 'pointrange', size = 0.75, shape = 21, alpha = ALPHA_LEVEL, color = 'black') + 
  #geom_point(size = 2.25, shape = 21, alpha = ALPHA_LEVEL, color = 'black', position = position_jitter(width = 0.1, height = 0.1)) + 
  stat_summary(fun.data = 'mean_se', geom = 'pointrange', size = 1, shape = 21, alpha = ALPHA_LEVEL, color = 'black') + 
  geom_smooth(method = 'lm', linewidth = LINE_WIDTH, color = 'black') + 
  scale_y_continuous(name = 'Pursue-rate [%]', limits = c(20, 100), breaks = c(20, 40, 60, 80, 100)) + 
  scale_x_continuous(name = 'EV [Z]', breaks = c(-1.0, 0, 1.0), labels = number_format(accuracy = 1)) + 
  scale_fill_manual(name = 'Richness [1-back]', labels = c('Low', 'High'), values = c('skyblue1', 'blue4')) + 
  theme_pubr() + 
  theme(
    text = element_text(size = TEXT_SIZE_MAIN),
    axis.text = element_text(size = TEXT_SIZE_AXIS),
    legend.text = element_text(size = TEXT_SIZE_AXIS),
    aspect.ratio = 1,
    legend.background = element_blank(), 
    legend.position = c(0.50, 0.95),
    legend.direction = 'vertical'
  ) +
  coord_cartesian(xlim=c(-1.5, 1.5))
figs1bi

figure_data <- d %>%
  mutate(
    ave_rw_10back_ind = ifelse(ave_rw_10back_z < 0, -1, ifelse(ave_rw_10back_z >= 0, 1, NA)),
  )

figs1bii <- ggplot(
  figure_data %>%
    drop_na(ave_rw_10back_ind) %>%
    mutate(ev_ind = cut(ev_z, breaks = unique(quantile(ev_z, probs = seq(0, 1, length.out = 12))), labels = FALSE, include.lowest = TRUE)) %>%
    group_by(ev_ind) %>%
    mutate(ev_graph = mean(ev_z, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(monkey, ev_graph, ave_rw_10back_ind) %>% 
    summarise(m = mean(response_mod, na.rm = TRUE) * 100), 
  aes(x = ev_graph, y = m, group = factor(ave_rw_10back_ind), fill = factor(ave_rw_10back_ind))
) +
  stat_summary(fun.data = 'mean_se', geom = 'pointrange', size = 0.75, shape = 21, alpha = ALPHA_LEVEL, color = 'black') + 
  #geom_point(size = 2.25, shape = 21, alpha = ALPHA_LEVEL, color = 'black', position = position_jitter(width = 0.1, height = 0.1)) + 
  stat_summary(fun.data = 'mean_se', geom = 'pointrange', size = 1, shape = 21, alpha = ALPHA_LEVEL, color = 'black') + 
  geom_smooth(method = 'lm', linewidth = LINE_WIDTH, color = 'black') + 
  scale_y_continuous(name = 'Pursue-rate [%]', limits = c(20, 100), breaks = c(20, 40, 60, 80, 100)) + 
  scale_x_continuous(name = 'EV [Z]', breaks = c(-1.0, 0, 1.0), labels = number_format(accuracy = 1)) + 
  scale_fill_manual(name = 'Richness [10-back]', labels = c('Low', 'High'), values = c('skyblue1', 'blue4')) + 
  theme_pubr() + 
  theme(
    text = element_text(size = TEXT_SIZE_MAIN),
    axis.text = element_text(size = TEXT_SIZE_AXIS),
    legend.text = element_text(size = TEXT_SIZE_AXIS),
    aspect.ratio = 1,
    legend.background = element_blank(), 
    legend.position = c(0.55, 0.95),
    legend.direction = 'vertical'
  ) +
  coord_cartesian(xlim=c(-1.5, 1.5))
figs1bii

# supplementary figures S1C&D: effect of richness and stochasticity of the environment
figs1c <-  ggplot(d %>%
                   drop_na(ave_rw_z) %>%
                   mutate(ave_rw_ind = cut(ave_rw_z, breaks = unique(quantile(ave_rw_z, probs = seq(0, 1, length.out = 11))), labels = FALSE, include.lowest = TRUE)) %>%
                   group_by(ave_rw_ind) %>%
                   mutate(ave_rw_graph = mean(ave_rw_z, na.rm = T)) %>%
                   ungroup() %>%
                   group_by(monkey, ave_rw_graph) %>% summarise(m = mean(response_mod, na.rm = T)*100), 
                 aes(x = ave_rw_graph, y = m)) + 
  stat_summary(fun.data = 'mean_se', geom = 'pointrange', size = 0.75, shape = 21, alpha = ALPHA_LEVEL, color = 'black', fill = 'blue4') + 
  geom_smooth(method = 'lm', linewidth = LINE_WIDTH, color = 'black', fill = 'blue4') + 
  scale_y_continuous(name = 'Pursue-rate [%]', limits = c(0, 100), breaks = seq(0, 100, 20)) + 
  scale_x_continuous(name = 'Env. richness [Z]', breaks = seq(-2, 3,1), labels = number_format(accuracy = 0.1)) + 
  theme_pubr() + 
  theme(
    text = element_text(size = 18),
    axis.text = element_text(size = 14),
    aspect.ratio = 1,
    legend.background = element_blank()
  )   + 
  coord_cartesian(ylim=c(22, 100)) + 
  border('black')
figs1c

figs1d <-  ggplot(d %>%
                   drop_na(sd_ave_rw_z) %>%
                   mutate(sd_ave_rw_ind = cut(sd_ave_rw_z, breaks = unique(quantile(sd_ave_rw_z, probs = seq(0, 1, length.out = 11))), labels = FALSE, include.lowest = TRUE)) %>%
                   group_by(sd_ave_rw_ind) %>%
                   mutate(sd_ave_rw_graph = mean(sd_ave_rw_z, na.rm = T)) %>%
                   ungroup() %>%
                   group_by(monkey, sd_ave_rw_graph) %>% summarise(m = mean(response_mod, na.rm = T)*100), 
                 aes(x = sd_ave_rw_graph, y = m)) + 
  stat_summary(fun.data = 'mean_se', geom = 'pointrange', size = 0.75, shape = 21, alpha = ALPHA_LEVEL, color = 'black', fill = '#6b0f2e') + 
  geom_smooth(method = 'lm', linewidth = LINE_WIDTH, color = 'black', fill = '#6b0f2e') + 
  scale_y_continuous(name = 'Pursue-rate [%]', limits = c(0, 100), breaks = seq(0, 100, 20)) + 
  scale_x_continuous(name = 'Env. stochasticity [Z]', breaks = seq(-2, 3,1), labels = number_format(accuracy = 0.1)) + 
  theme_pubr() + 
  theme(
    text = element_text(size = 18),
    axis.text = element_text(size = 14),
    aspect.ratio = 1,
    legend.background = element_blank(),
    plot.background = element_blank(),
  )   + 
  coord_cartesian(ylim=c(20, 100)) + 
  border('black')
figs1d
