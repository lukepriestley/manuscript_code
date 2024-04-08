# Code for all analysis involving GLM4, Figure 4D-E ("Non-invasive perturbation 
# of dorsal raphe nucleus shows that it is causally involved in the influence of 
# a monkey’s environment on its behaviour").

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
library(latex2exp)

################### Read data ################
behaviour_path <- "/Volumes/LaCie/cell_submission/data/behaviour/tus_experiment/tus_behaviour.csv"
d <- read.csv(behaviour_path)
d$stim_cond <- factor(d$stim_cond, levels = c('sham', 'sts', 'drn', 'vta'), labels = c('Sham', 'STS', 'DRN', 'VTA')) # code TUS-condition such that sham is the reference group

################### GLMS ################
GLM4.1 <- glmer(response_mod ~ mag_z*prob_z + ave_rw_z*stim_cond + sd_ave_rw_z + trial_z + (1|monkey),
                data = d, family = 'binomial')
summary(GLM4.1)

GLM4.1a <- glmer(response_mod ~ mag_z*prob_z + ave_rw_z*stim_cond + sd_ave_rw_z + trial_z + (1|monkey),
                data = d %>% filter(stim_cond == 'STS'|stim_cond == 'DRN'), family = 'binomial')
summary(GLM4.1a)

GLM4.1b <- glmer(response_mod ~ mag_z*prob_z + ave_rw_z*stim_cond + sd_ave_rw_z + trial_z + (1|monkey),
                 data = d %>% filter(stim_cond == 'VTA'|stim_cond == 'DRN'), family = 'binomial')
summary(GLM4.1b)

GLM4.2 <- glmer(response_mod ~ mag_z*prob_z + prev_response_mod*stim_cond + trial_z + (1|monkey),
                data = d, family = 'binomial')
summary(GLM4.2)

GLM4.2a <- glmer(response_mod ~ mag_z*prob_z + prev_response_mod*stim_cond + trial_z + (1|monkey),
                data = d %>% filter(stim_cond=='DRN'|stim_cond=='STS'), family = 'binomial')
summary(GLM4.2a)

GLM4.2b <- glmer(response_mod ~ mag_z*prob_z + prev_response_mod*stim_cond + trial_z + (1|monkey),
                 data = d %>% filter(stim_cond=='DRN'|stim_cond=='VTA'), family = 'binomial')
summary(GLM4.2b)

################### FIGURES ################

fig3d_i <- ggplot(d %>%
                    drop_na(ave_rw_z) %>%
                    mutate(ave_rw_ind = cut(ave_rw_z, breaks = unique(quantile(ave_rw_z, probs = seq(0, 1, length.out = 100))), labels = FALSE, include.lowest = TRUE)) %>%
                    group_by(ave_rw_ind) %>%
                    mutate(ave_rw_graph = mean(ave_rw_z, na.rm = T)) %>%
                    ungroup() %>%
                    filter(stim_cond == 'DRN'|stim_cond == 'Sham') %>%
                    group_by(stim_cond, monkey, ave_rw_graph) %>% 
                    group_by(monkey, stim_cond, ave_rw_graph) %>% summarise(m = mean(response_mod, na.rm = T)*100), 
                  aes(x = ave_rw_graph, y = m, group = stim_cond, fill = stim_cond)) + 
  geom_point(size = 1.5, shape = 21, alpha = 1, color = 'black', position = position_jitter(width = 0.1)) + 
  geom_smooth(method = 'lm', linewidth = 0.65, color = 'black') + 
  scale_y_continuous(name = 'Pursue-rate [%]', limits = c(0, 100), breaks = seq(0, 100, 20)) + 
  scale_x_continuous(name = 'Env. richness [Z]', breaks = seq(-2, 3,1)) + 
  scale_fill_manual(name = 'TUS condition', values = c('grey90', 'cadetblue1')) + 
  theme_pubr() + 
  theme(text = element_text(size = 16), axis.text = element_text(size = 14), legend.title = element_blank(), legend.background = element_blank(), legend.position = c(0.2, 0.9), aspect.ratio = 1)   + 
  coord_cartesian(ylim=c(15, 105)) + 
  border('black', size = 1.1)
fig3d_i

fig3d_ii <- ggplot(d %>%
                    drop_na(ave_rw_z) %>%
                    mutate(ave_rw_ind = cut(ave_rw_z, breaks = unique(quantile(ave_rw_z, probs = seq(0, 1, length.out = 100))), labels = FALSE, include.lowest = TRUE)) %>%
                    group_by(ave_rw_ind) %>%
                    mutate(ave_rw_graph = mean(ave_rw_z, na.rm = T)) %>%
                    ungroup() %>%
                    filter(stim_cond == 'STS'|stim_cond == 'DRN') %>%
                    group_by(stim_cond, monkey, ave_rw_graph) %>% 
                    group_by(monkey, stim_cond, ave_rw_graph) %>% summarise(m = mean(response_mod, na.rm = T)*100), 
                  aes(x = ave_rw_graph, y = m, group = stim_cond, fill = stim_cond)) + 
  geom_point(size = 1.5, shape = 21, alpha = 1, color = 'black', position = position_jitter(width = 0.1)) + 
  geom_smooth(method = 'lm', linewidth = 0.65, color = 'black') + 
  scale_y_continuous(name = 'Pursue-rate [%]', limits = c(0, 100), breaks = seq(0, 100, 20)) + 
  scale_x_continuous(name = 'Env. richness [Z]', breaks = seq(-2, 3,1)) + 
  scale_fill_manual(name = 'TUS condition', values = c('blue3', 'cadetblue1')) + 
  theme_pubr() + 
  theme(text = element_text(size = 16), axis.text = element_text(size = 14), legend.title = element_blank(), legend.background = element_blank(), legend.position = c(0.2, 0.9), aspect.ratio = 1)   + 
  coord_cartesian(ylim=c(15, 105)) + 
  border('black', size = 1.1)
fig3d_ii

fig3d_iii <- ggplot(d %>%
                     drop_na(ave_rw_z) %>%
                     mutate(ave_rw_ind = cut(ave_rw_z, breaks = unique(quantile(ave_rw_z, probs = seq(0, 1, length.out = 100))), labels = FALSE, include.lowest = TRUE)) %>%
                     group_by(ave_rw_ind) %>%
                     mutate(ave_rw_graph = mean(ave_rw_z, na.rm = T)) %>%
                     ungroup() %>%
                     filter(stim_cond == 'VTA'|stim_cond == 'DRN') %>%
                     group_by(stim_cond, monkey, ave_rw_graph) %>% 
                     group_by(monkey, stim_cond, ave_rw_graph) %>% summarise(m = mean(response_mod, na.rm = T)*100), 
                   aes(x = ave_rw_graph, y = m, group = stim_cond, fill = stim_cond)) + 
  geom_point(size = 1.5, shape = 21, alpha = 1, color = 'black', position = position_jitter(width = 0.1)) + 
  geom_smooth(method = 'lm', linewidth = 0.65, color = 'black') + 
  scale_y_continuous(name = 'Pursue-rate [%]', limits = c(0, 100), breaks = seq(0, 100, 20)) + 
  scale_x_continuous(name = 'Env. richness [Z]', breaks = seq(-2, 3,1)) + 
  scale_fill_manual(name = 'TUS condition', values = c('cadetblue1', '#0b2370')) + 
  theme_pubr() + 
  theme(text = element_text(size = 16), axis.text = element_text(size = 14), legend.title = element_blank(), legend.background = element_blank(), legend.position = c(0.2, 0.9), aspect.ratio = 1)   + 
  coord_cartesian(ylim=c(15, 105)) + 
  border('black', size = 1.1)
fig3d_iii

fig3e <- ggplot(d %>%
           drop_na(prev_response_mod) %>%
           group_by(monkey, session_id, stim_cond, prev_response_mod) %>% 
           summarise(m = mean(response_mod, na.rm = T)*100),
         aes(x = factor(stim_cond, levels = c('Sham', 'DRN', 'STS', 'VTA')), fill = factor(prev_response_mod), y = m)) + 
  stat_halfeye(position = position_dodge(width = 0.8), adjust = 1, justification = 0, .width = 0, point_colour = NA, alpha = 0.65) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.1, alpha = 0.60, outlier.color = NA) +
  geom_half_point(position = position_dodge(width = 0.8), side = "l", range_scale = .4, shape = 21, size = 1.5, alpha = 1, color = 'black') +
  scale_y_continuous(name = 'Pursue-rate[%]', limits = c(20, 100), breaks = seq(20, 100, 20)) + 
  scale_x_discrete(name = 'TUS-condition', labels = c('Sham', 'DRN', 'STS', 'VTA')) + 
  scale_fill_manual(name = 'Behavioural history', values = c('grey90', 'slateblue1'), labels = c('Reject', 'Pursue')) + 
  theme_pubr() + 
  theme(text = element_text(size = 16), axis.text = element_text(size = 14), legend.background = element_blank(), legend.position = 'none', aspect.ratio = 0.40)
fig3e
