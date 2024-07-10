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
behaviour_path <- "/Volumes/LaCie/nature_neuro_submission/data/behaviour/tus_experiment/tus_behaviour.csv"
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
                    mutate(plot_name = 'DRN-vs-Sham') %>%
                    mutate(ave_rw_ind = cut(ave_rw_z, breaks = unique(quantile(ave_rw_z, probs = seq(0, 1, length.out = 100))), labels = FALSE, include.lowest = TRUE)) %>%
                    group_by(ave_rw_ind) %>%
                    mutate(ave_rw_graph = mean(ave_rw_z, na.rm = T)) %>%
                    ungroup() %>%
                    filter(stim_cond == 'DRN'|stim_cond == 'Sham') %>%
                    group_by(plot_name, monkey, stim_cond, ave_rw_graph) %>% summarise(m = mean(response_mod, na.rm = T)*100), 
                  aes(x = ave_rw_graph, y = m, group = stim_cond, fill = stim_cond, color = stim_cond, linetype = stim_cond)) + 
  geom_smooth(method = 'lm', linewidth = 0.65, alpha = 0.25) + 
  stat_summary(data = d %>%
                 drop_na(ave_rw_z) %>%
                 mutate(ave_rw_ind = cut(ave_rw_z, breaks = unique(quantile(ave_rw_z, probs = seq(0, 1, length.out = 11))), labels = FALSE, include.lowest = TRUE)) %>%
                 group_by(ave_rw_ind) %>%
                 mutate(ave_rw_graph = mean(ave_rw_z, na.rm = T)) %>%
                 ungroup() %>%
                 filter(stim_cond == 'DRN'|stim_cond == 'Sham') %>%
                 group_by(stim_cond, monkey, ave_rw_graph) %>% 
                 group_by(session_id, stim_cond, ave_rw_graph) %>% summarise(m = mean(response_mod, na.rm = T)*100),
    fun.data = 'mean_se', geom = 'pointrange', shape = 21, size = 0.75, alpha = 1, show.legend=FALSE, color = 'black', linetype = 1) + 
  scale_y_continuous(name = 'Pursue-rate [%]', limits = c(0, 100), breaks = seq(0, 100, 20)) + 
  scale_x_continuous(name = 'Env. richness [Z]', breaks = seq(-2, 3,1)) + 
  scale_fill_manual(name = 'TUS condition', values = c('red3', 'blue')) + 
  scale_color_manual(name = 'TUS condition', values = c('red3', 'navy')) + 
  scale_linetype_manual(name = 'TUS condition', values = c(2, 1)) + 
  theme_pubr() + 
  theme(text = element_text(size = 16), axis.text = element_text(size = 14), legend.title = element_blank(), legend.background = element_blank(), legend.position = c(0.3, 0.9), aspect.ratio = 1)   + 
  coord_cartesian(ylim=c(30, 85)) + 
  border('black', size = 1.1) + 
  facet_grid(~plot_name)
fig3d_i

fig3d_ii <- ggplot(d %>%
                    drop_na(ave_rw_z) %>%
                    mutate(plot_name = 'DRN-vs-STS') %>%
                    mutate(ave_rw_ind = cut(ave_rw_z, breaks = unique(quantile(ave_rw_z, probs = seq(0, 1, length.out = 100))), labels = FALSE, include.lowest = TRUE)) %>%
                    group_by(ave_rw_ind) %>%
                    mutate(ave_rw_graph = mean(ave_rw_z, na.rm = T)) %>%
                    ungroup() %>%
                    filter(stim_cond == 'DRN'|stim_cond == 'STS') %>%
                    group_by(plot_name, monkey, stim_cond, ave_rw_graph) %>% summarise(m = mean(response_mod, na.rm = T)*100), 
                  aes(x = ave_rw_graph, y = m, group = stim_cond, fill = stim_cond, color = stim_cond, linetype = stim_cond)) + 
  geom_smooth(method = 'lm', linewidth = 0.65, alpha = 0.25) + 
  stat_summary(data = d %>%
                 drop_na(ave_rw_z) %>%
                 mutate(ave_rw_ind = cut(ave_rw_z, breaks = unique(quantile(ave_rw_z, probs = seq(0, 1, length.out = 11))), labels = FALSE, include.lowest = TRUE)) %>%
                 group_by(ave_rw_ind) %>%
                 mutate(ave_rw_graph = mean(ave_rw_z, na.rm = T)) %>%
                 ungroup() %>%
                 filter(stim_cond == 'DRN'|stim_cond == 'STS') %>%
                 group_by(stim_cond, monkey, ave_rw_graph) %>% 
                 group_by(session_id, stim_cond, ave_rw_graph) %>% summarise(m = mean(response_mod, na.rm = T)*100),
               fun.data = 'mean_se', geom = 'pointrange', shape = 21, size = 0.75, alpha = 1, show.legend=FALSE, color = 'black', linetype = 1) + 
  scale_y_continuous(name = 'Pursue-rate [%]', limits = c(0, 100), breaks = seq(0, 100, 20)) + 
  scale_x_continuous(name = 'Env. richness [Z]', breaks = seq(-2, 3,1)) + 
  scale_fill_manual(name = 'TUS condition', values = c('red3', 'blue')) + 
  scale_color_manual(name = 'TUS condition', values = c('red3', 'navy')) + 
  scale_linetype_manual(name = 'TUS condition', values = c(2, 1)) + 
  theme_pubr() + 
  theme(text = element_text(size = 16), axis.text = element_text(size = 14), legend.title = element_blank(), legend.background = element_blank(), legend.position = c(0.3, 0.9), aspect.ratio = 1)   + 
  coord_cartesian(ylim=c(30, 85)) + 
  border('black', size = 1.1) + 
  facet_grid(~plot_name)
fig3d_ii

fig3d_iii <- ggplot(d %>%
                     drop_na(ave_rw_z) %>%
                     mutate(plot_name = 'DRN-vs-VTA') %>%
                     mutate(ave_rw_ind = cut(ave_rw_z, breaks = unique(quantile(ave_rw_z, probs = seq(0, 1, length.out = 100))), labels = FALSE, include.lowest = TRUE)) %>%
                     group_by(ave_rw_ind) %>%
                     mutate(ave_rw_graph = mean(ave_rw_z, na.rm = T)) %>%
                     ungroup() %>%
                     filter(stim_cond == 'DRN'|stim_cond == 'VTA') %>%
                     mutate(stim_cond = factor(stim_cond, levels = c('VTA', 'DRN'))) %>%
                     group_by(plot_name, monkey, stim_cond, ave_rw_graph) %>% summarise(m = mean(response_mod, na.rm = T)*100), 
                   aes(x = ave_rw_graph, y = m, group = stim_cond, fill = stim_cond, color = stim_cond, linetype = stim_cond)) + 
  geom_smooth(method = 'lm', linewidth = 0.65, alpha = 0.25) + 
  stat_summary(data = d %>%
                 drop_na(ave_rw_z) %>%
                 mutate(ave_rw_ind = cut(ave_rw_z, breaks = unique(quantile(ave_rw_z, probs = seq(0, 1, length.out = 11))), labels = FALSE, include.lowest = TRUE)) %>%
                 group_by(ave_rw_ind) %>%
                 mutate(ave_rw_graph = mean(ave_rw_z, na.rm = T)) %>%
                 ungroup() %>%
                 filter(stim_cond == 'DRN'|stim_cond == 'VTA') %>%
                 group_by(stim_cond, monkey, ave_rw_graph) %>% 
                 group_by(session_id, stim_cond, ave_rw_graph) %>% summarise(m = mean(response_mod, na.rm = T)*100),
               fun.data = 'mean_se', geom = 'pointrange', shape = 21, size = 0.75, alpha = 1, show.legend=FALSE, color = 'black', linetype = 1) + 
  scale_y_continuous(name = 'Pursue-rate [%]', limits = c(0, 100), breaks = seq(0, 100, 20)) + 
  scale_x_continuous(name = 'Env. richness [Z]', breaks = seq(-2, 3,1)) + 
  scale_fill_manual(name = 'TUS condition', values = c('red3', 'blue')) + 
  scale_color_manual(name = 'TUS condition', values = c('red3', 'navy')) + 
  scale_linetype_manual(name = 'TUS condition', values = c(2, 1)) + 
  theme_pubr() + 
  theme(text = element_text(size = 16), axis.text = element_text(size = 14), legend.title = element_blank(), legend.background = element_blank(), legend.position = c(0.3, 0.9), aspect.ratio = 1)   + 
  coord_cartesian(ylim=c(30, 85)) + 
  border('black', size = 1.1) + 
  facet_grid(~plot_name)
fig3d_iii


fig3d_iii <- ggplot(d %>%
                     drop_na(ave_rw_z) %>%
                     mutate(ave_rw_ind = cut(ave_rw_z, breaks = unique(quantile(ave_rw_z, probs = seq(0, 1, length.out = 100))), labels = FALSE, include.lowest = TRUE)) %>%
                     group_by(ave_rw_ind) %>%
                     mutate(ave_rw_graph = mean(ave_rw_z, na.rm = T)) %>%
                     ungroup() %>%
                     filter(stim_cond == 'DRN'|stim_cond == 'VTA') %>%
                     group_by(stim_cond, monkey, ave_rw_graph) %>% 
                     group_by(monkey, stim_cond, ave_rw_graph) %>% summarise(m = mean(response_mod, na.rm = T)*100), 
                   aes(x = ave_rw_graph, y = m, group = stim_cond, fill = stim_cond, color = stim_cond, linetype = stim_cond)) + 
  stat_summary(data = d %>%
                 drop_na(ave_rw_z) %>%
                 mutate(ave_rw_ind = cut(ave_rw_z, breaks = unique(quantile(ave_rw_z, probs = seq(0, 1, length.out = 11))), labels = FALSE, include.lowest = TRUE)) %>%
                 group_by(ave_rw_ind) %>%
                 mutate(ave_rw_graph = mean(ave_rw_z, na.rm = T)) %>%
                 ungroup() %>%
                 filter(stim_cond == 'DRN'|stim_cond == 'VTA') %>%
                 group_by(stim_cond, monkey, ave_rw_graph) %>% 
                 group_by(session_id, stim_cond, ave_rw_graph) %>% summarise(m = mean(response_mod, na.rm = T)*100),
               fun.data = 'mean_se', geom = 'pointrange', shape = 21, color = 'black', size = 0.75, alpha = 0.75) + 
  geom_smooth(method = 'lm', linewidth = 0.65, alpha = 0.35) + 
  scale_y_continuous(name = 'Pursue-rate [%]', limits = c(0, 100), breaks = seq(0, 100, 20)) + 
  scale_x_continuous(name = 'Env. richness [Z]', breaks = seq(-2, 3,1)) + 
  scale_fill_manual(name = 'TUS condition', values = c('cadetblue1', '#0b2370')) + 
  scale_color_manual(name = 'TUS condition', values = c('grey5', 'grey5')) + 
  theme_pubr() + 
  theme(text = element_text(size = 16), axis.text = element_text(size = 14), legend.title = element_blank(), legend.background = element_blank(), legend.position = c(0.3, 0.9), aspect.ratio = 1)   + 
  coord_cartesian(ylim=c(30, 85)) + 
  border('black', size = 1.1)
fig3d_iii

ggplot(d %>%
         drop_na(ave_rw_z) %>%
         mutate(ave_rw_bin = ifelse(ave_rw_z > 0, 'high', 'low')) %>%
         group_by(monkey, session_id, stim_cond, ave_rw_bin, mag) %>%
         summarise(m = mean(response_mod, na.rm = T)), aes(x = stim_cond, y = m, group = ave_rw_bin, fill = ave_rw_bin)) + 
  geom_point(shape = 21, alpha = 0.50, color = 'black', position = position_dodge(width = 1)) + 
  stat_summary(fun = 'mean', geom = 'bar', position = position_dodge(width = 1), alpha = 0.50) +
  theme_pubr() + 
  theme(text = element_text(size = 16), axis.text = element_text(size = 14), legend.title = element_blank(), legend.background = element_blank(), legend.position = c(0.3, 0.9), aspect.ratio = 1.5) 
  
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

animal_data <- d %>%
  drop_na(prev_response_mod) %>%
  group_by(session_id, stim_cond, prev_response_mod) %>%
  summarise(m = mean(response_mod, na.rm = T)*100)
# Reshape data to wide format
animal_data_wide <- animal_data %>%
  pivot_wider(names_from = prev_response_mod, values_from = m, names_prefix = "prev_")

# Calculate the difference score
animal_data_diff <- animal_data_wide %>%
  mutate(diff_score = prev_1 - prev_0) %>% 
  mutate(stim_cond = factor(stim_cond, levels = c('Sham', 'STS', 'VTA', 'DRN')))# Adjust column names based on the actual levels of prev_response_mod

# View the resulting table with the difference score
print(animal_data_diff)

white_to_red_palette <- colorRampPalette(c('white', 'darkorchid4', "darkred"))(4)
fig3e <- ggplot(animal_data_diff, aes(x = factor(stim_cond), fill = factor(stim_cond), y = diff_score)
) + 
  stat_summary(fun = 'mean', geom = 'bar', position = position_dodge2(width = 0.8), color = 'black', width = 0.80, alpha = 0.75) + 
  geom_point(data = animal_data_diff, aes(x = factor(stim_cond), fill = factor(stim_cond), y = diff_score),
             size = 2.5, colour = 'grey50', shape = 21, alpha = 0.35, position = position_dodge(width = 0.8)) + 
  geom_hline(yintercept = 0, color = 'red', linetype = 2) + 
  scale_y_continuous(name = TeX("$\\Delta$ Pursue-rate [%]"), breaks = seq(20, 100, 20)) + 
  scale_x_discrete(name = 'TUS-condition') + 
  scale_fill_manual(name = 'TUS-condition', values = white_to_red_palette) + 
  theme_pubr() + 
  theme(text = element_text(size = 16), axis.text = element_text(size = 14), legend.background = element_blank(), legend.position = 'none', aspect.ratio = 1) + 
  border(color = 'black')
fig3e

figure_data <- d %>%
  drop_na(prev_response_mod)
fig3e <- ggplot(figure_data, aes(x = stim_cond, y = response_mod, group = prev_response_mod)) + 
  stat_summary(fun = 'mean', geom = 'bar', color = 'black', position = 'dodge')
fig3e

fig3e <- ggplot(animal_data, aes(x = stim_cond, y = m, group = prev_response_mod)) + 
  stat_summary(fun = 'mean', geom = 'bar', color = 'black', position = 'dodge')
fig3e

