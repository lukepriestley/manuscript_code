rm(list = ls())
setwd("/Volumes/LaCie/nature_neuro_submission/data/fmri/tc_output")
library(ggplot2)
library(tidyverse)
library(gghalves)
library(ggdist)
library(rstatix)
library(ggpubr)
library(latex2exp)
library(dplyr)
library(scales)

###################### Prepare TC data ##########################
# List CSV files matching timecourse pattern
file_paths <- dir(pattern = '*decision.csv')

# Read and process CSV files using lapply
d_list <- lapply(file_paths, function(path) {
  # Read the CSV file
  df <- read.csv(path, header = T)
  
  # Add a 'time' column
  df$timepoint <- 1:nrow(df)
  
  # Return the processed data frame
  return(df)
})

d <- bind_rows(d_list) %>%
  rename(
    m = tc_data1,
    se = tc_data2,
    region = tc_data3,
    glm = tc_data4,
    contrast = tc_data5,
  ) %>%
  mutate(
    region = dplyr::recode(region,
                    "DRN_decision" = "DRN",
                    "HB_LP_decision" = "Hb",
                    "LC_decision" = "LC",
                    "NB_decision" = "NB",
                    "BF_LP_decision" = "BF",
                    "SN_decision" = "SN",
                    "VTA_decision" = "VTA",
                    "SMA_decision" = "SMA",
                    "ACC_func_decision" = "ACC",
                    "AI_func_decision" = "AI",
                    "MRN_decision" = "MRN",
                    "vent_decision" = "Ventricle"),
    timepoint = rescale(timepoint, to = c(-1, 5))
  )

####################### FIG 3 ######################

# Figure 3A
# Filter the data
filtered_data <- d %>%
  filter(
    contrast == 'contr_3',
    glm %in% c('GLM_02', 'GLM_03', 'GLM_04'),
    region %in% c('DRN', 'Hb', 'VTA', 'SN', 'NB', 'LC')
  ) %>%
  mutate(fig_panel = ifelse(glm %in% c('GLM_03', 'GLM_04'), 1, 0))

# Create the base ggplot
fig3a <- ggplot(filtered_data, aes(
  x = timepoint, y = m, ymin = m - se, ymax = m + se,
  group = glm, fill = glm, color = glm
)) +
  geom_hline(yintercept = 0, linetype = 1, alpha = 1, color = 'grey80') +
  geom_vline(xintercept = 0, linetype = 1, alpha = 1, color = 'grey80') +
  geom_ribbon(alpha = 0.50, color = NA) +
  geom_line(linewidth = 0.5, color = 'black') +
  labs(
    x = 'Time (s) [0=decision]',
    y = TeX("$\\beta$-richness [a.u.]")
  ) +
  scale_y_continuous(limits = c(-0.085, 0.095), breaks = seq(-0.05, 0.05, 0.05)) +
  scale_fill_manual(
    name = 'Behavioural history',
    values = c('cadetblue2', '#B2D4FF', '#0b2370'),
    labels = c('All', 'Pursue', 'Reject')
  ) +
  theme_pubr() +
  border('black') +
  facet_grid(fig_panel ~ region) +
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.position = 'none',
    legend.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    aspect.ratio = 1.2,
    strip.text.y = element_blank(),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines")
  ) + 
  guides(fill = guide_legend(
    override.aes=list(size =0.25, alpha = 0.8)))
fig3a

# Figure 3C
# Filter the data
filtered_data <- d %>%
  filter(
    region %in% c('VTA', 'DRN', 'Hb'),
    contrast %in% c('contr_1', 'contr_2'),
    glm == 'GLM_05'
  )

# Create the ggplot
fig_3c <- ggplot(filtered_data, aes(
  x = timepoint, y = m, ymin = m - se, ymax = m + se, fill = contrast
)) +
  geom_hline(yintercept = 0, linetype = 1, alpha = 0.1) +
  geom_vline(xintercept = 0, linetype = 1, alpha = 0.1) +
  geom_ribbon(aes(group = contrast, fill = contrast), alpha = 0.60) +
  geom_line(aes(group = contrast), color = 'black', linewidth = 0.65) +
  labs(
    x = 'Time(s) [0 = decision]',
    y = TeX("$\\beta$-motivation [a.u.]")
  ) +
  scale_y_continuous(limits = c(-0.05, 0.065), breaks = c(-0.05, 0.00, 0.05)) +
  scale_fill_manual(
    values = c('darkred', 'purple4'),
    name = 'Motivation-state effect',
    labels = c('Transition', 'Level')
  ) +
  theme_pubr() +
  facet_grid(~region) + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.position = 'top',
    legend.background = element_blank(),
    legend.text = element_text(size = 14),
    aspect.ratio = 1.2,
    strip.text.y = element_blank(),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines")
  ) + 
  guides(fill = guide_legend(
    override.aes=list(size =0.25, alpha = 0.65)))
fig_3c

# Filter the data
filtered_data <- d %>%
  filter(
    region == 'DRN',
    contrast == 'contr_1',
    glm %in% c('GLM_06', 'GLM_07')
  )

# Figure 3E
# Create the ggplot
fig_3e <- ggplot(filtered_data, aes(
  x = timepoint, y = m, ymin = m - se, ymax = m + se, fill = glm
)) +
  geom_hline(yintercept = 0, linetype = 1, alpha = 0.1) +
  geom_vline(xintercept = 0, linetype = 1, alpha = 0.1) +
  geom_ribbon(aes(group = glm, fill = glm), alpha = 0.65) +
  geom_line(aes(group = glm), color = 'black', linewidth = 0.65) +
  labs(
    x = 'Time(s) [0 = decision]',
    y = TeX("$\\beta$-transition [a.u.]")
  ) +
  scale_y_continuous(limits = c(-0.05, 0.05), breaks = c(-0.05, 0.00, 0.05)) +
  scale_fill_manual(
    values = c('darkred', 'pink'),
    name = 'Transition-type',
    labels = c('High-to-Low', 'Low-to-High')
  ) +
  theme_pubr() +
    theme(
      axis.title.y = element_text(size = 18),
      axis.title.x = element_text(size = 18),
      strip.text.x = element_text(size = 14),
      axis.text = element_text(size = 14),
      legend.background = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      aspect.ratio = 1,
      strip.text.y = element_blank(),
      strip.background = element_rect(colour='black', fill='grey99'),
    legend.position = 'none'
  )
fig_3e

# Figure 3F
pink_to_red_palette <- colorRampPalette(c("darkred", "pink"))(3)
# Filter the data
filtered_data <- d %>%
  filter(
    (region == 'DRN'|region == 'MRN'|region == 'Ventricle'),
    contrast == 'contr_1',
    glm %in% c('GLM_05')
  )
fig_3f <- ggplot(filtered_data, aes(x = timepoint, y = m, ymin = m-(se), ymax = m+(se)), fill = region) + 
  geom_hline(yintercept = 0, linetype = 1, alpha = 0.1) +
  geom_vline(xintercept = 0, linetype = 1, alpha = 0.1) +
  geom_ribbon(aes(group = region, fill = region), alpha = 0.65) + 
  geom_line(aes(group = region), color='black', linewidth = 0.65) +
  labs(x = 'Time(s) [0 = decision]', 
       y = TeX("$\\beta$-transition [a.u.]")) + 
  scale_y_continuous(limits = c(-0.05, 0.05), breaks = c(-.05, 0.00, .05)) + 
  scale_fill_manual(values = pink_to_red_palette, name = 'ROI', labels = c('DRN', 'MRN', 'Ventricle')) + 
  theme_pubr() + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    aspect.ratio = 1.2,
    strip.text.y = element_blank(),
    strip.background = element_rect(colour='black', fill='grey99'),
    legend.position = 'none'
  ) + 
  facet_grid(~region)
fig_3f

############## FIG 6 #############
# List CSV files matching PPI pattern
file_paths <- dir(pattern = '*PPI.csv')

# Read and process CSV files using lapply
d_list <- lapply(file_paths, function(path) {
  # Read the CSV file
  df <- read.csv(path, header = T)
  
  # Add a 'time' column
  df$timepoint <- 1:nrow(df)
  
  # Return the processed data frame
  return(df)
})

d_ppi <- bind_rows(d_list) %>%
  rename(
    m = tc_data1,
    se = tc_data2,
    region = tc_data3,
    glm = tc_data4,
    contrast = tc_data5,
  ) %>%
  mutate(
    timepoint = rescale(timepoint, to = c(-1, 5)),
    region = case_when(
      region=='HB_LP_decision' & glm=='GLM_09' ~ 'Hb-AI [Decision]',
      region=='DRN_decision' & glm=='GLM_010' ~ 'DRN-AI [Richness]',
      region=='VTA_decision' & glm=='GLM_010' ~ 'VTA-AI [Richness]',
      region=='DRN_decision' & glm=='GLM_011' ~ 'DRN-AI [Richness]',
      region=='VTA_decision' & glm=='GLM_011' ~ 'VTA-AI [Richness]',
      region=='DRN_decision_to_AI_func_decision' & glm=='GLM_012' ~ 'DRN-AI [Transition]',
      region=='DRN_decision_to_AI_func_decision' & glm=='GLM_013' ~ 'DRN-AI [Transition]',
      region=='VTA_decision_to_AI_func_decision' & glm=='GLM_014' ~ 'VTA-AI [Level]'
    )
  )

# Filter the data
filtered_data <- d_ppi %>%
  filter(
    region == 'Hb-AI [Decision]',
    contrast == 'contr_1',
    glm == 'GLM_09'
  )

# Figure 6A(ii)
# Create the ggplot
fig6a2 <- ggplot(filtered_data, aes(
  x = timepoint, y = m, ymin = m - se, ymax = m + se
)) +
  geom_hline(yintercept = 0, linetype = 1, alpha = 0.1) +
  geom_vline(xintercept = 0, linetype = 1, alpha = 0.1) +
  geom_ribbon(fill = 'azure3', alpha = 0.65) +
  geom_line(color = 'black', linewidth = 0.75) +
  labs(
    x = 'Time(s) [decision = 0]',
    y = TeX("$\\beta$ PPI-Decision [a.u.]")
  ) + 
  scale_y_continuous(limits = c(-0.015, 0.04), breaks = seq(-0.03, 0.03, 0.03)) +
  theme_pubr() +
  border('black') +
  facet_grid(~region) + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.position = 'top',
    legend.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    aspect.ratio = 1,
    strip.text.y = element_text(),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines")
  )
fig6a2

# Figure 6A(i)
# Filter the data
filtered_data <- d_ppi %>%
  filter(
    region == 'DRN-AI [Richness]',
    contrast == 'contr_1',
    glm %in% c('GLM_010', 'GLM_011')
  )

# Create the ggplot
fig6a1_upper <- ggplot(filtered_data, aes(
  x = timepoint, y = m, group = glm, fill = glm, ymin = m - se, ymax = m + se
)) +
  geom_hline(yintercept = 0, linetype = 1, alpha = 0.1) +
  geom_vline(xintercept = 0, linetype = 1, alpha = 0.1) +
  geom_ribbon(alpha = 0.50) +
  geom_line(color = 'black', linewidth = 0.75) +
  labs(
    x = 'Time [s]',
    y = TeX("$\\beta$ PPI-Richness [a.u.]")
  ) + 
  scale_y_continuous(limits = c(-0.06, 0.06), breaks = seq(-0.05, 0.05, 0.05)) +
  scale_fill_manual(
    name = 'Behavioural history',
    labels = c('Pursue', 'Reject'),
    values = c('#B2D4FF', '#0b2370')
  ) +
  theme_pubr() +
  border('black') +
  facet_grid(~region) + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.position = 'none',
    legend.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    aspect.ratio = 1,
    strip.text.y = element_text(),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines"))
fig6a1_upper

# Filter the data
filtered_data <- d_ppi %>%
  filter(
    region == 'DRN-AI [Transition]',
    contrast == 'contr_1',
    glm %in% c('GLM_012', 'GLM_013')
  )

# Create the ggplot
fig6a1_lower <- ggplot(filtered_data, aes(
  x = timepoint, y = m, group = glm, fill = glm, ymin = m - se, ymax = m + se
)) +
  geom_hline(yintercept = 0, linetype = 1, alpha = 0.1) +
  geom_vline(xintercept = 0, linetype = 1, alpha = 0.1) +
  geom_ribbon(alpha = 0.50) +
  geom_line(color = 'black', linewidth = 0.75) +
  labs(
    x = 'Time(s) [0 = decision]',
    y = TeX("$\\beta$ PPI-transition [a.u.]")
  ) + 
  scale_y_continuous(limits = c(-0.035, 0.06), breaks = seq(-0.05, 0.05, 0.05)) +
  scale_fill_manual(
    name = 'Motivation-state transition',
    labels = c('High-to-Low', 'Low-to-High'),
    values = c('pink', 'darkred')
  ) +
  theme_pubr() +
  border('black') +
  facet_grid(~region) + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.position = c(0.40, 0.88),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    aspect.ratio = 1,
    strip.text.y = element_text(),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines"))
fig6a1_lower 

# Figure 6A(iii)
# Filter the data
filtered_data <- d_ppi %>%
  filter(
    region == 'VTA-AI [Richness]',
    contrast == 'contr_1',
    glm %in% c('GLM_010', 'GLM_011')
  )

# Create the ggplot
fig6a3_upper <- ggplot(filtered_data, aes(
  x = timepoint, y = m, group = glm, fill = glm, ymin = m - se, ymax = m + se
)) +
  geom_hline(yintercept = 0, linetype = 1, alpha = 0.1) +
  geom_vline(xintercept = 0, linetype = 1, alpha = 0.1) +
  geom_ribbon(alpha = 0.50) +
  geom_line(color = 'black', linewidth = 0.75) +
  labs(
    x = 'Time [s]',
    y = TeX("$\\beta$ PPI-Richness [a.u.]")
  ) + 
  scale_y_continuous(limits = c(-0.04, 0.06), breaks = seq(-0.05, 0.05, 0.05)) +
  scale_fill_manual(
    name = 'Behavioural history',
    labels = c('Pursue', 'Reject'),
    values = c('#B2D4FF', '#0b2370')
  ) +
  theme_pubr() +
  border('black') +
  facet_grid(~region) + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.position = 'none',
    legend.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    aspect.ratio = 1,
    strip.text.y = element_text(),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines"))
fig6a3_upper

# Filter the data
filtered_data <- d_ppi %>%
  filter(
    region == 'VTA-AI [Level]',
    contrast == 'contr_1',
    glm %in% c('GLM_014')
  )

# Create the ggplot
fig6a3_lower <- ggplot(filtered_data, aes(
  x = timepoint, y = m, ymin = m - se, ymax = m + se
)) +
  geom_hline(yintercept = 0, linetype = 1, alpha = 0.1) +
  geom_vline(xintercept = 0, linetype = 1, alpha = 0.1) +
  geom_ribbon(alpha = 0.65, fill = 'purple4') +
  geom_line(color = 'black', linewidth = 0.75) +
  labs(
    x = 'Time(s) [0 = decision]',
    y = TeX("$\\beta$ PPI-level [a.u.]")
  ) + 
  scale_y_continuous(limits = c(-0.0240, 0.04), breaks = seq(0, 0.03, 0.03)) +
  theme_pubr() +
  border('black') +
  facet_grid(~region) + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.position = c(0.40, 0.88),
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    aspect.ratio = 1,
    strip.text.y = element_text(),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines"))
fig6a3_lower

############## FIGURE S5 ###################
# Filter the data
filtered_data <- d %>%
  filter(
    contrast == 'contr_4',
    glm %in% c('GLM_02'),
    region %in% c('DRN', 'Hb', 'VTA', 'SN', 'NB', 'LC')
  )

# Create the base ggplot
figs5a <- ggplot(filtered_data, aes(
  x = timepoint, y = m, ymin = m - se, ymax = m + se,
  group = glm, fill = glm, color = glm
)) +
  geom_hline(yintercept = 0, linetype = 1, alpha = 1, color = 'grey80') +
  geom_vline(xintercept = 0, linetype = 1, alpha = 1, color = 'grey80') +
  geom_ribbon(alpha = 0.50, fill = 'slateblue', color = NA) + 
  geom_line(linewidth = 0.5, color = 'black') +
  labs(
    x = 'Time [s]',
    y = TeX("$\\beta$-stochasticity [a.u.]")
  ) +
  scale_y_continuous(breaks = seq(from = -0.1, to = 0.1, by = 0.05), limits = c(-0.05, 0.05)) + 
  theme_pubr() +
  border('black') +
  facet_grid(~region) +
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.position = 'top',
    legend.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    aspect.ratio = 1,
    strip.text.y = element_blank(),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines")
  ) + 
  guides(fill = guide_legend(
    override.aes=list(size =0.25, alpha = 0.8)))
figs5a

# Filter the data
filtered_data <- d %>%
  filter(
    contrast == 'contr_4',
    glm %in% c('GLM_02'),
    region %in% c('AI', 'ACC', 'SMA')
  )

# Create the base ggplot
figs5b <- ggplot(filtered_data, aes(
  x = timepoint, y = m, ymin = m - se, ymax = m + se,
  group = glm, fill = glm, color = glm
)) +
  geom_hline(yintercept = 0, linetype = 1, alpha = 1, color = 'grey80') +
  geom_vline(xintercept = 0, linetype = 1, alpha = 1, color = 'grey80') +
  geom_ribbon(alpha = 0.50, fill = 'slateblue', color = NA) + 
  geom_line(linewidth = 0.5, color = 'black') +
  labs(
    x = 'Time [s]',
    y = TeX("$\\beta$-stochasticity [a.u.]")
  ) +
  scale_y_continuous(breaks = seq(from = -0.1, to = 0.1, by = 0.05), limits = c(-0.05, 0.075)) + 
  theme_pubr() +
  border('black') +
  facet_grid(~region) +
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.position = 'top',
    legend.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    aspect.ratio = 1,
    strip.text.y = element_blank(),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines")
  ) + 
  guides(fill = guide_legend(
    override.aes=list(size =0.25, alpha = 0.8)))
figs5b

############## FIGURE S6 ###################

# Filter the data
filtered_data <- d %>%
  filter(
    contrast == 'contr_1',
    glm %in% c('GLM_01'),
    region %in% c('DRN', 'Hb', 'VTA', 'SN', 'NB', 'LC')
  )

# Create the base ggplot
figs6a <- ggplot(filtered_data, aes(
  x = timepoint, y = m, ymin = m - se, ymax = m + se,
  group = glm, fill = glm, color = glm
)) +
  geom_hline(yintercept = 0, linetype = 1, alpha = 1, color = 'grey80') +
  geom_vline(xintercept = 0, linetype = 1, alpha = 1, color = 'grey80') +
  geom_ribbon(alpha = 0.50, fill = 'seagreen3', color = NA) + 
  geom_line(linewidth = 0.5, color = 'black') +
  labs(
    x = 'Time [s]',
    y = TeX("$\\beta$-decision [a.u.]")
  ) +
  scale_y_continuous(breaks = seq(from = -0.1, to = 0.1, by = 0.05), limits = c(-0.05, 0.05)) + 
  theme_pubr() +
  border('black') +
  facet_grid(~region) +
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.position = 'top',
    legend.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    aspect.ratio = 1,
    strip.text.y = element_blank(),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines")
  ) + 
  guides(fill = guide_legend(
    override.aes=list(size =0.25, alpha = 0.8)))
figs6a

# Filter the data
filtered_data <- d %>%
  filter(
    contrast == 'contr_1',
    glm %in% c('GLM_01'),
    region %in% c('AI', 'ACC', 'SMA')
  )

# Create the base ggplot
figs6c <- ggplot(filtered_data, aes(
  x = timepoint, y = m, ymin = m - se, ymax = m + se,
  group = glm, fill = glm, color = glm
)) +
  geom_hline(yintercept = 0, linetype = 1, alpha = 1, color = 'grey80') +
  geom_vline(xintercept = 0, linetype = 1, alpha = 1, color = 'grey80') +
  geom_ribbon(alpha = 0.50, fill = 'seagreen3', color = NA) + 
  geom_line(linewidth = 0.5, color = 'black') +
  labs(
    x = 'Time [s]',
    y = TeX("$\\beta$-decision [a.u.]")
  ) +
  scale_y_continuous(breaks = seq(from = -0.1, to = 0.1, by = 0.05), limits = c(-0.05, 0.065)) + 
  theme_pubr() +
  border('black') +
  facet_grid(~region) +
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.position = 'top',
    legend.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    aspect.ratio = 1,
    strip.text.y = element_blank(),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines")
  ) + 
  guides(fill = guide_legend(
    override.aes=list(size =0.25, alpha = 0.8)))
figs6c

############## FIGURE S7 ###################

# Filter the data
filtered_data <- d %>%
  filter(
    contrast == 'contr_3',
    glm %in% c('GLM_02', 'GLM_03', 'GLM_04'),
    region %in% c('ACC', 'AI', 'SMA')
  ) %>%
  mutate(fig_panel = ifelse(glm %in% c('GLM_03', 'GLM_04'), 1, 0))

# Create the base ggplot
figs7a <- ggplot(filtered_data, aes(
  x = timepoint, y = m, ymin = m - se, ymax = m + se,
  group = glm, fill = glm, color = glm
)) +
  geom_hline(yintercept = 0, linetype = 1, alpha = 1, color = 'grey80') +
  geom_vline(xintercept = 0, linetype = 1, alpha = 1, color = 'grey80') +
  geom_ribbon(alpha = 0.50, color = NA) +
  geom_line(linewidth = 0.5, color = 'black') +
  labs(
    x = 'Time [s]',
    y = TeX("$\\beta$-richness [a.u.]")
  ) +
  scale_y_continuous(limits = c(-0.085, 0.095), breaks = seq(-0.05, 0.05, 0.05)) +
  scale_fill_manual(
    name = 'Behavioural history',
    values = c('cadetblue2', '#B2D4FF', '#0b2370'),
    labels = c('All', 'Pursue', 'Reject')
  ) +
  theme_pubr() +
  border('black') +
  facet_grid(fig_panel ~ region) +
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.position = 'none',
    legend.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    aspect.ratio = 1,
    strip.text.y = element_blank(),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines")
  ) + 
  guides(fill = guide_legend(
    override.aes=list(size =0.25, alpha = 0.8)))
figs7a

############## FIGURE S8 ###################

# Filter the data
filtered_data <- d %>%
  filter(
    contrast == 'contr_2',
    glm == c('GLM_08'),
    region == ('DRN')
  )

# Create the base ggplot
figs8a <- ggplot(filtered_data, aes(
  x = timepoint, y = m, ymin = m - se, ymax = m + se,
  group = glm, fill = glm, color = glm
)) +
  geom_hline(yintercept = 0, linetype = 1, alpha = 1, color = 'grey80') +
  geom_vline(xintercept = 0, linetype = 1, alpha = 1, color = 'grey80') +
  geom_ribbon(alpha = 0.60, color = NA, fill = 'orangered2') +
  geom_line(linewidth = 0.5, color = 'black') +
  labs(
    x = 'Time [s]',
    y = TeX("$\\beta$-Ave.EV [a.u.]")
  ) +
  scale_y_continuous(breaks = seq(-0.05, 0.05, 0.05), limits = c(-0.05, 0.075)) +
  theme_pubr() +
  border('black') +
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.position = 'none',
    legend.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    aspect.ratio = 1,
    strip.text.y = element_blank(),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines")
  ) + 
  facet_grid(~region) + 
  guides(fill = guide_legend(
    override.aes=list(size =0.25, alpha = 0.8)))
figs8a

# Filter the data
filtered_data <- d %>%
  filter(
    contrast == 'contr_1',
    glm == c('GLM_05'),
    region %in% c('DRN', 'MRN', 'Ventricle')
  )

fig_s8c <- ggplot(data = d %>%
                    filter(
                      region == 'DRN',
                      contrast == 'contr_1',
                      glm %in% c(
                        'GLM_050', 'GLM_051', 'GLM_052',
                        'GLM_053', 'GLM_054', 'GLM_055', 'GLM_056'
                      )
                    ), aes(x = timepoint, y = m, ymin = m - se, ymax = m + se, fill = glm)) + 
  geom_hline(yintercept = 0, linetype = 1, alpha = 0.1) +
  geom_vline(xintercept = 0, linetype = 1, alpha = 0.1) +
  geom_ribbon(aes(group = glm, fill = glm), alpha = 0.50) + 
  geom_line(aes(group = glm, color = glm), linewidth = 0.65) +
  labs(
    x = 'Time [s]',
    y = TeX("$\\beta$-transition [a.u.]")
  ) +
  scale_y_continuous(limits = c(-0.05, 0.025)) + 
  scale_fill_brewer(
    name = 'Time-window',
    palette = 'Reds',
    labels = c('(t-6):(t+0)', '(t-5):(t+1)', '(t-4):(t+2)', 
               '(t-3):(t+3)', '(t-2):(t+4)', '(t-1):(t+5)', '(t-0):(t+6)')
  ) + 
  scale_color_brewer(
    name = 'Time-window',
    palette = 'Reds',
    labels = c('(t-6):(t+0)', '(t-5):(t+1)', '(t-4):(t+2)', 
               '(t-3):(t+3)', '(t-2):(t+4)', '(t-1):(t+5)', '(t-0):(t+6)')
  ) + 
  theme_pubr() + 
  theme(
    text = element_text(size = 18),
    legend.background = element_blank(),
    legend.position = 'right',
    aspect.ratio = 1
  ) + 
  border(color = 'black', size = 1) + 
  guides(fill = 'none', color = 'none')
fig_s8c

############### PREPARE PEAK DATA #############

# List CSV files matching the pattern
list_files <- dir(pattern = '*peaks.csv')

# Read and process CSV files using lapply
d_peak_list <- lapply(list_files, function(file) {
  df <- read.csv(file)
  return(df)
})

# Combine data frames into one using dplyr
d_peak <- bind_rows(d_peak_list) %>%
  pivot_longer(cols = ends_with('decision'), names_to = 'region', values_to = 'peak') %>%
  mutate(region = dplyr::recode(region,
                         DRN_decision = 'DRN',
                         HB_LP_decision = 'Hb',
                         LC_decision = 'LC',
                         NB_decision = 'NB',
                         BF_LP_decision = 'BF',
                         SN_decision = 'SN',
                         VTA_decision = 'VTA',
                         SMA_decision = 'SMA',
                         ACC_func_decision = 'ACC',
                         AI_func_decision = 'AI',
                         MRN_decision = 'MRN',
                         vent_decision = 'Ventricle'
  ))


###################### FIGURE 3 ###################

fig3b <- ggplot(data = d_peak %>%
                     filter(
                       contrast == 'contr3',
                       GLM %in% c('GLM_02', 'GLM_03', 'GLM_04'),
                       region %in% c('DRN', 'Hb', 'VTA')
                     ), aes(x = region, y = peak, fill = GLM)) + 
  stat_summary(
    fun = 'mean',
    geom = 'bar', 
    color = 'black', 
    width = 0.8, 
    position = position_dodge(width = 0.8, preserve = 'total')
  ) + 
  geom_point(
    shape = 21, 
    color = 'grey50',
    alpha = 0.25,
    position = position_dodge(width = 0.8, preserve = 'total'),
    size = 2
  ) + 
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, color = 'red') + 
  scale_y_continuous(
    breaks = seq(-1, 1, 0.25), limits = c(-0.5, 0.5),
    labels = scales::number_format(accuracy = 0.01)
  ) + 
  scale_fill_manual(
    name = 'Behavioural history',
    values = c('cadetblue2', '#B2D4FF', '#0b2370'),
    labels = c('All', 'Pursue', 'Reject')
  ) + 
  labs(
    x = 'ROI', 
    y = TeX("$\\beta$-richness [a.u.]")
  ) + 
  theme_pubr() + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.position = 'none',
    legend.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    aspect.ratio = 0.50,
    strip.text.y = element_blank(),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines")
  )
fig3b

fig3d <- ggplot(data = d_peak %>%
                  filter(
                    contrast %in% c('contr1', 'contr2'),
                    GLM == 'GLM_05',
                    region %in% c('DRN', 'Hb', 'VTA')
                  ), aes(x = region, y = peak, fill = contrast)) + 
  stat_summary(
    fun = 'mean',
    geom = 'bar', 
    color = 'black', 
    width = 0.8, 
    alpha = 0.85,
    position = position_dodge(width = 0.8, preserve = 'total')
  ) + 
  geom_point(
    shape = 21, 
    color = 'grey50',
    alpha = 0.25,
    position = position_dodge(width = 0.8, preserve = 'total'),
    size = 2
  ) + 
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, color = 'red') + 
  scale_y_continuous(
    breaks = seq(-1, 1, 0.25),
    labels = scales::number_format(accuracy = 0.01)
  ) + 
  scale_fill_manual(
    name = NULL,
    values = c('darkred', 'purple4'),
    labels = c('Transition', 'Level')
  ) + 
  labs(
    x = 'ROI', 
    y = TeX("$\\beta$-motivation [a.u.]")
  ) + 
  theme_pubr() + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.position = 'none', 
    legend.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    aspect.ratio = 0.75,
    strip.text.y = element_blank(),
    #strip.background = element_rect(colour='black', fill='grey95'),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines")
  )
fig3d

t.test(d_peak$peak[d_peak$contrast=='contr1' & d_peak$GLM=='GLM_05' & d_peak$region=='DRN'],
       d_peak$peak[d_peak$contrast=='contr1' & d_peak$GLM=='GLM_05' & d_peak$region=='VTA'], paired = T)

t.test(d_peak$peak[d_peak$contrast=='contr2' & d_peak$GLM=='GLM_05' & d_peak$region=='DRN'],
       d_peak$peak[d_peak$contrast=='contr2' & d_peak$GLM=='GLM_05' & d_peak$region=='VTA'], paired = T)

fig3e_inset <- ggplot(data = d_peak %>%
                                filter(
                                  region == 'DRN',
                                  contrast == 'contr1',
                                  GLM %in% c('GLM_06', 'GLM_07')
                                ), aes(x = GLM, y = peak, fill = GLM)) + 
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, color = 'red') + 
  stat_summary(
    fun = 'mean',
    geom = 'bar', 
    color = 'black', 
    width = 0.5, 
    alpha = 0.75,
    position = position_dodge2(width = 0)
  ) +   
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', color = 'black', width = 0.20, linewidth = 1) + 
  scale_y_continuous(breaks = c(-0.04, -0.02, 0.00)) + 
  scale_x_discrete(labels = NULL) + 
  scale_fill_manual(values = c('darkred', 'pink')) + 
  labs(x = NULL, y = NULL) + 
  theme_pubr() + 
  theme(
    axis.text = element_text(size = 35),
    axis.ticks.y = element_blank(),
    legend.position = 'none',
    legend.background = element_blank(),
    plot.background = element_blank(),
    aspect.ratio = 0.50
  ) + 
  border(size = 2.5) +
  coord_flip(ylim = c(-0.041, 0.00))
fig3e_inset

pink_to_red_palette <- colorRampPalette(c("darkred", "pink", "white"))(3)
fig3g <- ggplot(data = d_peak %>%
                  filter(
                    contrast %in% c('contr1'),
                    GLM == 'GLM_05',
                    region %in% c('DRN', 'MRN', 'Ventricle')
                  ), aes(x = region, y = peak, fill = region)) + 
  stat_summary(
    fun = 'mean',
    geom = 'bar', 
    color = 'black', 
    width = 0.8, 
    position = position_dodge(width = 0.8, preserve = 'total')
  ) + 
  geom_point(
    shape = 21, 
    color = 'grey50',
    alpha = 0.25,
    position = position_dodge(width = 0.8, preserve = 'total'),
    size = 2
  ) + 
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, color = 'red') + 
  scale_y_continuous(
    breaks = seq(-0.2, 0.2, 0.2),
    labels = scales::number_format(accuracy = 0.01)
  ) + 
  scale_x_discrete(
    labels = c('DRN', 'MRN', 'Vent.')
    ) + 
  scale_fill_manual(
    name = NULL,
    values = pink_to_red_palette,
    labels = c('DRN', 'MRN', 'Vent.')
  ) + 
  labs(
    x = 'ROI', 
    y = TeX("$\\beta$-transition [a.u.]")
  ) + 
  theme_pubr() + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.position = 'none', 
    legend.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.direction = 'horizontal',
    aspect.ratio = 0.8,
    strip.text.y = element_blank(),
    #strip.background = element_rect(colour='black', fill='grey95'),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines")
  )
fig3g

t.test(d_peak$peak[d_peak$contrast=='contr1' & d_peak$GLM=='GLM_05' & d_peak$region=='DRN'],
       d_peak$peak[d_peak$contrast=='contr1' & d_peak$GLM=='GLM_05' & d_peak$region=='MRN'], paired = T)
t.test(d_peak$peak[d_peak$contrast=='contr1' & d_peak$GLM=='GLM_05' & d_peak$region=='DRN'],
       d_peak$peak[d_peak$contrast=='contr1' & d_peak$GLM=='GLM_05' & d_peak$region=='Ventricle'], paired = T)

x_label <- c(
  '(t-6):(t+0)',
  '(t-5):(t+1)',
  '(t-4):(t+2)',
  '(t-3):(t+3)',
  '(t-2):(t+4)',
  '(t-1):(t+5)',
  '(t-0):(t+6)')
pink_to_red_palette <- colorRampPalette(c("lightpink", "darkred"))(7)

fig3g <- ggplot(data = d_peak %>%
                  filter(contrast == 'contr1' & region == 'DRN' & GLM %in% c('GLM_050', 'GLM_051', 'GLM_052', 'GLM_053', 'GLM_054', 'GLM_055', 'GLM_056')), aes(x = GLM, y = peak, fill = GLM)) + 
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, color = 'red') + 
  stat_summary(fun = mean, geom = 'line', group = 1) + 
  stat_summary(fun.data = mean_ci, geom = 'pointrange', shape = 21, size = 1) + 
  scale_x_discrete(name = 'Time-window', labels = x_label) + 
  scale_fill_manual(name = 'Time-window', values = pink_to_red_palette, labels = x_label) + 
  scale_color_manual(name = 'Time-window', values = pink_to_red_palette, labels = x_label) + 
  labs(x = 'ROI', 
       y = TeX("$\\beta$-transition [a.u.]")) + 
  theme_pubr() + 
  theme(    axis.title.y = element_text(size = 18),
            axis.title.x = element_text(size = 18),
            strip.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14),
            axis.text.x = element_text(size = 14, angle = 35, vjust = 1, hjust=1),
            legend.background = element_blank(),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 14),
            aspect.ratio = 0.50,
            legend.position = 'none')
fig3g

###################### FIGURE S4 ###################

figs4b <- ggplot(data = d_peak %>%
                  filter(
                    contrast == 'contr1',
                    GLM == 'GLM_01',
                    region %in% c('Hb', 'AI')
                  ), aes(x = region, y = peak, fill = contrast)) + 
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, color = 'red') + 
  stat_summary(
    fun = 'mean',
    fill = 'seagreen3',
    geom = 'bar', 
    color = 'black', 
    width = 0.8, 
    position = position_dodge(width = 0.8, preserve = 'total')
  ) + 
  geom_point(
    shape = 21, 
    fill = 'seagreen1',
    alpha = 0.25,
    position = position_dodge(width = 0.8, preserve = 'total'),
    size = 2
  ) + 
  scale_y_continuous(
    limits = c(-0.20, 0.20),
    breaks = seq(-1, 1, 0.10),
    labels = scales::number_format(accuracy = 0.01)
  ) + 
  labs(
    x = 'ROI', 
    y = TeX("$\\beta$-decision [a.u.]")
  ) + 
  theme_pubr() + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.position = c(0.3, 0.95), 
    legend.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    aspect.ratio = 1.2,
    strip.text.y = element_blank(),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines")
  )
figs4b

###################### FIGURE S5 ################
figs5b <- ggplot(data = d_peak %>%
                  filter(
                    contrast == 'contr3',
                    GLM %in% c('GLM_02', 'GLM_03', 'GLM_04'),
                    region == c('AI')
                  ), aes(x = region, y = peak, fill = GLM)) + 
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, color = 'red') + 
  stat_summary(
    fun = 'mean',
    geom = 'bar', 
    color = 'black', 
    width = 0.8, 
    position = position_dodge(width = 1, preserve = 'total')
  ) + 
  geom_point(
    shape = 21, 
    alpha = 0.25,
    position = position_dodge(width = 1, preserve = 'total'),
    size = 2
  ) + 
  scale_y_continuous(
    breaks = seq(-1, 1, 0.25), limits = c(-0.5, 0.5),
    labels = scales::number_format(accuracy = 0.01)
  ) + 
  scale_fill_manual(
    name = 'Behavioural history',
    values = c('cadetblue2', '#B2D4FF', '#0b2370'),
    labels = c('All', 'Pursue', 'Reject')
  ) + 
  labs(
    x = 'ROI', 
    y = TeX("$\\beta$-richness [a.u.]")
  ) + 
  theme_pubr() + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_blank(),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.text.x = element_blank(),
    legend.position = 'none',
    legend.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    aspect.ratio = 1.5,
    strip.text.y = element_blank(),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines")
  ) + 
  facet_grid(~region)
figs5b

###################### FIG s6 ###################
figs6b <- ggplot(data = d_peak %>%
                   filter(
                     contrast == 'contr2',
                     GLM == c('GLM_08'),
                     region == c('DRN')
                   ), aes(x = region, y = peak, fill = GLM)) + 
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, color = 'red') + 
  geom_hline(yintercept = 0, size = 1, linetype = 2, alpha = 1, color = 'red') + 
  stat_summary(
    fun = 'mean',
    geom = 'bar', 
    color = 'black', 
    width = 0.5, 
    alpha = 0.85,
    position = position_dodge2(width = 0)
  ) +   
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', color = 'black', width = 0.20, size = 2.25) + 
  scale_y_continuous(
    breaks = seq(-1, 1, 0.02),
    labels = scales::number_format(accuracy = 0.01)
  ) + 
  labs(
    x = 'ROI', 
    y = TeX("$\\beta$-Ave.EV [a.u.]")
  ) + 
  theme_pubr() + 
  border('black', size = 2) + 
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 30),
    axis.text.y = element_blank(),
    legend.position = 'none',
    legend.background = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    aspect.ratio = 0.5,
    strip.text.y = element_blank(),
    strip.background = element_rect(colour='black', fill='grey99'),
    panel.spacing.x = unit(0.2, "lines")
  ) + 
  coord_flip()
figs6b

figs6d <- ggplot(data = d_peak %>%
                            filter(
                              contrast == 'contr1',
                              GLM %in% c(
                                'GLM_050', 'GLM_051', 'GLM_052',
                                'GLM_053', 'GLM_054', 'GLM_055', 'GLM_056'
                              ),
                              region == 'DRN'
                            ), aes(x = GLM, y = peak, fill = GLM)) + 
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, color = 'red') + 
  stat_summary(
    fun = 'mean',
    geom = 'bar', 
    color = 'black', 
    width = 0.8, 
    position = position_dodge(width = 1, preserve = 'total')
  ) + 
  geom_point(
    shape = 21, 
    alpha = 0.25,
    position = position_dodge(width = 1, preserve = 'total'),
    size = 2
  ) + 
  scale_x_discrete(
    name = 'Time-window',
    labels = c(
      '(t-6):(t+0)', '(t-5):(t+1)', '(t-4):(t+2)',
      '(t-3):(t+3)', '(t-2):(t+4)', '(t-1):(t+5)', '(t-0):(t+6)'
    )
  ) + 
  scale_y_continuous(breaks = seq(-0.2, 0.2, 0.2)) + 
  scale_fill_brewer(
    name = 'Time-window',
    palette = 'Reds',
    labels = c(
      '(t-6):(t+0)', '(t-5):(t+1)', '(t-4):(t+2)',
      '(t-3):(t+3)', '(t-2):(t+4)', '(t-1):(t+5)', '(t-0):(t+6)'
    )
  ) + 
  scale_color_brewer(
    name = 'Time-window',
    palette = 'Reds',
    labels = c(
      '(t-6):(t+0)', '(t-5):(t+1)', '(t-4):(t+2)',
      '(t-3):(t+3)', '(t-2):(t+4)', '(t-1):(t+5)', '(t-0):(t+6)'
    )
  ) + 
  labs(
    x = 'ROI', 
    y = TeX("$\\beta$-transition [a.u.]")
  ) + 
  theme_pubr() + 
  theme(
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.position = 'none',
    legend.background = element_blank(),
    plot.background = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 0.7, hjust = 0.6),
    aspect.ratio = 0.5
  )
figs6d

figs6d_inset <- ggplot(data = d_peak %>%
                   filter(
                     contrast == 'contr1',
                     GLM %in% c(
                       'GLM_050', 'GLM_051', 'GLM_052',
                       'GLM_053', 'GLM_054', 'GLM_055', 'GLM_056'
                     ),
                     region == 'DRN'
                   ), aes(x = GLM, y = peak, fill = GLM, group = GLM)) + 
  geom_hline(yintercept = 0, linetype = 2, alpha = 1, color = 'red') + 
  stat_summary(
    aes(group = 1),
    fun = mean,
    geom = 'line',
    position = position_dodge(width = 0.8),
    size = 0.65,
    alpha = 1,
    color = 'black'
  ) +
  stat_summary(
    fun.data = mean_ci,
    geom = 'pointrange',
    position = position_dodge(width = 0.8),
    size = 3,
    shape = 21, 
    alpha = 1,
    color = 'black',
    fatten = 1.5
  ) +
  scale_x_discrete(
    name = 'Time-window',
    labels = c(
      '(t-6):(t+0)', '(t-5):(t+1)', '(t-4):(t+2)',
      '(t-3):(t+3)', '(t-2):(t+4)', '(t-1):(t+5)', '(t-0):(t+6)'
    )
  ) + 
  scale_y_continuous(breaks = c(-0.04, -0.02, 0)) + 
  scale_fill_brewer(
    name = 'Time-window',
    palette = 'Reds',
    labels = c(
      '(t-6):(t+0)', '(t-5):(t+1)', '(t-4):(t+2)',
      '(t-3):(t+3)', '(t-2):(t+4)', '(t-1):(t+5)', '(t-0):(t+6)'
    )
  ) + 
  scale_color_brewer(
    name = 'Time-window',
    palette = 'Reds',
    labels = c(
      '(t-6):(t+0)', '(t-5):(t+1)', '(t-4):(t+2)',
      '(t-3):(t+3)', '(t-2):(t+4)', '(t-1):(t+5)', '(t-0):(t+6)'
    )
  ) + 
  labs(
    x = 'ROI', 
    y = TeX("$\\beta$-transition [a.u.]")
  ) + 
  theme_pubr() + 
  theme(
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    legend.position = 'none',
    aspect.ratio = 0.5
  ) + 
  border(color = 'black', size = 1.5)
figs6d_inset

###################### FIG4C ###################
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(scales)
library(ggplot2)
library(fmsb)

d_rsfMRI <- read.csv('rsfMRI_results.csv', header = T)

DRN_TUS <- d_rsfMRI %>% filter(TUS_condition=='DRN_TUS') %>% select(-TUS_condition)
VTA_TUS <- d_rsfMRI %>% filter(TUS_condition=='VTA_TUS') %>% select(-TUS_condition)

# create DRN-TUS plot
rownames(DRN_TUS) <- DRN_TUS$Seed_area
DRN_TUS <- DRN_TUS %>% select(-Seed_area)
DRN_TUS <- rbind(rep(0, ncol(DRN_TUS)), DRN_TUS)
DRN_TUS <- rbind(rep(max(DRN_TUS), ncol(DRN_TUS)), DRN_TUS)

radarchart(
  DRN_TUS, axistype = 1,
  # Customize the polygon
  plty = c(1,1),
  plwd = 2,
  pcol = c('red', 'purple3'),
  pfcol = scales::alpha(c('red', 'purple3'), 0.3),
  # Customize the grid
  cglcol = "black", cglty = 1, cglwd = 0.5,
  # Customize the axis
  axislabcol = "black", 
  # Variable labels
  vlcex = 1.2, vlabels = colnames(DRN_TUS),
  caxislabels = seq(0, 0.25, 0.05), calcex = 1
)
par(mar = rep(0.5, 4))


# create VTA-TUA plot
rownames(VTA_TUS) <- VTA_TUS$Seed_area
VTA_TUS <- VTA_TUS %>% select(-Seed_area)
VTA_TUS <- rbind(rep(0, ncol(VTA_TUS)), VTA_TUS)
VTA_TUS <- rbind(rep(max(VTA_TUS), ncol(VTA_TUS)), VTA_TUS)

radarchart(
  VTA_TUS, axistype = 1,
  # Customize the polygon
  plty = c(1,1),
  plwd = 2,
  pcol = c('red', 'purple3'),
  pfcol = scales::alpha(c('red', 'purple3'), 0.3),
  # Customize the grid
  cglcol = "black", cglty = 1, cglwd = 0.5,
  # Customize the axis
  axislabcol = "black", 
  # Variable labels
  vlcex = 1.2, vlabels = colnames(DRN_TUS),
  caxislabels = seq(0, 0.25, 0.05), calcex = 1
)
par(mar = rep(0.5, 4))


