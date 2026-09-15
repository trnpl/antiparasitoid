## Pancrustacea PPO Copy Number Analyses ##
## Last modified by RLT 260913 ##

## Load packages

library(dplyr)
library(ggplot2)
library(multcompView)

## Load data

all_ppo <- read.csv('PancrustaceaPPO.csv')
states <- read.csv('alltips_new.csv')

## Add PPO3/4 character state data 

ppo34 <- states %>% filter(PPO3 > 0 | PPO4 > 0)
ppo34 <- ppo34[,c("species", "PPO3", "PPO4")]
ppo34 <- ppo34 %>% 
  mutate(PPO = case_when(
    PPO3 > 0 ~ "PPO3",
    PPO4 > 0 ~ "PPO4"
  ))
ppo34 <- ppo34[,c("species","PPO")]

all_ppo <- left_join(all_ppo, ppo34, by = "species")

ppocolor <- c("PPO3" = "#E76355", "PPO4" = "#FCC166", na.value = "#DCDCDC")

## Calculate means 

# means by 
mean_ppo <- all_ppo %>% group_by(alt_taxon) %>% 
  summarize(mean = mean(ppo_count), sd = sd(ppo_count), N = table(alt_taxon))

# Set larger taxonomic groups
all_ppo$groups <- NA
all_ppo <- all_ppo %>%
  mutate(groups = case_when(
    alt_taxon %in% c("Diptera_other", "Chironomidae") ~ "Diptera_other",
    alt_taxon == "clade_1" ~ "clade_1",
    alt_taxon == "clade_2" ~ "clade_2", 
    alt_taxon == "drosophilidae_basal" ~ "drosophilidae_basal", 
    alt_taxon == "Crustacean" ~ "Crustacean",
    TRUE ~ "Other_insect"
  ))

all_ppo <- all_ppo %>%
  mutate(orders = case_when(
    alt_taxon %in% c("Diptera_other", "clade_1", "clade_2", "drosophilidae_basal") ~ "Diptera",
    TRUE ~ all_ppo$alt_taxon
  ))

# Calculate mean by taxonomic group
mean_ppo4 <- all_ppo %>% group_by(groups) %>% 
  summarize(mean = mean(ppo_count), sd = sd(ppo_count), N = table(groups))

# Relevel groups variable for plotting
all_ppo$groups <- factor(all_ppo$groups,
                         levels = c("clade_1", "clade_2", "drosophilidae_basal", "Diptera_other", "Other_insect", "Crustacean"))

# Plot PPOs by groups
ppo_plot <- ggplot(all_ppo, aes(x = groups, y = ppo_count), na.omit = TRUE) + geom_jitter(aes(color = PPO), width = 0.2) +
  stat_summary(fun = "mean", geom = "crossbar", width = 0.4) +
  geom_hline(yintercept = 4.41, linetype = "dashed") + 
  ylim(0, 18) +
  scale_color_manual(values = ppocolor) + 
  theme_bw() +
  theme(legend.position = 'none') 

ggsave('ppo_insects.pdf', ppo_plot, width = 6, height = 4, dpi = 300)

## Statistical testing
# Run ANOVA
test <- aov(ppo_count ~ groups, data = all_ppo)

# Tukey Test
tukey_output <- TukeyHSD(test)

# Set significance groups
significance_groups <- multcompLetters4(test, tukey_output)
