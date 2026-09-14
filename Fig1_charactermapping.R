### Code to generate Figure 1 — Decorated Species Tree ###
### Last modified by RLT 260830 ###

## Load packages

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggnewscale)
library(ape)
library(phytools)
library(ggtree)
library(ggtreeExtra)
library(ggpp)

## Load files 
newtree <- read.tree('species_rescaled_4d_secCalib_renamed.tree')
states <- read.csv('alltips_new.csv')

## Plot tree

t <- ggtree(newtree, layout = "circular", branch.length = 'none')
t$data <- left_join(t$data, states, by = "label")

t2 <- t + geom_tiplab(size = 0,
                      align = TRUE,
                      linetype = "dashed",
                      color = "grey",
                      linesize = 0.2) 

# Add tiles to tree 

linetype <- c("intact" = "solid", "inc" = "dotted", 
              "pseud_part" = "solid", "pseud_all" = "solid")
point_type <- c("intact" = NA, "inc" = NA, "pseud_part" = 1, "pseud_all" = 19)


t3t <- t2 +   
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, 
                  color = factor(code),
                  fill = factor(code)),
    offset = 0.01,
    linewidth = 0.5,
    width = 0.5) +
  scale_color_manual(values = setNames(rep("black", nlevels(factor(states$code))), 
                                       levels(factor(states$code))), 
                     na.value = NA) +
  scale_fill_discrete(na.value = NA) +
  new_scale_fill() + 
  new_scale_color() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, 
                  color = as.numeric(PPO1),
                  fill = as.numeric(PPO1),
                  linetype = PPO1_qc),
    width = 0.5,
    linewidth = 0.5,
    offset = 0.016) + 
  scale_linetype_manual(values = linetype) +
  scale_color_continuous(low = "black", high = "black", na.value = NA) +
  scale_fill_continuous(low = "white", high = "grey60", na.value = NA) +
  new_scale_fill() + 
  new_scale_color() +
  geom_fruit(
    geom = geom_point,
    mapping = aes(y = label, 
                  shape = as.factor(PPO1_qc)),
    color = "white",
    fill = "white",
    size = 1,
    offset = 0.0) + 
  scale_shape_manual(values = point_type, na.value = NA) +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, 
                  color = as.numeric(PPO2),
                  fill = as.numeric(PPO2),
                  linetype = PPO2_qc),
    width = 0.5,
    linewidth = 0.5,
    offset = 0.016) + 
  scale_color_continuous(low = "black", high = "black", na.value = NA) +
  scale_fill_continuous(low = "white", high = "black", na.value = NA) +
  new_scale_fill() + 
  new_scale_color() +
  new_scale("shape") +
  geom_fruit(
    geom = geom_point,
    mapping = aes(y = label, 
                  shape = as.factor(PPO2_qc)),
    color = "white",
    fill = "white",
    size = 1,
    offset = 0.0) + 
  scale_shape_manual(values = point_type, na.value = NA) +
  new_scale_fill() + 
  new_scale_color() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, 
                  alpha = as.factor(PPO3),
                  color = as.factor(PPO3),
                  fill = as.numeric(PPO3),
                  linetype = PPO3_qc),
    width = 0.5,
    linewidth = 0.5,
    offset = 0.016) +
  scale_alpha_manual(values = c("0" = 0, "1" = 1, "2" = 1, "3" = 1)) +
  scale_color_manual(values = c("0" = NA, "1" = "black", "2" = "black", "3" = "black"), na.value = NA) +
  scale_fill_continuous(low = "white", high = "#E76355", na.value = NA) + 
  new_scale_fill() + 
  new_scale_color() +
  new_scale("shape") +
  geom_fruit(
    geom = geom_point,
    mapping = aes(y = label, 
                  shape = as.factor(PPO3_qc)),
    color = "white",
    fill = "white",
    size = 1,
    offset = 0.0) + 
  scale_shape_manual(values = point_type, na.value = NA) + 
  new_scale_fill() + 
  new_scale("alpha") +
  new_scale_color() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, 
                  alpha = as.factor(PPO4),
                  color = as.factor(PPO4),
                  fill = as.numeric(PPO4),
                  linetype = PPO4_qc),
    width = 0.5,
    linewidth = 0.5,
    offset = 0.0) + 
  scale_alpha_manual(values = c("0" = 0, "1" = 1)) +
  scale_color_manual(values = c("0" = NA, "1" = "black"), na.value = NA) +
  scale_fill_continuous(low = "white", high = "#FCC166", na.value = NA) +
  new_scale_fill() + 
  new_scale_color() +
  new_scale("shape") +
  geom_fruit(
    geom = geom_point,
    mapping = aes(y = label, 
                  shape = as.factor(PPO4_qc)),
    color = "white",
    fill = "white",
    size = 1,
    offset = 0.0) + 
  scale_shape_manual(values = point_type, na.value = NA) +
  new_scale_fill() + 
  new_scale("alpha") +
  new_scale_color() +
  new_scale("pattern") +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, 
                  alpha = as.numeric(cdtb_count),
                  color = as.factor(cdtb_count),
                  fill = as.numeric(cdtb_count),
                  linetype = cdtb_qc),
    width = 0.5,
    linewidth = 0.5,
    offset = 0.016) + 
  scale_color_manual(values = c("0" = NA, "1" = "black", "2" = "black", "3" = "black", "4" = "black", "7" = "black"), na.value = NA) +
  scale_fill_continuous(low = "#72BCD5", high = "#1E466E", na.value = NA) +
  scale_alpha_continuous(range = c(0, 1),
                         rescaler = function(x, ...) ifelse(x == 0, 0, 1)) +
  new_scale_fill() + 
  new_scale_color() +
  new_scale("shape") +
  geom_fruit(
    geom = geom_point,
    mapping = aes(y = label, 
                  shape = as.factor(cdtb_qc)),
    color = "white",
    fill = "white",
    size = 1,
    offset = 0.0) + 
  scale_shape_manual(values = point_type, na.value = NA) 

## Add clade labels

tlab2 <- t3t +
  geom_cladelab(node = 718, label = "melanogaster subgroup", offset = 3) + 
  geom_cladelab(node = 674, label = "montium group", offset = 3) + 
  geom_cladelab(node = 745, label = "ananassae subgroup", offset = 3) + 
  geom_cladelab(node = 760, label = "obscura group", offset = 3) +
  geom_cladelab(node = 779, label = "saltans group", offset = 3) + 
  geom_cladelab(node = 793, label = "willistoni group", offset = 3) +
  geom_cladelab(node = 486, label = "virilis group", offset = 3) + 
  geom_cladelab(node = 454, label = "repleta group", offset = 3) +
  geom_cladelab(node = 598, label = "quinaria group", offset = 3) + 
  geom_cladelab(node = 612, label = "testacea group", offset = 3) +
  geom_cladelab(node = 619, label = "cardini group", offset = 3) + 
  geom_cladelab(node = 629, label = "immigrans group", offset = 3) +
  geom_cladelab(node = 646, label = "Zaprionus", offset = 3) + 
  geom_cladelab(node = 496, label = "robusta + melanica groups", offset = 3) +
  geom_cladelab(node = 510, label = "Hawaiian Drosophila", offset = 3) + 
  geom_cladelab(node = 566, label = "Scaptomyza", offset = 3) + 
  geom_cladelab(node = 415, label = "Steganinae", offset = 3) + 
  geom_cladelab(node = 809, label = "non-Drosophilina Drosophilinae", offset = 3) +
  geom_cladelab(node = 800, label = "Lordiphosa", offset = 3)

tlab3 <- tlab2 + geom_cladelab(node = 444, label = "subgenus Drosophila", offset = 5) + 
  geom_cladelab(node = 670, label = "subgenus Sophophora", offset = 5)

ggsave('cladelabels.pdf', tlab3, width = 20, height = 40)