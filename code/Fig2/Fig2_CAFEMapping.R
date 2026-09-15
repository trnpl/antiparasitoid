### Figure 2 — CAFE tree mapping ###
### Last modified by RLT 260830. Claude Sonnet 5 assisted with converting CAFE data to mappable states. ###

## Load packages

library(ape)
library(dplyr)
library(ggplot2)
library(ggtree)

## Read trees 

ppo1 <- read.tree('Gamma_asr-PPO1.tre')
ppo234 <- read.tree('Gamma_asr-PPO234.tre')

## Get CAFE values from each node

ppo1   <- ppo1$`TREEPPO1=`
ppo234 <- ppo234$`TREEPPO234=`

## Get CAFE values from each node
ppo_nodedata <- data.frame(
  ppo1label = c(ppo1$node.label, ppo1$tip.label),
  ppo2label = c(ppo234$node.label, ppo234$tip.label)
) %>%
  mutate(ppo1_copy = as.numeric(sapply(strsplit(ppo1label, ">_"), function(x) x[2])))
ppo_nodedata <- ppo_nodedata %>% 
  mutate(ppo2_copy = as.numeric(sapply(strsplit(ppo2label, ">_"), function(x) x[2])))
ppo_nodedata <- ppo_nodedata %>%
  mutate(node = as.integer(gsub(".*<([^>]+)>.*", "\\1", ppo1label))) %>%  # renamed label -> node
  select(node, ppo1_copy, ppo2_copy)

# Map state data to branches 

ppo2_colors <- c("0" = "#DCDCDC", "1" = "#A2A2A2", "2" = "#4E4E4E", "3" = "#000000", "4" = "#FFD070")

# PPO1 data
cafetest <- revts(ggtree(ppo1)) %<+% ppo_nodedata + 
  geom_tree(aes(color = as.factor(ppo1_copy), size = ppo1_copy)) +
  scale_size_continuous(range = c(0.25, 2), breaks = c(0,1,2,3,4)) +
  scale_color_manual(values = ppo2_colors) +
  theme_tree2(legend.position = 'none') 

# PPO2 data

cafetest2 <- revts(ggtree(ppo234)) %<+% ppo_nodedata + 
  geom_tree(aes(color = as.factor(ppo2_copy), size = ppo2_copy)) +
  scale_size_continuous(range = c(0.25, 2), breaks = c(0,1,2,3,4)) +
  scale_color_manual(values = ppo2_colors) +
  theme_tree2(legend.position = 'none') 

# Find nodes to collapse 

getMRCA(newtree, c('DROSOPHILA_HYDEI', 'DROSOPHILA_POLYCHAETA'))
#447

getMRCA(newtree, c('DROSOPHILA_SETOSIMENTUM', 'DROSOPHILA_PRIMAEVA'))
#511

getMRCA(newtree, c('HIRTODROSOPHILA_TRIVITTATA', 'MULGRAVEA_PARASIATICA'))
#580

getMRCA(newtree, c('DROSOPHILA_ANANASSAE', 'DROSOPHILA_SETIFEMUR'))
#744

getMRCA(newtree, c('DROSOPHILA_HISTRIO', 'DROSOPHILA_ORNATIFRONS'))
#594

getMRCA(newtree, c('DROSOPHILA_IMMIGRANS', 'DROSOPHILA_PRUINOSA'))
#628 

getMRCA(newtree, c('ZAPRIONUS_NIGRANUS', 'SPHAEROGASTRELLA_JAVANA'))
#645

getMRCA(newtree, c('DROSOPHILA_TRISTIS', 'DROSOPHILA_GUANCHE'))
#761

# Collapse nodes and plot collapsed tree

cafe2collapsed <- scaleClade(cafetest2, 447, 0.2) %>% 
  scaleClade(511, 0.2) %>%
  scaleClade(580, 0.2) %>%
  scaleClade(594, 0.2) %>%
  scaleClade(628, 0.2) %>% 
  scaleClade(645, 0.2) %>% 
  scaleClade(761, 0.2) %>% 
  scaleClade(674, 0.2) %>% 
  scaleClade(744, 0.2) %>% 
  #scaleClade(800, 0.2) %>% 
  #scaleClade(810, 0.2) %>% 
  ggtree::collapse(447, 'max', color = "#B3B3B3", fill= "#D9D9D9") %>%
  ggtree::collapse(511, 'max', color = "#B3B3B3", fill= "#D9D9D9") %>%
  ggtree::collapse(580, 'max', color = "#B3B3B3", fill= "#D9D9D9") %>%
  ggtree::collapse(594, 'max', color = "#B3B3B3", fill= "#D9D9D9") %>%
  ggtree::collapse(628, 'max', color = "#B3B3B3", fill= "#D9D9D9") %>% 
  ggtree::collapse(645, 'max', color = "#B3B3B3", fill= "#D9D9D9") %>% 
  ggtree::collapse(744, 'max', color = "#B3B3B3", fill= "#D9D9D9") %>%
  ggtree::collapse(761, 'max', color = "#B3B3B3", fill= "#D9D9D9") %>%
  #ggtree::collapse(800, 'max', color = "#B3B3B3", fill= "#D9D9D9") %>%
  #ggtree::collapse(810, 'max', color = "#B3B3B3", fill= "#D9D9D9") %>%
  ggtree::collapse(674, 'max', color = "#B3B3B3", fill= "#D9D9D9") 

ggsave('cafetestcollapsed_new.pdf', cafe2collapsed, height = 20, width = 10, dpi = 300)
