### PPO gene tree figures ###
## Last modified by RLT 260725, Claude Sonnet 5 assisted with some code generation ###

## Load packages

library(ape)
library(ggplot2)
library(ggtree)
library(dplyr)
library(phytools)

## Load data

# PPO nucleotide tree

ppont <- read.tree('260724_PPONT.treefile')

## Reroot PPO tree

getMRCA(ppont, c("Daphnia_magna_PPO", "Homarus_grammarus_PPO", "Metacarcinus_magister_PPO"))
#1680

target_node <- 1680
branch_len <- ppont$edge.length[ppont$edge[,2] == target_node]

ppont_rooted <- reroot(ppont, node.number = target_node, position = branch_len / 2)

ppont_midpoint <- midpoint.root(ppont)

getMRCA(ppont_midpoint, c("Daphnia_magna_PPO", "Homarus_grammarus_PPO"))
ppont_root <- root(ppont_midpoint, node = 1681, resolve.root = FALSE)
plot(ppont_rooted)
plot(ppont_midpoint)

# Labels are "alrt/bootstrap" e.g. "93.9/100"
# Split on "/" and take the second value
bootstrap <- sapply(strsplit(ppont_rooted$node.label, "/"), function(x) as.numeric(x[2]))
bootstrap

# Identify which internal nodes have bootstrap < 70
low_support <- which(bootstrap < 70) + Ntip(ppont_rooted)

# Set those edge lengths to 0
for (node in low_support) {
  edge_idx <- which(ppont_rooted$edge[,2] == node)
  ppont_rooted$edge.length[edge_idx] <- 0
}

# Collapse zero-length branches into polytomies
ppont_collapsed <- di2multi(ppont_rooted, tol = 1e-8)

# Verify
plot(ppont_collapsed, cex = 0.2)
ppont_collapsed$tip.label <- gsub("'", "", ppont_collapsed$tip.label)

# save bootstrap values to map to nodes

node_data <- data.frame(
  node = (Ntip(ppont_collapsed) + 1):(Ntip(ppont_collapsed) + Nnode(ppont_collapsed)),
  label = ppont_collapsed$node.label
) %>%
  mutate(bootstrap = as.numeric(sapply(strsplit(label, "/"), function(x) x[2])))

## Plot tree

p <- ggtree(ppont_collapsed, layout = 'rectangular')  %<+% node_data

# Collapse nodes of major clades -- first need to find them
#crustaceans
getMRCA(ppont_collapsed, c("Daphnia_magna_PPO", "Homarus_grammarus_PPO"))
#1658 
# ppo1
getMRCA(ppont_collapsed, c("Rhinoleucophenga_americana_-_17124.1_-_PPO1_-_g2.t1", "Chymomyza_caudatula_-_42_1_-_PPO1_-_CDS"))
#854 
# ppo4
getMRCA(ppont_collapsed, c("Drosophila_neosaltans_NODE_3450_length_5582_cov_36_-_PPO4", "DROSOPHILA_SP._14030-0761.01_90.1_-_PPO4_-_g2.t1"))
#1390
#ppo3 
getMRCA(ppont_collapsed, c("Drosophila_takahashii_c3R_-_NC_091681_-_PPO3_-_LOC108067041_CDS", "Drosophila_rhopaloa_-_NW_025335034_-_PPO3_-_LOC108048366_CDS"))
#1240 
# mel ppo2
getMRCA(ppont_collapsed, c("Drosophila_eugracilis-_024572605_-_PPO2_-_LOC108104532_CDS",
                           "Drosophila_ficusphila_-_025064328_-_PPO2_-_LOC108090478_CDS"))
#1272
# montium ppo2
getMRCA(ppont_collapsed, c("Drosophila_barbarae_-_JAEIIT010046403.1_-_PPO2_-_CDS",
                           "Drosophila_pectinifera_VNKC01000282.1_-_PPO2_-_CDS"))
#1295
#ana ppo2
getMRCA(ppont_collapsed, c("Drosophila_parabipectinata_-_282.1_-_PPO2_-_CDS",
                           "Drosophila_ironensis_-_8764.1_-_PPO2_-_CDS"))
#1337
#obscura ppo2
getMRCA(ppont_collapsed, c("Drosophila_imaii_NODE_1323_length_9261_cov_9.32687_-_PPO2",
                           "Drosophila_subobscura_chrE_-_048531_-_PPO2_-_LOC117890056_CDS"))
#1352
#saltans ppo2
getMRCA(ppont_collapsed, c("Drosophila_nigrosaltans_NODE_15_length_317025_cov__-_PPO2",
                           "Drosophila_sucinea_-_1.1_-_PPO2_-_CDS"))
#1371
#lordiphosa ppo2
getMRCA(ppont_collapsed, c("Lordiphosa_mommai_5.1_-_PPO2_2_-_CDS_2",
                           "Lordiphosa_magnipectinata_-_30.1_-_PPO2_-_CDS"))
#1604
# subgenus drosophila ppo2
getMRCA(ppont_collapsed, c("Drosophila_maculinotata_-_contig_387_-_PPO2_-_CDS",
                           "Drosophila_innubila_chr2R_-_022995374_-_PPO2_-_PPO2_CDS"))
#1402
# duplicated ppo2s in steganinae
getMRCA(ppont_collapsed, c("GITONA_DISTIGMA_21742.1_-_PPO2_-_PPO2_-_g2.t1",
                           "Colocasiomyia_alocasiae_-_59.1_-_PPO2_-_CDS"))
#1633
# chymomyza ppo2
getMRCA(ppont_collapsed, c("Chymomyza_costata_-_379.1_-_PPO2_2_-_CDS_2",
                           "SCAPTODROSOPHILA_INORNATA_19160.1_-_PPO2_-_PPO2_-_CDS"))
#1609
# amiota and friends
getMRCA(ppont_collapsed, c("Cacoxenus_indagator_-_1222_1_-_PPO2_-_CDS",
                           "Amiota_communis_-__1667_1_-_PPO2_-_CDS"))
#1619 
# leucophenga and friends
getMRCA(ppont_collapsed, c("LEUCOPHENGA_SP.J_90046.1_-_PPO2_2_-_PPO2_-_g5.t1",
                           "GITONA_DISTIGMA_21742.1_-_PPO2_2_g3.t1"))
#1626

#Now we can collapse them
p_collapsed <- p %>% 
  ggtree::collapse(node = 1658, mode = "max", color = "black", fill = "black")  %>% #crustaceans
  ggtree::collapse(node = 854, mode = "max", color = "black", fill = "#DCDCDC") %>% # PPO1
  ggtree::collapse(node = 1390, mode = "max", color = "black", fill = "#F39A4F") %>% # PPO4
  ggtree::collapse(node = 1240, mode = "max", color = "black", fill = "#E76355") %>% # PPO3
  ggtree::collapse(node = 1272, mode = "max", color = "black", fill = "#A2A2A2") %>% #mel ppo2
  ggtree::collapse(node = 1295, mode = "max", color = "black", fill = "#A2A2A2") %>% # montium ppo2
  ggtree::collapse(node = 1337, mode = "max", color = "black", fill = "#A2A2A2") %>% #ana ppo2
  ggtree::collapse(node = 1371, mode = "max", color = "black", fill = "#A2A2A2") %>% # saltans ppo2
  ggtree::collapse(node = 1604, mode = "max", color = "black", fill = "#A2A2A2") %>% # lordiphosa ppo2
  ggtree::collapse(node = 1352, mode = "max", color = "black", fill = "#A2A2A2") %>% # obscura ppo2
  ggtree::collapse(node = 1402, mode = "max", color = "black", fill = "#A2A2A2") %>% #subgenus drosophila ppo2
  ggtree::collapse(node = 1633, mode = "max", color = "black", fill = "#4E4E4E") %>% # duplicated steganinae ppo2
  ggtree::collapse(node = 1609, mode = "max", color = "black", fill = "#4E4E4E") %>% # chymomyza ppo2
  ggtree::collapse(node = 1619, mode = "max", color = "black", fill = "#4E4E4E") %>% # amiota ppo2
  ggtree::collapse(node = 1626, mode = "max", color = "black", fill = "#4E4E4E") # leucophenga ppo2

p2 <- p_collapsed +
  geom_nodepoint(aes(subset = !isTip & as.numeric(bootstrap), size = as.numeric(bootstrap)), shape = 21, fill = "white") + 
  scale_size_continuous(range = c(0.1, 1.5)) + 
  geom_treescale(y = 300, width = 0.2)

ggsave('ppont.pdf', p2, width = 8, height = 18, dpi = 300)
