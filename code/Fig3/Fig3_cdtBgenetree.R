## Mapping characters to CdtB gene tree
## Last modified by RLT 260914 ##

# Load packages
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggnewscale)
library(ape)
library(phytools)
library(ggtree)
library(ggtreeExtra)
library(ggpp)
library(phangorn)
library(MetBrewer)
library(treeio)

## Load data 
labels <- read.csv('cdtb_tiplabels.csv')
tree <- ape::read.nexus("fixed.nexus")
tree$node.label 

# Fix tip labels to match 
tips <- tree$tip.label
tree$tip.label <- gsub("'", "", tree$tip.label)
tree$node.label
tree$tip.label[tree$tip.label == "XCA34573"] <- "XCA34573.1_MAG__hypothetical_protein_ABS861_04190__Wolbachia_endosymbiont_of_Oeneis_ivallda_"
mismatches <- tips[!(tips %in% labels$label)]

## Collapse nodes with <50 bootstrap support

# Labels are "alrt/bootstrap" e.g. "93.9/100"
# Split on "/" and take the second value
bootstrap <- sapply(strsplit(tree$node.label, "/"), function(x) as.numeric(x[2]))
bootstrap

# Identify which internal nodes have bootstrap < 70
low_support <- which(bootstrap < 70) + Ntip(tree)

# Set those edge lengths to 0
for (node in low_support) {
  edge_idx <- which(tree$edge[,2] == node)
  tree$edge.length[edge_idx] <- 0
}

# Collapse zero-length branches into polytomies
tree_collapsed <- di2multi(tree, tol = 1e-8)

# Verify
plot(tree_collapsed, cex = 0.2)
tree_collapsed$tip.label <- gsub("'", "", tree_collapsed$tip.label)

## Reorder synteny labels 
labels$locus <- factor(labels$locus, levels = c("1", "1*", "2", "3", "4", "5", "6", "7", "8", "9"))

## Merge tree data with tip state label 

# First, save bootstrap values to map to nodes

node_data <- data.frame(
  node = (Ntip(tree_collapsed) + 1):(Ntip(tree_collapsed) + Nnode(tree_collapsed)),
  label = tree_collapsed$node.label
) %>%
  mutate(bootstrap = as.numeric(sapply(strsplit(label, "/"), function(x) x[2])))

# Plot tree
p <- ggtree(tree_collapsed, layout = 'rectangular')  %<+% node_data

plab <- p + geom_tiplab()
p$data <- left_join(p$data, labels, by = "label")

## Collapse nodes on the cartoon tree 

# First, need to find the MRCA nodes for each clade

# Drosophila subgenus node: 
getMRCA(tree_collapsed, tip = c("_HIRTODROSOPHILA_TRIVITTATA.nanopore_contig_1026_1_-_CdtB_translation", "_ZAPRIONUS_NIGRANUS.GCA_018903425_JAEIGD010000012._extraction_Copy_-_putative_cdtb_translation"))
# 613

# ananassae-like node: 
getMRCA(tree_collapsed, tip = c("_RHINOLEUCOPHENGA_AMERICANA_JBNMFV010001541_-_putative_cdtb_translation", "_DROSOPHILA_ATRIPEX_cdtB_translation"))
# 568

# microbial node: 
getMRCA(tree_collapsed, tip = c("ENG9823380.1_hypothetical_protein__Campylobacter_coli_", "HDL6963311.1_cytolethal_distending_toxin_subunit_B_family_protein__Yersinia_enterocolitica_"))
# 355

# yersinia node
getMRCA(tree_collapsed, tip = c("'HDL6963311.1_cytolethal_distending_toxin_subunit_B_family_protein__Yersinia_enterocolitica_'", "'RDH39918.1_MAG__cytolethal_distending_toxin_subunit_B_family_protein__Candidatus_Aquirickettsiella_gammari_'"))
#671

## Plot characters onto tree

# Plot character data 
p2 <- p + geom_rootedge(rootedge = 0.1) +
  geom_fruit(
    geom = geom_point,
    mapping = aes(y = label, size = exon_count),
    offset = 0.05
  ) +
  scale_size_continuous(breaks = c(1, 2, 3, 5), range = c(0.5,1.5)) +
  new_scale_fill() + 
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, fill = factor(taxonomy)),
    width = 0.05
  ) +
  scale_fill_discrete() +
  scale_fill_manual(values = c("Actinobacteria" = "grey95", "Proteobacteria" = "grey70", 
                    "Rhodothermota" = "grey60", "Firmicutes" = "grey75", "Bacteriodetes" = "grey87",
                    "Chlorobiota" = "grey82", "Archaea" = "grey50", "Endosymbiont" = "#FFD06FFF",
                    "Insect" = "#EF8A47FF", "Drosophilidae" = "#E76254FF")) + 
  new_scale_fill() + 
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, fill = locus),
    width = 0.05) +
  scale_fill_manual(values = met.brewer("Hiroshige", 24, direction = -1), na.translate = FALSE, drop = FALSE)

p3 <- p2 %>% collapse(node = 355, mode = "max", color = "black", fill = "grey70") %>%
  collapse(node = 613, mode = "max", color = "black", fill = "#1E466E") %>%
  collapse(node = 568, mode = "max", color = "black", fill = "#305E8B") 

p4 <- p3 + new_scale("size") +
  geom_nodepoint(aes(subset = !isTip & as.numeric(bootstrap), size = as.numeric(bootstrap), col = as.numeric(bootstrap)), shape = 21, fill = "white") + 
  scale_size_continuous(range = c(0.1, 1.5)) + 
  scale_color_continuous(low = "white", high = "black") +
  geom_treescale()

ggsave('260611_mappedtree_dots2.svg', p4, width = 8, height = 18)

## full size tree

newtree <- tree_collapsed
newtree$edge.length <- newtree$edge.length * 0.25

pfull <- ggtree(newtree, layout = "rectangular")  %<+% node_data
pfull$data <- left_join(pfull$data, labels, by = "label")

p2full <- p +
  geom_rootedge(rootedge = 0.1) +
  geom_tiplab(align = TRUE, offset = 0.15, size = 1.5) +
  geom_fruit(
    geom = geom_point,
    mapping = aes(y = label, size = exon_count),
    offset = 1
  ) +
  scale_size_continuous(breaks = c(1, 2, 3, 5), range = c(0.5,1.5)) +
  new_scale_fill() + 
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, fill = factor(taxonomy)),
    width = 0.05
  ) +
  scale_fill_discrete() +
  scale_fill_manual(values = c("Actinobacteria" = "grey95", "Proteobacteria" = "grey70", 
                               "Rhodothermota" = "grey60", "Firmicutes" = "grey75", "Bacteriodetes" = "grey87",
                               "Chlorobiota" = "grey82", "Archaea" = "grey50", "Endosymbiont" = "#FFD06FFF",
                               "Insect" = "#EF8A47FF", "Drosophilidae" = "#E76254FF")) + 
  new_scale_fill() + 
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, fill = locus),
    width = 0.05) +
  scale_fill_manual(values = met.brewer("Hiroshige", 24, direction = -1), na.translate = FALSE, drop = FALSE) +
  new_scale("size") +
  geom_nodepoint(aes(subset = !isTip & as.numeric(bootstrap), size = as.numeric(bootstrap), col = as.numeric(bootstrap)), shape = 21, fill = "white") + 
  scale_size_continuous(range = c(0.2, 2.5)) + 
  scale_color_continuous(low = "white", high = "black") +
  geom_treescale()
  
  ggsave('testfull_new.pdf', p2full, width = 10, height = 40)
  
# Highlight insect CdtB nodes (for SI figure)
  
getMRCA(tree_collapsed, tip = c("_LEUCOPHENGA_MACULATA.nanopore_contig_3768_1_-_putative_cdtb_translation", "_DROSOPHILA_ATRIPEX_cdtB_translation"))
# 550
  
dros <- viewClade(p2full, node = 550)

ggsave('drosinset.pdf', dros, width = 20, height = 25)
