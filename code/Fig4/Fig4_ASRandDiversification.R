### Code for Figure 4: cdtB ancestral state reconstruction + diversification rate analyses ### 
### Claude Sonnet 5 assisted in writing code for these analyses. ###
### Last modified by RLT 260830 ###

## Load packages

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
library(ggpattern)
library(ggimage)
library(deeptime)
library(geiger)
library(caper)

## Load data 

newtree <- read.tree('species_rescaled_4d_secCalib_renamed.tree')
states <- read.csv('alltips_new.csv')

## ASR 

# Define states 
states$cdtb_locus[is.na(states$cdtb_locus)] <- "0"
states$cdtb_locus[states$cdtb_locus == "2,4"] <- "2"

# Extract tip states
tip_states <- setNames(as.character(states$cdtb_locus), states$label)

# --- 3. Convert tip states to phyDat object ---

tip_states_mat <- matrix(tip_states, nrow = length(tip_states), ncol = 1,
                          dimnames = list(names(tip_states), "site1"))


tip_states_phy <- phyDat(tip_states_mat, type = "USER", levels = 
                            c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))

tip_states_phy <- as.character(tip_states_phy)
# For fitMk, you need to create an INDEX matrix, not a rate matrix
# The index matrix tells fitMk which transitions share the same rate parameter. We are fitting a model
# such that each independent cdtB acquisition has its own gain/loss rate parameter.

n_locations <- 1
n_states <- 13

# Set values for index matrix
Q_index <- matrix(0, n_states, n_states)

# Gains (absent -> any location): all get unique rate parameter
Q_index[1, 2:n_states] <- 1:(n_states - 1)

# Losses (any location -> absent): all get unique rate parameter
Q_index[2:n_states, 1] <- 13:24

## Edit the cost matrix such that 0->3 transitions cannot occur since we have priors that this wasn't an indepedent acquisition
Q_index[1, 7] <- 0
# Make new rate for 1 -> 3 transition
Q_index[2,7] <- 25 

state_labels <- as.character(0:12)  # or whatever matches your data
rownames(Q_index) <- colnames(Q_index) <- state_labels

# Vectorize tip states
tip_states_vec <- as.vector(tip_states_phy)
names(tip_states_vec) <- rownames(tip_states_phy)

# Fit index matrix. pi sets root prior to 100% probability of no cdtB at root.
fit <- fitMk(newtree, tip_states_vec, model = Q_index, pi = c(1,0,0,0,0,0,0,0,0,0,0,0,0))

# Stochastic character mapping using this fit model
Q_matrix <- as.Qmatrix(fit)
simmap <- make.simmap(newtree, tip_states_vec, model= fit, Q = Q_matrix, pi = c(1,0,0,0,0,0,0,0,0,0,0,0,0), nsim=500)

# Count state changes from simulation data 
change_counts <- countSimmap(simmap)

# Aggregate gain/loss events for each state
change_counts <- as.data.frame(change_counts) 
change_counts <- change_counts %>%
  mutate(gains = rowSums(across(starts_with("Tr.0.")), na.rm = TRUE))
change_counts <- change_counts %>%
  mutate(losses = rowSums(across(ends_with(".0")), na.rm = TRUE))
change_counts_min <- change_counts %>% select(Tr.N, gains, losses)
write.csv(change_counts_min, 'ASRpi0_mixedratesimulations.csv')

# Plot histogram of total gains and losses across all 500 simulations
gain_loss <- ggplot(change_counts_min) + 
  geom_histogram(aes(x = gains, y = after_stat(density)), binwidth = 1, col = "black", fill = "#79C0D6") +
  geom_histogram(aes(x = losses, y = after_stat(density)), binwidth = 1, col = "black", fill = "grey90") + 
  geom_vline(aes(xintercept = mean(gains)), col = "black", linetype = "dashed") +
  geom_vline(aes(xintercept = mean(losses)), col = "black", linetype = "dashed") +
  geom_vline(aes(xintercept = median(gains)), col = "black", linetype = "solid") +
  geom_vline(aes(xintercept = median(losses)), col = "black", linetype = "solid") +
  xlab("events") +
  xlim(0, 80) +
  theme_bw()

ggsave('~/Desktop/cdtb/gainloss_mixedrates_pi0.pdf', gain_loss_30, width = 6, height = 3, dpi = 300)

# Summarise simmaps
summary_maps <- describe.simmap(simmap)

## To color branches: Get dominant state per edge averaged across simmaps
edge_states <- do.call(cbind, lapply(simmap, function(s) {
  sapply(s$maps, function(x) names(which.max(x)))
}))

confidence <- apply(edge_states, 1, function(x) {
  max(table(x)) / length(x)
})

hist(confidence)

# Set color scheme for each state
color_full <- c("0" = "#B3B3B3", 
                 "1" = "#1E466E",
                 "2" = "#315F8C", 
                 "3" = "#3B6D99", 
                 "4" = "#467DA2",
                 "5" = "#508DAB",
                 "6" = "#5D9EBA",
                 "7" = "#69B0CA",
                 "8" = "#79C0D6",
                 "9" = "#8FCCDA",
                 "10" = "#9BD5DD",
                 "11" = "#B4DDD9",
                 "12" = "#D4E0CB")

# Most common state per edge across all simulations
dominant_branch_state <- apply(edge_states, 1, function(x) {
  names(which.max(table(x)))
})

edge_df <- data.frame(
  node  = simmap[[1]]$edge[, 2],
  state = dominant_branch_state
)

## Ancestral state reconstruction at nodes 
ace_result <- ancr(fit, tips = TRUE)

# Plot
plot(ace_result, args.plotTree = list(ftype = "i"), type = "fan")

ancestral_probs <- ace_result$ace

# Convert to data frame with node column
ancestral_df <- as.data.frame(ancestral_probs)
ancestral_df$node <- as.numeric(rownames(ancestral_probs))
write.csv(ancestral_df, '~/Desktop/cdtb/ancestralprobs_synteny_mixedrates_pi0.csv')
ancestral_df <- read.csv('~/Desktop/cdtb/ancestralprobs_synteny_mixedrates_pi0.csv')

# Reorder to put node column first
ancestral_df <- ancestral_df[, c("node", colnames(ancestral_probs))]

print(head(ancestral_df))

# Convert to data frame and add node column
n_tips <- length(newtree$tip.label)
trait_anc <- as.data.frame(ancestral_probs)
trait_anc$node <- suppressWarnings(as.numeric(rownames(ancestral_probs)))

# The ancestral_probs includes BOTH tips and internal nodes
# Separate them

trait_anc_filtered <- trait_anc %>%
  filter(rownames(.) %in% newtree$tip.label | `0` <= 0.90)

trait_anc_filtered$node <- ifelse(is.na(trait_anc_filtered$node),
                                   match(rownames(trait_anc_filtered), newtree$tip.label),
                                   trait_anc_filtered$node)

print(paste("Original nodes:", nrow(trait_anc)))
print(paste("Filtered nodes:", nrow(trait_anc_filtered)))


# Create pies only for filtered nodes
n_states <- ncol(trait_anc_filtered)   # excluding node column
pies_filtered <- nodepie(trait_anc_filtered, 
                          cols = 1:13,
                          color = color_full, 
                          alpha = 1,
                          outline.color = "black",
                          outline.size = 0.1)

# Create df with only filtered nodes
df_filtered <- tibble(node = as.numeric(trait_anc_filtered$node), 
                       pies = pies_filtered)

p <- revts(ggtree(newtree, layout = "rectangular", open.angle = 180)) %<+% df_filtered %<+% edge_df +
  geom_tree(aes(colour = state)) +
  scale_colour_manual(values = color_full) +
  geom_plot(data = td_filter(!isTip), 
            mapping = aes(x = x, y = y, label = pies), 
            vp.width = 0.02, 
            vp.height = 0.02,
            hjust = 0.6,
            vjust = 0.55) +
  geom_plot(data = td_filter(isTip), 
            mapping = aes(x = x, y = y, label = pies),
            vp.width = 0.013, 
            vp.height = 0.013,
            hjust = 0.8,
            vjust = 0.55) +
  theme_tree2(legend.position = 'none') 

ggsave('~/Desktop/cdtb/260830_ASRpi0_mixedrates.pdf', p, height = 20, width = 10, dpi = 300)

# Collapse nodes without cdtB acquisition events 
pcollapsed <- scaleClade(p, 512, 0.2) %>% 
  scaleClade(537, 0.2) %>% 
  scaleClade(454, 0.2) %>% 
  scaleClade(418, 0.2) %>% 
  scaleClade(486, 0.2) %>% 
  scaleClade(629, 0.2) %>% 
  scaleClade(793, 0.2) %>% 
  scaleClade(760, 0.2) %>% 
  scaleClade(674, 0.2) %>% 
  scaleClade(722, 0.2) %>% 
  scaleClade(738, 0.2) %>% 
  scaleClade(800, 0.2) %>% 
  scaleClade(810, 0.2) %>% 
  ggtree::collapse(512, 'max', color = "#B3B3B3", fill= "#D9D9D9") %>%
  ggtree::collapse(537, 'max', color = "#B3B3B3", fill= "#D9D9D9") %>%
  ggtree::collapse(454, 'max', color = "#B3B3B3", fill= "#D9D9D9") %>%
  ggtree::collapse(418, 'max', color = "#B3B3B3", fill= "#D9D9D9") %>%
  ggtree::collapse(486, 'max', color = "#B3B3B3", fill= "#D9D9D9") %>%
  ggtree::collapse(629, 'max', color = "#B3B3B3", fill= "#D9D9D9") %>% 
  ggtree::collapse(793, 'max', color = "#B3B3B3", fill= "#D9D9D9") %>%
  ggtree::collapse(760, 'max', color = "#B3B3B3", fill= "#D9D9D9") %>%
  ggtree::collapse(722, 'max', color = "#B3B3B3", fill= "#D9D9D9") %>%
  ggtree::collapse(738, 'max', color = "#B3B3B3", fill= "#D9D9D9") %>%
  ggtree::collapse(800, 'max', color = "#B3B3B3", fill= "#D9D9D9") %>%
  ggtree::collapse(810, 'max', color = "#B3B3B3", fill= "#D9D9D9") %>%
  ggtree::collapse(674, 'max', color = "#B3B3B3", fill= "#D9D9D9") 

ggsave('~/Desktop/cdtb/260724_ASRpi0_mixedrates_collapsed.pdf', p6collapsed, height = 20, width = 10, dpi = 300)

## Estimating acquisition ages for all 0 -> N transitions

# Calculate transition times across tree

get_transitions <- function(s, to_state = NULL, from_state = NULL) {
  transitions <- data.frame()
  node_ages <- branching.times(s)
  
  for (i in 1:nrow(s$edge)) {
    parent <- s$edge[i, 1]
    child  <- s$edge[i, 2]
    maps   <- s$maps[[i]]
    states <- names(maps)
    
    if (length(states) < 2) next  # no transition on this branch
    
    for (j in 1:(length(states) - 1)) {
      from <- states[j]
      to   <- states[j + 1]
      
      if (!is.null(from_state) && from != from_state) next
      if (!is.null(to_state)   && to   != to_state)   next
      
      parent_age     <- node_ages[as.character(parent)]
      child_age      <- ifelse(child <= Ntip(s), 0, node_ages[as.character(child)])
      time_into_branch <- sum(maps[1:j])
      transition_age <- parent_age - time_into_branch
      
      transitions <- rbind(transitions, data.frame(
        parent         = parent,
        child          = child,
        from           = from,
        to             = to,
        transition_age = transition_age,
        parent_age     = parent_age,
        child_age      = child_age
      ))
    }
  }
  return(transitions)
}

all_events <- do.call(rbind, lapply(seq_along(simmap), function(i) {
  gains  <- get_transitions(simmap[[i]], from_state = "0")
  losses <- get_transitions(simmap[[i]], to_state   = "0")
  gains$event_type  <- "gain"
  losses$event_type <- "loss"
  df <- rbind(gains, losses)
  df$sim <- i
  df
}))

event_summary <- all_events %>%
  group_by(parent, child, event_type, to) %>%  # add to here
  summarise(
    mean_age  = mean(transition_age),
    lower_90 = quantile(transition_age, 0.05),
    upper_90 = quantile(transition_age, 0.950),
    lower_95  = quantile(transition_age, 0.025),
    upper_95  = quantile(transition_age, 0.975),
    n_sims    = n(),
    .groups   = "drop"
  ) %>%
  mutate(node = child)

# For each simulation, find the oldest gain per state
earliest_per_sim <- all_events %>%
  filter(event_type == "gain") %>%
  group_by(sim, to) %>%
  slice_max(transition_age, n = 1, with_ties = FALSE) %>%
  ungroup()

# The earliest age supported by ALL sims = minimum of these oldest ages
earliest_supported <- earliest_per_sim %>%
  group_by(to) %>%
  summarise(
    min_age  = max(transition_age),   # oldest age present in every sim
    mean_age = mean(transition_age),  # for reference
    lower_95 = quantile(transition_age, 0.025),
    upper_95 = quantile(transition_age, 0.975),
    n_sims   = n(),                   # should equal 100 if state gained in all sims
    .groups  = "drop"
  )

# Reorder such that data is displayed oldest -> youngest transition
earliest_supported <- earliest_supported %>%
  mutate(to = reorder(to, mean_age)) 

# Plot acquisition times
ages <- ggplot(earliest_supported,
                aes(x = -(mean_age), y = to, colour = to)) +
  geom_errorbarh(
    aes(xmin = -(lower_95), xmax = -(upper_95)),
    height = 0.3, linewidth = 0.5
  ) +
  geom_point(aes(fill = to),
             shape = 21, colour = "white", size = 3, stroke = 0.3
  ) +
  scale_fill_manual(values = color_full, name = "State") +
  scale_colour_manual(values = color_full, name = "State") +
  labs(x = "Age (Ma)", y = "State") +
  xlim(-80, 0) + 
  theme(legend.position = 'none') +
  theme_classic()

ggsave('~/Desktop/cdtb/ASRpi0_mixedrates_ages.pdf', ages4, width = 10, height = 4, dpi = 300)

## Diversification rate analyses 

# Sister clade analysis, see Table SX for group size, stem age, and references for each group. 
# cdtB comparisons
# All sister clade comparisons
derived_species <- c(45, 25, 10, 329, 3079, 35)
ancestral_species <- c(4, 32, 8, 482, 107, 22)
pairs <- cbind(derived_species, ancestral_species)
rownames(pairs) <- c("ananassae/setifemur", "cardini/guarani", 
                     "saltans_saltans/saltans_cordataelliptica", 
                     "Leucophenga/Amiota_Cacoxenus_Phortica", 
                     "SubgenusDrosophila+/SubgenusDorsilopha_Microdrosophila",
                     "Rhinoleucophenga/Gitona")
div_times <- c(24.0026, 13.336, 8.2272, 49, 39.123, 39.3) 
richness.yule.test(pairs, div_times)

# drop ananassae
derived_species2 <- c(25, 10, 329, 3079, 35)
ancestral_species2 <- c(32, 8, 482, 107, 22)
pair26 <- cbind(derived_species2, ancestral_species2)
div_times2 <- c(13.336, 8.2272, 49, 39.123, 39.3) 
richness.yule.test(pairs2, div_times2)

# drop cardini
derived_species3 <- c(45, 10, 329, 3079, 35)
ancestral_species3 <- c(4, 8, 482, 107, 22)
pairs3 <- cbind(derived_species3, ancestral_species3)
div_times3 <- c(24.0026, 8.2272, 49, 39.123, 39.3) 
richness.yule.test(pairs3, div_times3)

# drop saltans
derived_species4 <- c(45, 25, 329, 3079, 35)
ancestral_species4 <- c(4, 32, 482, 107, 22)
pairs4 <- cbind(derived_species4, ancestral_species4)
div_times4 <- c(24.0026, 13.336, 49, 39.123, 39.3) 
richness.yule.test(pairs4, div_times4)

# drop leucophenga
derived_species5<- c(45, 25, 10, 3079, 35)
ancestral_species5 <- c(4, 32, 8, 107, 22)
pairs5 <- cbind(derived_species5, ancestral_species5)
div_times5 <- c(24.0026, 13.336, 8.2272, 39.123, 39.3) 
richness.yule.test(pairs5, div_times5)

#drop subgenus drosophila
derived_species6 <- c(45, 25, 10, 329, 35)
ancestral_species6 <- c(4, 32, 8, 482, 22)
pairs6 <- cbind(derived_species6, ancestral_species6)
div_times6 <- c(24.0026, 13.336, 8.2272, 49, 39.3) 
richness.yule.test(pairs6, div_times6)

# drop rhinoleucophenga
derived_species7 <- c(45, 25, 10, 329, 3079)
ancestral_species7 <- c(4, 32, 8, 482, 107)
pairs7 <- cbind(derived_species7, ancestral_species7)
div_times7 <- c(24.0026, 13.336, 8.2272, 49, 39.123) 
richness.yule.test(pairs7, div_times7)

#PPO pair comparison 
derived_ppo <- c(102, 74, 24, 1485, 8)
ancestral_ppo <- c(126, 99, 123, 5024, 45)
pairs_ppo <- cbind(derived_ppo2, ancestral_ppo2)
rownames(pairs_ppo) <- c("ananassae/setifemur", "cardini/guarani", 
                     "saltans_saltans/saltans_cordataelliptica", 
                     "Leucophenga/Amiota_Cacoxenus_Phortica", 
                     "SubgenusDrosophila+/SubgenusDorsilopha_Microdrosophila",
                     "Rhinoleucophenga/Gitona")
div_times_ppo <- c(22.173, 31.297, 21.8046, 67.3, 16.6204) 
richness.yule.test(pairs_ppo, div_times_ppo)
# Not significant

# PGLS analysis for cdtB clades

# Estimate diversification rates for each clade under three different extinction scenarios
## bd.ms estimated diversification rates
# ananassae 
bd.ms(time=24.0026, n=45, missing = 30, crown = FALSE, epsilon = 0)
bd.ms(time=24.0026, n=45, missing = 30, crown = FALSE, epsilon = 0.5)
bd.ms(time=24.0026, n=45, missing = 30, crown = FALSE, epsilon = 0.9)
# setifemur
bd.ms(time=24.0026, n=4, missing = 3, crown = FALSE, epsilon = 0)
bd.ms(time=24.0026, n=4, missing = 3, crown = FALSE, epsilon = 0.5)
bd.ms(time=24.0026, n=4, missing = 3, crown = FALSE, epsilon = 0.9)
# drosophila + others
bd.ms(time=39.0717, n=3079, missing = 2851, crown = FALSE, epsilon = 0)
bd.ms(time=39.0717, n=3079, missing = 2851, crown = FALSE, epsilon = 0.5)
bd.ms(time=39.0717, n=3079, missing = 2851, crown = FALSE, epsilon = 0.9)
# dorsilopha + others
bd.ms(time=39.0717, n=107, missing = 105, crown = FALSE, epsilon = 0)
bd.ms(time=39.0717, n=107, missing = 105, crown = FALSE, epsilon = 0.5)
bd.ms(time=39.0717, n=107, missing = 105, crown = FALSE, epsilon = 0.9)
# cardini
bd.ms(time=13.3254, n=25, missing = 18, crown = FALSE, epsilon = 0)
bd.ms(time=13.3254, n=25, missing = 18, crown = FALSE, epsilon = 0.5)
bd.ms(time=13.3254, n=25, missing = 18, crown = FALSE, epsilon = 0.9)
# guarani
bd.ms(time=13.3254, n=32, missing = 30, crown = FALSE, epsilon = 0)
bd.ms(time=13.3254, n=32, missing = 30, crown = FALSE, epsilon = 0.5)
bd.ms(time=13.3254, n=32, missing = 30, crown = FALSE, epsilon = 0.9)
# saltans_saltans
bd.ms(time=8.2272, n=10, missing = 5, crown = FALSE, epsilon = 0)
bd.ms(time=8.2272, n=10, missing = 5, crown = FALSE, epsilon = 0.5)
bd.ms(time=8.2272, n=10, missing = 5, crown = FALSE, epsilon = 0.9)
# saltans_neocordataetc
bd.ms(time=8.2272, n=8, missing = 4, crown = FALSE, epsilon = 0)
bd.ms(time=8.2272, n=8, missing = 4, crown = FALSE, epsilon = 0.5)
bd.ms(time=8.2272, n=8, missing = 4, crown = FALSE, epsilon = 0.9)
# leucophenga
bd.ms(time=49, n=329, missing = 324, crown = FALSE, epsilon = 0)
bd.ms(time=49, n=329, missing = 324, crown = FALSE, epsilon = 0.5)
bd.ms(time=49, n=329, missing = 324, crown = FALSE, epsilon = 0.9)
# leucophenga sister
bd.ms(time=49, n=482, missing = 474, crown = FALSE, epsilon = 0)
bd.ms(time=49, n=482, missing = 474, crown = FALSE, epsilon = 0.5)
bd.ms(time=49, n=482, missing = 474, crown = FALSE, epsilon = 0.9)
# rhinoleucophenga
bd.ms(time=36.1139, n=35, missing = 33, crown = FALSE, epsilon = 0)
bd.ms(time=36.1139, n=35, missing = 33, crown = FALSE, epsilon = 0.5)
bd.ms(time=36.1139, n=35, missing = 33, crown = FALSE, epsilon = 0.9)
# gitona
bd.ms(time=36.1139, n=22, missing = 21, crown = FALSE, epsilon = 0)
bd.ms(time=36.1139, n=22, missing = 21, crown = FALSE, epsilon = 0.5)
bd.ms(time=36.1139, n=22, missing = 21, crown = FALSE, epsilon = 0.9)

# Prune tree for PGLS
pgls_tips <- c('DROSOPHILA_PRIMAEVA', 'DROSOPHILA_DUNNI', "DROSOPHILA_SUBBADIA", "DROSOPHILA_BUSCKII", 
               "DROSOPHILA_ANANASSAE", "DROSOPHILA_SETIFEMUR", "DROSOPHILA_SALTANS", "DROSOPHILA_NEOCORDATA",
               "LEUCOPHENGA_VARIA", "AMIOTA_MINOR", "RHINOLEUCOPHENGA_AMERICANA", "GITONA_DISTIGMA")

pglstree <- keep.tip(newtree, pgls_tips)
plot(pglstree)

# Load data 
pglsdata <- read.csv('PGLS_input.csv') # data calculated from bd.ms above

# Load packages 
library(caper)

# Perform PGLS across all three extinction rate classes 
pglsall <- comparative.data(pglstree, pglsdata, label, vcv=FALSE, vcv.dim=2, na.omit=TRUE, 
                            force.root=FALSE, warn.dropped=FALSE, scope=NULL)

# Extinction = 0 (no extinction)
f_0 <- pgls(bd.ms_0 ~ cdtb_species, pglsall, lambda = "ML", kappa = 1, delta = 1, param.CI = 0.95)
anova.pgls(f_0)
summary.pgls(f_0)

# Moderate extinction (extinction = 0.5)
f_0.5 <- pgls(bd.ms_0.5 ~ cdtb_species, pglsall, lambda = "ML", kappa = 1, delta = 1, param.CI = 0.95)
anova.pgls(f_0.5)
summary.pgls(f_0.5)

# High extinction (extinction = 0.9)
f_0.9 <- pgls(bd.ms_0.9 ~ cdtb_species, pglsall, lambda = "ML", kappa = 1, delta = 1, param.CI = 0.95)
anova.pgls(f_0.9)
summary.pgls(f_0.9)

# Lineage through time analyses, adapted from Liam Revell

# Mapping lineages in a given state over time 
# Factorize tip states
tip_states_fac <- as.factor(tip_states_vec)

foo<-function(newtree,x){
  tt<-map.to.singleton(newtree)
  H<-nodeHeights(tt)
  h<-max(H)-branching.times(tt)
  ss<-setNames(as.factor(names(tt$edge.length)),
               tt$edge[,2])
  lineages<-matrix(0,length(h),length(levels(x)),
                   dimnames=list(names(h),levels(x)))
  for(i in 1:length(h)){
    ii<-intersect(which(h[i]>H[,1]),which(h[i]<=H[,2]))
    lineages[i,]<-summary(ss[ii])
  }
  ii<-order(h)
  times<-h[ii]
  lineages<-lineages[ii,]
  list(times=times,ltt=lineages)
}

ltts<-lapply(simmap,foo,x=tip_states_fac)

group_binary <- function(ltt_result) {
  ltt_mat <- ltt_result$ltt
  state0_col <- which(colnames(ltt_mat) == "0")
  
  state0    <- ltt_mat[, state0_col]
  state_pos <- rowSums(ltt_mat[, -state0_col, drop = FALSE])   # sum all non-"0" columns
  
  list(times = ltt_result$times,
       ltt   = cbind("0" = state0, ">0" = state_pos))
}

ltts_binary <- lapply(ltts, group_binary)

# Convert to dataframe to plot with ggplot2 
ltt_long <- do.call(rbind, lapply(seq_along(ltts_binary), function(i) {
  rbind(
    data.frame(tree = i, time = ltts_binary[[i]]$times,
               lineages = ltts_binary[[i]]$ltt[, "0"],  state = "0"),
    data.frame(tree = i, time = ltts_binary[[i]]$times,
               lineages = ltts_binary[[i]]$ltt[, ">0"], state = ">0")
  )
}))

# Define a common grid across the full time range
time_grid <- seq(min(sapply(ltts, function(x) min(x$times))),
                 max(sapply(ltts, function(x) max(x$times))),
                 length.out = 500)

# For each tree, build a step function and evaluate it at every grid point
interp_state <- function(x, col_name, grid) {
  sf <- stepfun(x$times[-1], x$ltt[, col_name])
  sf(grid)
}

mat1 <- sapply(ltts_binary, interp_state, col_name = "0",  grid = time_grid)
mat2 <- sapply(ltts_binary, interp_state, col_name = ">0", grid = time_grid)

# Calculate total lineages per tree = state1 + state2, computed BEFORE summarizing
mat_all <- mat1 + mat2

# Summarize all three states

summarize_ci <- function(mat, grid, state_label) {
  data.frame(
    time  = grid,
    mean  = rowMeans(mat, na.rm = TRUE),
    lower = apply(mat, 1, quantile, probs = 0.025, na.rm = TRUE),
    upper = apply(mat, 1, quantile, probs = 0.975, na.rm = TRUE),
    state = state_label
  )
}

ci_df <- rbind(
  summarize_ci(mat1,   time_grid, "cdtb-"),
  summarize_ci(mat2,   time_grid, "cdtb+"),
  summarize_ci(mat_all, time_grid, "all")
)

present <- max(sapply(ltts, function(x) max(x$times)))

ci_df$time_bp <- ci_df$time - present   # 0 = present, negative = further in the past

# Guard against zeros before log-transforming (log(0) = -Inf)
ci_df$mean[ci_df$mean <= 0]   <- NA
ci_df$lower[ci_df$lower <= 0] <- 1   # lower CI bound of 0 is common/expected; NA lets ggplot drop it cleanly on log scale


# Catch NAs too, not just <= 0
true_zero <- is.na(ci_df$upper) | ci_df$upper <= 0

# Only floor the LOWER bound where there's actually a positive mean/upper 
# (i.e., real data exists, just the lower CI quantile dipped to/below 0)
ci_df$lower[!true_zero & (is.na(ci_df$lower) | ci_df$lower <= 0)] <- 1

# Leave true-zero rows as NA for both bounds — there's genuinely nothing to plot there
ci_df$lower[true_zero] <- NA
ci_df$upper[true_zero] <- NA

# Plot
state_colors <- c("cdtb-" = "#B3B3B3", "cdtb+" = "#72BCD5", "all" = "black")

ltt_states <- ggplot(ci_df, aes(x = time_bp, y = mean, color = state, fill = state)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_step(linewidth = 1.2) +
  geom_vline(xintercept = -44.32, color = "black", linewidth = 0.5, linetype = "dashed") + #cdtb acquisition in subgenus drosophila +
  scale_color_manual(values = state_colors, breaks = c("cdtb-", "cdtb+", "all")) +
  scale_fill_manual(values = state_colors, breaks = c("cdtb-", "cdtb+", "all")) +
  scale_y_log10(breaks = c(1, 2, 5, 10, 20, 50, 100, 200, 500)) +
  coord_cartesian(ylim = c(1, 500)) +
  labs(x = "time", y = "lineages", color = NULL, fill = NULL) +
  theme_classic() +
  theme(legend.position = c(0.15, 0.85))

ggsave('LTTstates_260831.pdf', ltt_states, width = 10, height = 4, dpi = 300)
