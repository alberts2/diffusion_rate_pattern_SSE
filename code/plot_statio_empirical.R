# This script plots relative clade age, clade size, and distance to stationary across phylogenies 
# using inference result (for parameter estimates) from BiSSE analysis 

# NOTE: 26 Sep
# Change Euclidean distance to absolute distance in one of the states
# added color-code according to evolutionary scenarios

#### LOAD LIBRARIES ####
library(ggplot2)
library(ape)
library(scatterplot3d)

#### FILE SETTING ####
fp     = "/Users/albertsoewongsono/Documents/Code\ Testing/rate_pattern_diffusion_SSE/"
in_fp  = paste0(fp,"data/Inference/")
dat_fp = paste0(fp,"data/Empirical/")
out_fp = paste0(fp,"plot/")

#### LOAD THE TREES FILES ####
# BATS
bats_tre     = read.tree(paste0(dat_fp,"bats_phylogeny/bats_mcc.tre"))
bats_size    = length(bats_tre$tip.label)
bats_age     = max(node.depth.edgelength(bats_tre))
# SQUAMATES
squamate_tre  = read.tree(paste0(dat_fp,"squamate_phylogeny/squamate.txt"))  
squamate_size = length(squamate_tre$tip.label)
squamate_age  = max(node.depth.edgelength(squamate_tre))
# GRASS
grass_tre  = read.tree(paste0(dat_fp,"grass_phylogeny/Poaceae_hyp1_allC4_maxcred.phy"))  
grass_size = length(grass_tre$tip.label)
grass_age  = max(node.depth.edgelength(grass_tre))

#### LOAD INFERENCE FILE ####
# BATS
bats_inf     = read.table(paste0(in_fp,"bats_summary.csv"),sep=";",header = T)
bats_inf     = cbind(bats_inf,rep("bats_tree",nrow(bats_inf)))
colnames(bats_inf)[15] = "phylogeny"
# SQUAMATE 
squamate_inf = read.table(paste0(in_fp,"squamate_summary.csv"),sep=";",header = T)
squamate_inf = cbind(squamate_inf,rep("squamate_tree",nrow(squamate_inf)))
colnames(squamate_inf)[15] = "phylogeny"
# GRASS 
grass_inf = read.table(paste0(in_fp,"grass_summary.csv"),sep=";",header = T)
grass_inf = cbind(grass_inf,rep("grass_tree",nrow(grass_inf)))
colnames(grass_inf)[15] = "phylogeny"

#### ANALYSIS ####
# Combine the inference file
comb_dat = rbind(bats_inf,squamate_inf,grass_inf)
# add relative clade age to root age (note: small number <<1 means it's a young clade)
comb_dat$rel_age  = NA
# add relative clade size to tree size (note: small number <<1 means it's a small clade)
comb_dat$rel_size = NA
# fill in the entries for the newly added columns 
for (i in 1:nrow(comb_dat)){
  if (comb_dat[i,]$phylogeny == "bats_tree"){
    comb_dat[i,]$rel_age  = comb_dat[i,]$root_age/bats_age
    comb_dat[i,]$rel_size = comb_dat[i,]$tip_count/bats_size
  } else if (comb_dat[i,]$phylogeny == "squamate_tree"){
    comb_dat[i,]$rel_age  = comb_dat[i,]$root_age/squamate_age
    comb_dat[i,]$rel_size = comb_dat[i,]$tip_count/squamate_size
  } else if (comb_dat[i,]$phylogeny == "grass_tree"){
    comb_dat[i,]$rel_age  = comb_dat[i,]$root_age/grass_age
    comb_dat[i,]$rel_size = comb_dat[i,]$tip_count/grass_size
  }
}

# Using absolute distance for one of the states instead of euclidean distance 
comb_dat$distance = NULL
comb_dat$distance = rep(NA,nrow(comb_dat))
for (i in 1:nrow(comb_dat)){
  comb_dat$distance[i] = abs(comb_dat$emp_0[i]-comb_dat$statio_0[i])
}

#### CLASSIFYING ####
## For statio_A < statio_B
# non_trivial cases (under our approach)
# evo_scenario 1a: species loss in B, species gain in A, more transition into B
# evo_scenario 2a: species gain in B, species loss in A
# trivial cases
# evo_scenario 3a: equal net diversification in A&B, more transition into B 
# evo_scenario 4a: species gain in B, species gain in A, more transition into B 
# evo_scenario 5a: species loss in B, species loss in A, more transition into B
# evo_scenario 6a: species gain in B, species gain in A, more transition into A
# evo_scenario 7a: species loss in B, species loss in A, more transition into A

## For statio_A > statio_B
# non_trivial cases (under our approach)
# evo_scenario 1b: species gain in B, species loss in A, more transition into A
# evo_scenario 2b: species loss in B, species gain in A 
# trivial cases
# evo_scenario 3b: equal net diversification in A&B, more transition into A 
# evo_scenario 4b: species gain in B, species gain in A, more transition into A 
# evo_scenario 5b: species loss in B, species loss in A, more transition into A
# evo_scenario 6b: species gain in B, species gain in B, more transition into B
# evo_scenario 7b: species loss in B, species loss in A, more transition into B 

## For statio_A = statio_B 
# evo_scenario 1c: species gain in B, species loss in A, more transition into A 
# evo_scenario 2c: species loss in B, species gain in A, more transition into B
# evo_scenario 3c: equal speciation, equal extinction, equal transition rates 

comb_dat$summary = rep(NA,nrow(comb_dat))

for (i in 1:nrow(comb_dat)){ #loop over each clade
  sign_obs    = sign(comb_dat$emp_0[i] - comb_dat$emp_1[i]) #sign for observed frequency (1 means \Pi_A > \Pi_B)
  sign_statio = sign(comb_dat$statio_0[i] - comb_dat$statio_1[i]) #sign for statio frequency (\Pi_hat_A - \Pi_hat_B)
  #
  dir_net_diver_A = comb_dat$lambda_0[i] - comb_dat$mu_0[i] # species loss or gain in A
  dir_net_diver_B = comb_dat$lambda_1[i] - comb_dat$mu_1[i] # species loss or gain in B
  dir_transition  = comb_dat$q_01[i] - comb_dat$q_10[i] # q_01 > q_10 or q_01 < q_10
  #
  if (sign_obs == sign_statio){ #if the direction of difference in observed frequencies matches direction of difference in statio freq
    # then it is an outlier, and we assign which of the evolutionary scenario it corresponds to
    if (dir_net_diver_B < 0 & dir_net_diver_A > 0 & dir_transition > 0 & sign_statio == -1){
      comb_dat$summary[i] = "not outlier & scenario 1a"
    } else if (dir_net_diver_B > 0 & dir_net_diver_A < 0 & sign_statio == -1){
      comb_dat$summary[i] = "not outlier & scenario 2a"
    } else if (dir_net_diver_A == dir_net_diver_B & dir_transition > 0 & sign_statio == -1){
      comb_dat$summary[i] = "not outlier & scenario 3a"
    } else if (dir_net_diver_B > 0 & dir_net_diver_A > 0 & dir_transition > 0 & sign_statio == -1){
      comb_dat$summary[i] = "not outlier & scenario 4a"
    } else if (dir_net_diver_B < 0 & dir_net_diver_A < 0 & dir_transition > 0 & sign_statio == -1){
      comb_dat$summary[i] = "not outlier & scenario 5a"
    } else if (dir_net_diver_B > 0 & dir_net_diver_A > 0 & dir_transition < 0 & sign_statio == -1){
      comb_dat$summary[i] = "not outlier & scenario 6a"
    } else if (dir_net_diver_B < 0 & dir_net_diver_A < 0 & dir_transition < 0 & sign_statio == -1){
      comb_dat$summary[i] = "not outlier & scenario 7a"
    } else if (dir_net_diver_B > 0 & dir_net_diver_A < 0 & dir_transition < 0 & sign_statio == 1){
      comb_dat$summary[i] = "not outlier & scenario 1b"
    } else if (dir_net_diver_B < 0 & dir_net_diver_A > 0 & sign_statio == 1){
      comb_dat$summary[i] = "not outlier & scenario 2b"
    } else if (dir_net_diver_A == dir_net_diver_B & dir_transition < 0 & sign_statio == 1){
      comb_dat$summary[i] = "not outlier & scenario 3b"
    } else if (dir_net_diver_B > 0 & dir_net_diver_A > 0 & dir_transition < 0 & sign_statio == 1){
      comb_dat$summary[i] = "not outlier & scenario 4b"
    } else if (dir_net_diver_B < 0 & dir_net_diver_A < 0 & dir_transition < 0 & sign_statio == 1){
      comb_dat$summary[i] = "not outlier & scenario 5b"
    } else if (dir_net_diver_B > 0 & dir_net_diver_A > 0 & dir_transition > 0 & sign_statio == 1){
      comb_dat$summary[i] = "not outlier & scenario 6b"
    } else if (dir_net_diver_B < 0 & dir_net_diver_A < 0 & dir_transition > 0 & sign_statio == 1){
      comb_dat$summary[i] = "not outlier & scenario 7b"
    } else if (dir_net_diver_B > 0 & dir_net_diver_A < 0 & dir_transition < 0 & sign_statio == 0){
      comb_dat$summary[i] = "not outlier & scenario 1c"
    } else if (dir_net_diver_B < 0 & dir_net_diver_A > 0 & dir_transition > 0 & sign_statio == 0){
      comb_dat$summary[i] = "not outlier & scenario 2c"
    } else if (dir_net_diver_B == dir_net_diver_A & dir_transition == 0 &sign_statio == 0){
      comb_dat$summary[i] = "not outlier & scenario 3c"
    }
  } else { #if the directions in obs and statio frequencies do not match -> outliers 
    if (dir_net_diver_B < 0 & dir_net_diver_A > 0 & dir_transition > 0 & sign_statio == -1){
      comb_dat$summary[i] = "outlier & scenario 1a"
    } else if (dir_net_diver_B > 0 & dir_net_diver_A < 0 & sign_statio == -1){
      comb_dat$summary[i] = "outlier & scenario 2a"
    } else if (dir_net_diver_A == dir_net_diver_B & dir_transition > 0 & sign_statio == -1){
      comb_dat$summary[i] = "outlier & scenario 3a"
    } else if (dir_net_diver_B > 0 & dir_net_diver_A > 0 & dir_transition > 0 & sign_statio == -1){
      comb_dat$summary[i] = "outlier & scenario 4a"
    } else if (dir_net_diver_B < 0 & dir_net_diver_A < 0 & dir_transition > 0 & sign_statio == -1){
      comb_dat$summary[i] = "outlier & scenario 5a"
    } else if (dir_net_diver_B > 0 & dir_net_diver_A > 0 & dir_transition < 0 & sign_statio == -1){
      comb_dat$summary[i] = "outlier & scenario 6a"
    } else if (dir_net_diver_B < 0 & dir_net_diver_A < 0 & dir_transition < 0 & sign_statio == -1){
      comb_dat$summary[i] = "outlier & scenario 7a"
    } else if (dir_net_diver_B > 0 & dir_net_diver_A < 0 & dir_transition < 0 & sign_statio == 1){
      comb_dat$summary[i] = "outlier & scenario 1b"
    } else if (dir_net_diver_B < 0 & dir_net_diver_A > 0 & sign_statio == 1){
      comb_dat$summary[i] = "outlier & scenario 2b"
    } else if (dir_net_diver_A == dir_net_diver_B & dir_transition < 0 & sign_statio == 1){
      comb_dat$summary[i] = "outlier & scenario 3b"
    } else if (dir_net_diver_B > 0 & dir_net_diver_A > 0 & dir_transition < 0 & sign_statio == 1){
      comb_dat$summary[i] = "outlier & scenario 4b"
    } else if (dir_net_diver_B < 0 & dir_net_diver_A < 0 & dir_transition < 0 & sign_statio == 1){
      comb_dat$summary[i] = "outlier & scenario 5b"
    } else if (dir_net_diver_B > 0 & dir_net_diver_A > 0 & dir_transition > 0 & sign_statio == 1){
      comb_dat$summary[i] = "outlier & scenario 6b"
    } else if (dir_net_diver_B < 0 & dir_net_diver_A < 0 & dir_transition > 0 & sign_statio == 1){
      comb_dat$summary[i] = "outlier & scenario 7b"
    } else if (dir_net_diver_B > 0 & dir_net_diver_A < 0 & dir_transition < 0 & sign_statio == 0){
      comb_dat$summary[i] = "outlier & scenario 1c"
    } else if (dir_net_diver_B < 0 & dir_net_diver_A > 0 & dir_transition > 0 & sign_statio == 0){
      comb_dat$summary[i] = "outlier & scenario 2c"
    } else if (dir_net_diver_B == dir_net_diver_A & dir_transition == 0 &sign_statio == 0){
      comb_dat$summary[i] = "outlier & scenario 3c"
    }
  } 
}
#### PLOT ####

# group the clades according to their phylogeny 
comb_dat$phylogeny <- factor(
  comb_dat$phylogeny,
  levels = c("bats_tree", "squamate_tree", "grass_tree")  
)

# plot_ly(x = ~relative_age, 
#         y = ~relative_size, 
#         z = ~distance_from_statio,
#         type = "scatter3d", 
#         mode = "markers",
#         color = ~col_group,
#         colors = colors) %>%    # color by group
#   layout(scene = list(
#     xaxis = list(title = "rel.clade age to root", showbackground = FALSE),
#     yaxis = list(title = "rel.clade size to tree", showbackground = FALSE),
#     zaxis = list(title = "distance from statio", showbackground = FALSE)
#   ))

plot_ly(comb_dat,
        x = ~rel_age, 
        y = ~rel_size, 
        z = ~distance,
        type   = "scatter3d", 
        mode   = "markers",
        color  = ~summary,
        symbol = ~phylogeny,
        symbols = c("circle", "square", "diamond")) %>%    # color by group
  layout(
    scene = list(
      xaxis = list(title = "rel.clade age to root", showbackground = FALSE),
      yaxis = list(title = "rel.clade size to tree", showbackground = FALSE),
      zaxis = list(title = "distance from statio", showbackground = FALSE)
    ),
    legend = list(traceorder = "normal")  # preserve the factor level order
  )

