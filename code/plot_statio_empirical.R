# This script plots relative clade age, clade size, and distance to stationary across phylogenies 
# using inference result (for parameter estimates) from BiSSE analysis 

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

#### PLOT ####
relative_age          = comb_dat$rel_age
relative_size         = comb_dat$rel_size
distance_from_statio  = comb_dat$distance

col_group = rep(c("bats","squamates","grass"),times = c(nrow(comb_dat[comb_dat$phylogeny=="bats_tree",]),
                                                        nrow(comb_dat[comb_dat$phylogeny=="squamate_tree",]),
                                                        nrow(comb_dat[comb_dat$phylogeny=="grass_tree",])))

colors <- c("blue", "red","orange")  # one color per group 


plot_ly(x = ~relative_age, 
        y = ~relative_size, 
        z = ~distance_from_statio,
        type = "scatter3d", 
        mode = "markers",
        color = ~col_group,
        colors = colors) %>%    # color by group
  layout(scene = list(
    xaxis = list(title = "rel.clade age to root", showbackground = FALSE),
    yaxis = list(title = "rel.clade size to tree", showbackground = FALSE),
    zaxis = list(title = "distance from statio", showbackground = FALSE)
  ))
