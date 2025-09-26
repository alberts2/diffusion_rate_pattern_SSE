#### LOAD LIBRARIES ####
library(diversitree)
library(ape)
library(phangorn)

#### FILE SETTING ####
# fp     = "/Users/albertsoewongsono/Documents/Code\ Testing/rate_pattern_diffusion_SSE/"
fp     = "/storage/albert/rate_pattern_diffusion_SSE/"
in_fp  = paste0(fp,"data/Empirical/bats_phylogeny/")
out_fp = paste0(fp,"plot/")

setwd(in_fp)

num_starts = 50  #number of likelihood searches
# num_starts = 2  #number of likelihood searches
set.seed(1)

# pars (lambda0,lambda1, mu0, mu1, q01, q10)
jitters_func <- function(pars){
  rand <- runif(length(pars), 0, 1) + 0.1  # one random factor per element
  pars * rand
}

compute_statio <- function(pars){
  g = (pars[1]-pars[3])-(pars[2]-pars[4])
  ans1 = (-(g-pars[5]-pars[6])+ sqrt((g-pars[5]-pars[6])^2+4*g*pars[6]))/(-2*g)
  ans2 = (-(g-pars[5]-pars[6])- sqrt((g-pars[5]-pars[6])^2+4*g*pars[6]))/(-2*g)
  if (ans1 > 0 & ans1<=1){
    return(ans1)
  } else {
    return(ans2)
  }
}

# load tree 
### BATS
Bats  = read.tree("bats_mcc.tre") # mcc tree from posterior tree distribution of bats that Gabriel sent 

# posterior tree distribution across clades 
Emballonuoidea    = read.tree("Emballonuoidea.tre")
Noctilionoidea    = read.tree("Noctilionoidea.tre") 
Pteropodoidea     = read.tree("Pteropodoidea.tre")
Rhinolophoidea    = read.tree("Rhinolophoidea.tre")
Vespertilionoidea = read.tree("Vespertilionoidea.tre")

# get mcc tree from bats clade tree posterior
Emballonuoidea_mcc    = maxCladeCred(Emballonuoidea)
Noctilionoidea_mcc    = maxCladeCred(Noctilionoidea)
Pteropodoidea_mcc     = maxCladeCred(Pteropodoidea)
Rhinolophoidea_mcc    = maxCladeCred(Rhinolophoidea)
Vespertilionoidea_mcc = maxCladeCred(Vespertilionoidea)
# get the species list from each clade 
Emballonuoidea_spec    = Emballonuoidea_mcc$tip.label
Noctilionoidea_spec    = Noctilionoidea_mcc$tip.label
Pteropodoidea_spec     = Pteropodoidea_mcc$tip.label
Rhinolophoidea_spec    = Rhinolophoidea_mcc$tip.label
Vespertilionoidea_spec = Vespertilionoidea_mcc$tip.label
# get the clades from mcc tree of bats phylogeny
Emballonuoidea_mcc     = drop.tip(Bats,setdiff(Bats$tip.label,Emballonuoidea_spec))
Noctilionoidea_mcc     = drop.tip(Bats,setdiff(Bats$tip.label,Noctilionoidea_spec))
Pteropodoidea_mcc      = drop.tip(Bats,setdiff(Bats$tip.label,Pteropodoidea_spec))
Rhinolophoidea_mcc     = drop.tip(Bats,setdiff(Bats$tip.label,Rhinolophoidea_spec))
Vespertilionoidea_mcc  = drop.tip(Bats,setdiff(Bats$tip.label,Vespertilionoidea_spec))

# load bats character data 
# plant       = 1 (included), 0 (excluded)
# litter size = 0 (=1 litter), 1 (>1 litters)
char_dat_bat = read.table("Data_Bisse_Albert.csv",sep=";",header = T)
char_dat_bat = char_dat_bat[-1]

# get the character data for each clade 
Emballonuoidea_char_dat = char_dat_bat[char_dat_bat$species %in% Emballonuoidea_mcc$tip.label,]
Noctilionoidea_char_dat = char_dat_bat[char_dat_bat$species %in% Noctilionoidea_mcc$tip.label,]
Pteropodoidea_char_dat = char_dat_bat[char_dat_bat$species %in% Pteropodoidea_mcc$tip.label,]
Rhinolophoidea_char_dat = char_dat_bat[char_dat_bat$species %in% Rhinolophoidea_mcc$tip.label,]
Vespertilionoidea_char_dat = char_dat_bat[char_dat_bat$species %in% Vespertilionoidea_mcc$tip.label,]

# Associate states with tips 
# species names as row names
rownames(Emballonuoidea_char_dat)    = Emballonuoidea_char_dat$species
rownames(Noctilionoidea_char_dat)    = Noctilionoidea_char_dat$species
rownames(Pteropodoidea_char_dat)     = Pteropodoidea_char_dat$species
rownames(Rhinolophoidea_char_dat)    = Rhinolophoidea_char_dat$species
rownames(Vespertilionoidea_char_dat) = Vespertilionoidea_char_dat$species
# reorder rows to match the order of species on trees 
Emballonuoidea_char_dat    = Emballonuoidea_char_dat[Emballonuoidea_mcc$tip.label,]
Noctilionoidea_char_dat    = Noctilionoidea_char_dat[Noctilionoidea_mcc$tip.label,]
Pteropodoidea_char_dat     = Pteropodoidea_char_dat[Pteropodoidea_mcc$tip.label,]
Rhinolophoidea_char_dat    = Rhinolophoidea_char_dat[Rhinolophoidea_mcc$tip.label,]
Vespertilionoidea_char_dat = Vespertilionoidea_char_dat[Vespertilionoidea_mcc$tip.label,]

# create a vector of tip labels and their states for each character
##### CHARACTER: PLANT IN DIET 
Emballonuoidea_states_plant <- Emballonuoidea_char_dat$plant
names(Emballonuoidea_states_plant) = Emballonuoidea_char_dat$species
#
Noctilionoidea_states_plant <- Noctilionoidea_char_dat$plant
names(Noctilionoidea_states_plant) = Noctilionoidea_char_dat$species
#
Pteropodoidea_states_plant <- Pteropodoidea_char_dat$plant
names(Pteropodoidea_states_plant) = Pteropodoidea_char_dat$species
#
Rhinolophoidea_states_plant <- Rhinolophoidea_char_dat$plant
names(Rhinolophoidea_states_plant) = Rhinolophoidea_char_dat$species
#
Vespertilionoidea_states_plant <- Vespertilionoidea_char_dat$plant
names(Vespertilionoidea_states_plant) = Vespertilionoidea_char_dat$species

##### CHARACTER: LITTER SIZE 
Emballonuoidea_states_litter <- Emballonuoidea_char_dat$litter
names(Emballonuoidea_states_litter) = Emballonuoidea_char_dat$species
#
Noctilionoidea_states_litter <- Noctilionoidea_char_dat$litter
names(Noctilionoidea_states_litter) = Noctilionoidea_char_dat$species
#
Pteropodoidea_states_litter <- Pteropodoidea_char_dat$litter
names(Pteropodoidea_states_litter) = Pteropodoidea_char_dat$species
#
Rhinolophoidea_states_litter <- Rhinolophoidea_char_dat$litter
names(Rhinolophoidea_states_litter) = Rhinolophoidea_char_dat$species
#
Vespertilionoidea_states_litter <- Vespertilionoidea_char_dat$litter
names(Vespertilionoidea_states_litter) = Vespertilionoidea_char_dat$species

# create data frame
bats_df = data.frame(clade = NA, tip_count = NA, root_age = NA)
#
bats_df[1,] = c("Emballonuoidea_plant",length(Emballonuoidea_mcc$tip.label),round(max(node.depth.edgelength(Emballonuoidea_mcc)),2))
bats_df[2,] = c("Noctilionoidea_plant",length(Noctilionoidea_mcc$tip.label),round(max(node.depth.edgelength(Noctilionoidea_mcc)),2))

bats_df[3,] = c("Emballonuoidea_litter",length(Emballonuoidea_mcc$tip.label),round(max(node.depth.edgelength(Emballonuoidea_mcc)),2))
bats_df[4,] = c("Noctilionoidea_litter",length(Noctilionoidea_mcc$tip.label),round(max(node.depth.edgelength(Noctilionoidea_mcc)),2))
bats_df[5,] = c("Pteropodoidea_litter",length(Pteropodoidea_mcc$tip.label),round(max(node.depth.edgelength(Pteropodoidea_mcc)),2))
bats_df[6,] = c("Rhinolophoidea_litter",length(Rhinolophoidea_mcc$tip.label),round(max(node.depth.edgelength(Rhinolophoidea_mcc)),2))
bats_df[7,] = c("Vespertilionoidea_litter",length(Vespertilionoidea_mcc$tip.label),round(max(node.depth.edgelength(Vespertilionoidea_mcc)),2))

# create bisse class for each clade 
# note: Pteropodoidea, Rhinolophoidea and Vespertilionoidea dont have tip variation so I dont use them for plant character 
bisse_Emballonuoidea_plant     = make.bisse(Emballonuoidea_mcc,Emballonuoidea_states_plant)
bisse_Noctilionoidea_plant     = make.bisse(Noctilionoidea_mcc,Noctilionoidea_states_plant)

bisse_Emballonuoidea_litter    = make.bisse(Emballonuoidea_mcc,Emballonuoidea_states_litter)
bisse_Noctilionoidea_litter    = make.bisse(Noctilionoidea_mcc,Noctilionoidea_states_litter)
bisse_Pteropodoidea_litter     = make.bisse(Pteropodoidea_mcc,Pteropodoidea_states_litter)
bisse_Rhinolophoidea_litter    = make.bisse(Rhinolophoidea_mcc,Rhinolophoidea_states_litter)
bisse_Vespertilionoidea_litter = make.bisse(Vespertilionoidea_mcc,Vespertilionoidea_states_litter)

# fit bisse 
# draw initial starting point
lambda_0 = runif(1,0,0.5)
lambda_1 = runif(1,0,0.5)
mu_0     = runif(1,0,lambda_0)
mu_1     = runif(1,0,lambda_1)
q_01     = runif(1,0,0.5)
q_10     = runif(1,0,0.5)

params_Emballonuoidea_plant = c(lambda_0,lambda_1,
                                mu_0,mu_1,
                                q_01,q_10)

params_Noctilionoidea_plant  <-
  params_Emballonuoidea_litter  <- params_Noctilionoidea_litter <- 
  params_Pteropodoidea_litter   <- params_Rhinolophoidea_litter <- 
  params_Vespertilionoidea_litter <- params_Emballonuoidea_plant 



# create the likelihood functions
best_fit_Emballonuoidea_plant = find.mle(bisse_Emballonuoidea_plant,params_Emballonuoidea_plant,"subplex")
best_lik_Emballonuoidea_plant = best_fit_Emballonuoidea_plant$lnLik
params_Emballonuoidea_plant   = best_fit_Emballonuoidea_plant$par # update the param values

best_fit_Noctilionoidea_plant = find.mle(bisse_Noctilionoidea_plant,params_Noctilionoidea_plant ,"subplex")
best_lik_Noctilionoidea_plant = best_fit_Noctilionoidea_plant$lnLik
params_Noctilionoidea_plant   = best_fit_Noctilionoidea_plant$par # update the param values

best_fit_Emballonuoidea_litter = find.mle(bisse_Emballonuoidea_litter,params_Emballonuoidea_litter,"subplex")
best_lik_Emballonuoidea_litter = best_fit_Emballonuoidea_litter$lnLik
params_Emballonuoidea_litter   = best_fit_Emballonuoidea_litter$par # update the param values

best_fit_Noctilionoidea_litter = find.mle(bisse_Noctilionoidea_litter,params_Noctilionoidea_litter ,"subplex")
best_lik_Noctilionoidea_litter = best_fit_Noctilionoidea_litter$lnLik
params_Noctilionoidea_litter   = best_fit_Noctilionoidea_litter$par # update the param values

best_fit_Pteropodoidea_litter = find.mle(bisse_Pteropodoidea_litter,params_Pteropodoidea_litter,"subplex")
best_lik_Pteropodoidea_litter = best_fit_Pteropodoidea_litter$lnLik
params_Pteropodoidea_litter   = best_fit_Pteropodoidea_litter$par # update the param values

best_fit_Rhinolophoidea_litter = find.mle(bisse_Rhinolophoidea_litter,params_Rhinolophoidea_litter ,"subplex")
best_lik_Rhinolophoidea_litter = best_fit_Rhinolophoidea_litter$lnLik
params_Rhinolophoidea_litter   = best_fit_Rhinolophoidea_litter$par # update the param values

best_fit_Vespertilionoidea_litter = find.mle(bisse_Vespertilionoidea_litter,params_Vespertilionoidea_litter ,"subplex")
best_lik_Vespertilionoidea_litter = best_fit_Vespertilionoidea_litter$lnLik
params_Vespertilionoidea_litter   = best_fit_Vespertilionoidea_litter$par # update the param values

for (j in 1:num_starts){
  print(paste0("running iteration ",j))
  # draw new starting point for both clades for each character
  new_params_Emballonuoidea_plant = jitters_func(params_Emballonuoidea_plant)
  new_params_Noctilionoidea_plant = jitters_func(params_Noctilionoidea_plant)
  #
  new_params_Emballonuoidea_litter = jitters_func(params_Emballonuoidea_litter)
  new_params_Noctilionoidea_litter = jitters_func(params_Noctilionoidea_litter)
  new_params_Pteropodoidea_litter = jitters_func(params_Pteropodoidea_litter)
  new_params_Rhinolophoidea_litter = jitters_func(params_Rhinolophoidea_litter)
  new_params_Vespertilionoidea_litter = jitters_func(params_Vespertilionoidea_litter)
  # Compute the new likelihood
  fit_Emballonuoidea_plant = find.mle(bisse_Emballonuoidea_plant,new_params_Emballonuoidea_plant,"subplex")
  fit_Noctilionoidea_plant = find.mle(bisse_Noctilionoidea_plant,new_params_Noctilionoidea_plant,"subplex")
  #
  fit_Emballonuoidea_litter = find.mle(bisse_Emballonuoidea_litter,new_params_Emballonuoidea_litter,"subplex")
  fit_Noctilionoidea_litter = find.mle(bisse_Noctilionoidea_litter,new_params_Noctilionoidea_litter,"subplex")
  fit_Pteropodoidea_litter = find.mle(bisse_Pteropodoidea_litter,new_params_Pteropodoidea_litter,"subplex")
  fit_Rhinolophoidea_litter = find.mle(bisse_Rhinolophoidea_litter,new_params_Rhinolophoidea_litter,"subplex")
  fit_Vespertilionoidea_litter = find.mle(bisse_Vespertilionoidea_litter,new_params_Vespertilionoidea_litter,"subplex")
  #
  if (fit_Emballonuoidea_plant$lnLik > best_lik_Emballonuoidea_plant){
    best_lik_Emballonuoidea_plant = fit_Emballonuoidea_plant$lnLik
    best_fit_Emballonuoidea_plant = fit_Emballonuoidea_plant
    print("accept new param for Emballonuoidea")
    params_Emballonuoidea_plant  = best_fit_Emballonuoidea_plant$par #accept the new values and use it as the new starting values
  }
  if (fit_Noctilionoidea_plant$lnLik > best_lik_Noctilionoidea_plant){
    best_lik_Noctilionoidea_plant = fit_Noctilionoidea_plant$lnLik
    best_fit_Noctilionoidea_plant = fit_Noctilionoidea_plant
    print("accept new param for Noctilionoidea")
    params_Noctilionoidea_plant = best_fit_Noctilionoidea_plant$par
  }
  if (fit_Emballonuoidea_litter$lnLik > best_lik_Emballonuoidea_litter){
    best_lik_Emballonuoidea_litter = fit_Emballonuoidea_litter$lnLik
    best_fit_Emballonuoidea_litter = fit_Emballonuoidea_litter
    print("accept new param for Emballonuoidea litter")
    params_Emballonuoidea_litter = best_fit_Emballonuoidea_litter$par
  }
  if (fit_Noctilionoidea_litter$lnLik > best_lik_Noctilionoidea_litter){
    best_lik_Noctilionoidea_litter = fit_Noctilionoidea_litter$lnLik
    best_fit_Noctilionoidea_litter = fit_Noctilionoidea_litter
    print("accept new param for Noctilionoidea litter")
    params_Noctilionoidea_litter = best_fit_Noctilionoidea_litter$par
  }
  if (fit_Pteropodoidea_litter$lnLik > best_lik_Pteropodoidea_litter){
    best_lik_Pteropodoidea_litter = fit_Pteropodoidea_litter$lnLik
    best_fit_Pteropodoidea_litter = fit_Pteropodoidea_litter
    print("accept new param for Pteropodoidea litter")
    params_Pteropodoidea_litter = best_fit_Pteropodoidea_litter$par
  }
  if (fit_Rhinolophoidea_litter$lnLik > best_lik_Rhinolophoidea_litter){
    best_lik_Rhinolophoidea_litter = fit_Rhinolophoidea_litter$lnLik
    best_fit_Rhinolophoidea_litter = fit_Rhinolophoidea_litter
    print("accept new param for Rhinolophoidea litter")
    params_Rhinolophoidea_litter = best_fit_Rhinolophoidea_litter$par
  }
  if (fit_Vespertilionoidea_litter$lnLik > best_lik_Vespertilionoidea_litter){
    best_lik_Vespertilionoidea_litter = fit_Vespertilionoidea_litter$lnLik
    best_fit_Vespertilionoidea_litter = fit_Vespertilionoidea_litter
    print("accept new param for Vespertilionoidea litter")
    params_Vespertilionoidea_litter = best_fit_Vespertilionoidea_litter$par
  }
}

# add the fitted parameter
bats_df$lambda_0 = c(params_Emballonuoidea_plant[1],params_Noctilionoidea_plant[1],
                     params_Emballonuoidea_litter[1],params_Noctilionoidea_litter[1],
                     params_Pteropodoidea_litter[1],params_Rhinolophoidea_litter[1],
                     params_Vespertilionoidea_litter[1])

bats_df$lambda_1 = c(params_Emballonuoidea_plant[2],params_Noctilionoidea_plant[2],
                     params_Emballonuoidea_litter[2],params_Noctilionoidea_litter[2],
                     params_Pteropodoidea_litter[2],params_Rhinolophoidea_litter[2],
                     params_Vespertilionoidea_litter[2])

bats_df$mu_0 = c(params_Emballonuoidea_plant[3],params_Noctilionoidea_plant[3],
                 params_Emballonuoidea_litter[3],params_Noctilionoidea_litter[3],
                 params_Pteropodoidea_litter[3],params_Rhinolophoidea_litter[3],
                 params_Vespertilionoidea_litter[3])

bats_df$mu_1 = c(params_Emballonuoidea_plant[4],params_Noctilionoidea_plant[4],
                 params_Emballonuoidea_litter[4],params_Noctilionoidea_litter[4],
                 params_Pteropodoidea_litter[4],params_Rhinolophoidea_litter[4],
                 params_Vespertilionoidea_litter[4])

bats_df$q_01 = c(params_Emballonuoidea_plant[5],params_Noctilionoidea_plant[5],
                 params_Emballonuoidea_litter[5],params_Noctilionoidea_litter[5],
                 params_Pteropodoidea_litter[5],params_Rhinolophoidea_litter[5],
                 params_Vespertilionoidea_litter[5])

bats_df$q_10 = c(params_Emballonuoidea_plant[6],params_Noctilionoidea_plant[6],
                 params_Emballonuoidea_litter[6],params_Noctilionoidea_litter[6],
                 params_Pteropodoidea_litter[6],params_Rhinolophoidea_litter[6],
                 params_Vespertilionoidea_litter[6])

# theoretical stationary dist based on best fitting parameters
bats_df$statio_0 = c(compute_statio(params_Emballonuoidea_plant),compute_statio(params_Noctilionoidea_plant),
                     compute_statio(params_Emballonuoidea_litter),compute_statio(params_Noctilionoidea_litter),
                     compute_statio(params_Pteropodoidea_litter),compute_statio(params_Rhinolophoidea_litter),
                     compute_statio(params_Vespertilionoidea_litter))

bats_df$statio_1 = 1-bats_df$statio_0

# observed frequency
bats_df$emp_0 = c(length(which(Emballonuoidea_char_dat$plant==0))/length(Emballonuoidea_char_dat$plant),
                  length(which(Noctilionoidea_char_dat$plant==0))/length(Noctilionoidea_char_dat$plant),
                  length(which(Emballonuoidea_char_dat$litter==0))/length(Emballonuoidea_char_dat$litter),
                  length(which(Noctilionoidea_char_dat$litter==0))/length(Noctilionoidea_char_dat$litter),
                  length(which(Pteropodoidea_char_dat$litter==0))/length(Pteropodoidea_char_dat$litter),
                  length(which(Rhinolophoidea_char_dat$litter==0))/length(Rhinolophoidea_char_dat$litter),
                  length(which(Vespertilionoidea_char_dat$litter==0))/length(Vespertilionoidea_char_dat$litter))

bats_df$emp_1 = 1-bats_df$emp_0

# compute Euclidean distance between obs freq and statio freq across clades
bats_df$distance = sqrt((bats_df$emp_0-bats_df$statio_0)^2+(bats_df$emp_1-bats_df$statio_1)^2)

# save output
write.table(bats_df,paste0(out_fp,"bats_summary.csv"),sep = ";",row.names = F)

#then draw 3D plot