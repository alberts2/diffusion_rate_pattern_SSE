# NOTE: don't use extract.clade because it will also include species that share common ancestor but does not belong
# to our group of interest 

#### LOAD LIBRARIES ####
library(diversitree)
library(ape)
library(phangorn)
library(R.utils)

#### FILE SETTING ####
# fp     = "/Users/albertsoewongsono/Documents/Code\ Testing/rate_pattern_diffusion_SSE/"
fp     = "/storage/albert/rate_pattern_diffusion_SSE/"
in_fp  = paste0(fp,"data/Empirical/grass_phylogeny/")
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

# add failsafe if happen to draw bad parameters as starting value to compute the likelihood
safe_find_mle <- function(model, params, method = "subplex",
                          timeout = 100, max_attempts = 50) {
  
  for (i in seq_len(max_attempts)) {
    cat("Attempt", i, "\n")
    
    result <- tryCatch({
      withTimeout({
        nan_flag <- FALSE
        
        fit <- withCallingHandlers({
          diversitree::find.mle(model, params, method = method,
                                control = list(reltol = 1e-6, maxit = 1000))
        }, warning = function(w) {
          if (grepl("NaNs produced", w$message)) nan_flag <<- TRUE
          invokeRestart("muffleWarning")
        })
        
        # If NaN occurred, treat as failure
        if (nan_flag || any(!is.finite(fit$lnLik))) stop("NaN or invalid likelihood")
        
        fit
      }, timeout = timeout, onTimeout = "error")
    }, error = function(e) {
      cat("Attempt", i, "failed, timed out, or NaN produced:", e$message, "\n")
      NULL
    })
    
    if (!is.null(result)) {
      cat("Success on attempt", i, "\n")
      return(result)
    } else {
      # redraw parameters
      params <- jitters_func(params)
      cat("Redrawing params:", params, "\n")
    }
    
  }
  stop("Failed to find MLE after max attempts")
}


# load tree 
### GRASS MCC phylogeny (3595 species)
grass                       = read.tree("Poaceae_hyp1_allC4_maxcred.phy")
# load each clade nexus file 
Andropogoneae               = read.nexus.data("nexus_files/Andropogoneae.nex")
Aristidoideae               = read.nexus.data("nexus_files/Aristidoideae.nex")
Arundinoideae_Micrairoideae = read.nexus.data("nexus_files/Arundinoideae_Micrairoideae.nex")
Bambusiodeae                = read.nexus.data("nexus_files/Bambusiodeae.nex")
Centotheaceae               = read.nexus.data("nexus_files/Centotheaceae.nex")
Chloridoideae               = read.nexus.data("nexus_files/Chloridoideae.nex")
Danthonioideae              = read.nexus.data("nexus_files/Danthonioideae.nex")
Ehrhartoideae               = read.nexus.data("nexus_files/Ehrhartoideae.nex")
Paniceae                    = read.nexus.data("nexus_files/Paniceae.nex")
Paspaleae                   = read.nexus.data("nexus_files/Paspaleae.nex")
Pooideae                    = read.nexus.data("nexus_files/Pooideae.nex")

# get species name from each clade
Andropogoneae_vec               = names(Andropogoneae)
Aristidoideae_vec               = names(Aristidoideae)
Arundinoideae_Micrairoideae_vec = names(Arundinoideae_Micrairoideae)
#
Bambusiodeae_vec                = names(Bambusiodeae)
# remove taxa that are not in the  grass tree
Bambusiodeae_vec                = Bambusiodeae_vec[Bambusiodeae_vec != "Bonia_saxatilis"]
#
Centotheaceae_vec               = names(Centotheaceae)
Chloridoideae_vec               = names(Chloridoideae)
Danthonioideae_vec              = names(Danthonioideae)
Ehrhartoideae_vec               = names(Ehrhartoideae)
Paniceae_vec                    = names(Paniceae)
Paspaleae_vec                   = names(Paspaleae)
#
Pooideae_vec                    = names(Pooideae)
# remove taxa that are not in the  grass tree
Pooideae_vec                    = Pooideae_vec[Pooideae_vec != "Aegilops_neglecta" & Pooideae_vec != "Calamagrostis_emodensis_Peterson__Saarela_18749"]

# get the clades
Andropogoneae_mcc                = drop.tip(grass,setdiff(grass$tip.label,Andropogoneae_vec))
Aristidoideae_mcc                = drop.tip(grass,setdiff(grass$tip.label,Aristidoideae_vec))
Arundinoideae_Micrairoideae_mcc  = drop.tip(grass,setdiff(grass$tip.label,Arundinoideae_Micrairoideae_vec))
Bambusiodeae_mcc                 = drop.tip(grass,setdiff(grass$tip.label,Bambusiodeae_vec))
Centotheaceae_mcc                = drop.tip(grass,setdiff(grass$tip.label,Centotheaceae_vec))
Chloridoideae_mcc                = drop.tip(grass,setdiff(grass$tip.label,Chloridoideae_vec))
Danthonioideae_mcc               = drop.tip(grass,setdiff(grass$tip.label,Danthonioideae_vec))
Ehrhartoideae_mcc                = drop.tip(grass,setdiff(grass$tip.label,Ehrhartoideae_vec))
Paniceae_mcc                     = drop.tip(grass,setdiff(grass$tip.label,Paniceae_vec))
Paspaleae_mcc                    = drop.tip(grass,setdiff(grass$tip.label,Paspaleae_vec))
Pooideae_mcc                     = drop.tip(grass,setdiff(grass$tip.label,Pooideae_vec))

# load bats character data 
# photosynthesis type       = 0 (C3), 1 (C4)
char_dat_grass           = read.table("C3_C4_coding.txt",sep="\t",header = T)
colnames(char_dat_grass) = c("species","tip_states")

# get the character data for each clade 
# also, ignore all those clades that do not have tip variation (see commented lines below)

# Andropogoneae_char_dat    = char_dat_grass[char_dat_grass$species %in% Andropogoneae_mcc$tip.label,]
Aristidoideae_char_dat    = char_dat_grass[char_dat_grass$species %in% Aristidoideae_mcc$tip.label,]
Arundinoideae_Micrairoideae_char_dat = char_dat_grass[char_dat_grass$species %in% Arundinoideae_Micrairoideae_mcc$tip.label,]
# Bambusiodeae_char_dat     = char_dat_grass[char_dat_grass$species %in% Bambusiodeae_mcc$tip.label,]
Centotheaceae_char_dat    = char_dat_grass[char_dat_grass$species %in% Centotheaceae_mcc$tip.label,]
Chloridoideae_char_dat    = char_dat_grass[char_dat_grass$species %in% Chloridoideae_mcc$tip.label,]
# Danthonioideae_char_dat   = char_dat_grass[char_dat_grass$species %in% Danthonioideae_mcc$tip.label,]
# Ehrhartoideae_char_dat    = char_dat_grass[char_dat_grass$species %in% Ehrhartoideae_mcc$tip.label,]
Paniceae_char_dat         = char_dat_grass[char_dat_grass$species %in% Paniceae_mcc$tip.label,]
Paspaleae_char_dat        = char_dat_grass[char_dat_grass$species %in% Paspaleae_mcc$tip.label,]
# Pooideae_char_dat         = char_dat_grass[char_dat_grass$species %in% Pooideae_mcc$tip.label,]

# name the rows using species name
rownames(Aristidoideae_char_dat)               = Aristidoideae_char_dat$species
rownames(Arundinoideae_Micrairoideae_char_dat) = Arundinoideae_Micrairoideae_char_dat$species
rownames(Centotheaceae_char_dat)               = Centotheaceae_char_dat$species
rownames(Chloridoideae_char_dat)               = Chloridoideae_char_dat$species
rownames(Paniceae_char_dat)                    = Paniceae_char_dat$species
rownames(Paspaleae_char_dat )                  = Paspaleae_char_dat $species

# reorder rows in character data to match the order of species on trees 
Aristidoideae_char_dat    = Aristidoideae_char_dat[Aristidoideae_mcc$tip.label,]
Arundinoideae_Micrairoideae_char_dat    = Arundinoideae_Micrairoideae_char_dat[Arundinoideae_Micrairoideae_mcc$tip.label,]
Centotheaceae_char_dat    = Centotheaceae_char_dat[Centotheaceae_mcc$tip.label,]
Chloridoideae_char_dat    = Chloridoideae_char_dat[Chloridoideae_mcc$tip.label,]
Paniceae_char_dat         = Paniceae_char_dat[Paniceae_mcc$tip.label,]
Paspaleae_char_dat        = Paspaleae_char_dat[Paspaleae_mcc$tip.label,]


# create a vector of tip labels and their states for each character
##### CHARACTER: PHOTOSYNTHESIS TYPES
Aristidoideae_states <- Aristidoideae_char_dat$tip_states
names(Aristidoideae_states) = Aristidoideae_char_dat$species

Arundinoideae_Micrairoideae_states <- Arundinoideae_Micrairoideae_char_dat$tip_states
names(Arundinoideae_Micrairoideae_states) = Arundinoideae_Micrairoideae_char_dat$species

Centotheaceae_states <- Centotheaceae_char_dat$tip_states
names(Centotheaceae_states) = Centotheaceae_char_dat$species

Chloridoideae_states <- Chloridoideae_char_dat$tip_states
names(Chloridoideae_states) = Chloridoideae_char_dat$species

Paniceae_states <- Paniceae_char_dat$tip_states
names(Paniceae_states) = Paniceae_char_dat$species

Paspaleae_states <- Paspaleae_char_dat$tip_states
names(Paspaleae_states) = Paspaleae_char_dat$species

# create data frame
grass_df = data.frame(clade = NA, tip_count = NA, root_age = NA)
#
grass_df[1,] = c("Aristidoideae",length(Aristidoideae_mcc$tip.label),round(max(node.depth.edgelength(Aristidoideae_mcc)),2))
grass_df[2,] = c("Arundinoideae_Micrairoideae",length(Arundinoideae_Micrairoideae_mcc$tip.label),round(max(node.depth.edgelength(Arundinoideae_Micrairoideae_mcc)),2))
grass_df[3,] = c("Centotheaceae",length(Centotheaceae_mcc$tip.label),round(max(node.depth.edgelength(Centotheaceae_mcc)),2))
grass_df[4,] = c("Chloridoideae",length(Chloridoideae_mcc$tip.label),round(max(node.depth.edgelength(Chloridoideae_mcc)),2))
grass_df[5,] = c("Paniceae",length(Paniceae_mcc$tip.label),round(max(node.depth.edgelength(Paniceae_mcc)),2))
grass_df[6,] = c("Paspaleae",length(Paspaleae_mcc$tip.label),round(max(node.depth.edgelength(Paspaleae_mcc)),2))

# create bisse class for each clade 
bisse_Aristidoideae                  = make.bisse(Aristidoideae_mcc,Aristidoideae_states)
bisse_Arundinoideae_Micrairoideae    = make.bisse(Arundinoideae_Micrairoideae_mcc,Arundinoideae_Micrairoideae_states)
bisse_Centotheaceae                  = make.bisse(Centotheaceae_mcc,Centotheaceae_states)
bisse_Chloridoideae                  = make.bisse(Chloridoideae_mcc,Chloridoideae_states)
bisse_Paniceae                       = make.bisse(Paniceae_mcc,Paniceae_states)
bisse_Paspaleae                      = make.bisse(Paspaleae_mcc,Paspaleae_states)

# fit bisse 
# draw initial starting point
lambda_0 = runif(1,0,0.5)
lambda_1 = runif(1,0,0.5)
mu_0     = runif(1,0,lambda_0)
mu_1     = runif(1,0,lambda_1)
q_01     = runif(1,0,0.5)
q_10     = runif(1,0,0.5)

params_Aristidoideae = c(lambda_0,lambda_1,
                         mu_0,mu_1,
                         q_01,q_10)

params_Arundinoideae_Micrairoideae  <-
  params_Centotheaceae  <- params_Chloridoideae <- 
  params_Paniceae   <- params_Paspaleae <- 
  params_Aristidoideae

# create the likelihood functions
best_fit_Aristidoideae = safe_find_mle(bisse_Aristidoideae,params_Aristidoideae)
best_lik_Aristidoideae = best_fit_Aristidoideae$lnLik
params_Aristidoideae = best_fit_Aristidoideae$par # start with the optimized param after the intiial draw

best_fit_Arundinoideae_Micrairoideae = safe_find_mle(bisse_Arundinoideae_Micrairoideae,params_Arundinoideae_Micrairoideae)
best_lik_Arundinoideae_Micrairoideae = best_fit_Arundinoideae_Micrairoideae$lnLik
params_Arundinoideae_Micrairoideae = best_fit_Arundinoideae_Micrairoideae$par

best_fit_Centotheaceae = safe_find_mle(bisse_Centotheaceae,params_Centotheaceae)
best_lik_Centotheaceae = best_fit_Centotheaceae$lnLik
params_Centotheaceae = best_fit_Centotheaceae$par

best_fit_Chloridoideae = safe_find_mle(bisse_Chloridoideae,params_Chloridoideae)
best_lik_Chloridoideae = best_fit_Chloridoideae$lnLik
params_Chloridoideae = best_fit_Chloridoideae$par 

best_fit_Paniceae = safe_find_mle(bisse_Paniceae,params_Paniceae)
best_lik_Paniceae = best_fit_Paniceae$lnLik
params_Paniceae = best_fit_Paniceae$par
  
best_fit_Paspaleae = safe_find_mle(bisse_Paspaleae,params_Paspaleae)
best_lik_Paspaleae = best_fit_Paspaleae$lnLik
params_Paspaleae = best_fit_Paspaleae$par # start with param values that actually work

for (j in 1:num_starts){
  print(paste0("running iteration ",j))
  # draw new starting point for both clades for each character
  new_params_Aristidoideae = jitters_func(params_Aristidoideae)
  new_params_Arundinoideae_Micrairoideae = jitters_func(params_Arundinoideae_Micrairoideae)
  new_params_Centotheaceae = jitters_func(params_Centotheaceae)
  new_params_Chloridoideae = jitters_func(params_Chloridoideae)
  new_params_Paniceae = jitters_func(params_Paniceae)
  new_params_Paspaleae = jitters_func(params_Paspaleae)
  # Compute the new likelihood
  fit_Aristidoideae = safe_find_mle(bisse_Aristidoideae,new_params_Aristidoideae)
  fit_Arundinoideae_Micrairoideae = safe_find_mle(bisse_Arundinoideae_Micrairoideae,new_params_Arundinoideae_Micrairoideae)
  fit_Centotheaceae = safe_find_mle(bisse_Centotheaceae,new_params_Centotheaceae)
  fit_Chloridoideae = safe_find_mle(bisse_Chloridoideae,new_params_Chloridoideae)
  fit_Paniceae = safe_find_mle(bisse_Paniceae,new_params_Paniceae)
  fit_Paspaleae = safe_find_mle(bisse_Paspaleae,new_params_Paspaleae)
  #
  if (fit_Aristidoideae$lnLik > best_lik_Aristidoideae){
    best_lik_Aristidoideae = fit_Aristidoideae$lnLik
    best_fit_Aristidoideae = fit_Aristidoideae
    print("accept new param for Aristidoideae")
    params_Aristidoideae  = best_fit_Aristidoideae$par #accept the new values and use it as the new starting values
  }
  #
  if (fit_Arundinoideae_Micrairoideae$lnLik > best_lik_Arundinoideae_Micrairoideae){
    best_lik_Arundinoideae_Micrairoideae = fit_Arundinoideae_Micrairoideae$lnLik
    best_fit_Arundinoideae_Micrairoideae = fit_Arundinoideae_Micrairoideae
    print("accept new param for Arundinoideae_Micrairoideae")
    params_Arundinoideae_Micrairoideae  = best_fit_Arundinoideae_Micrairoideae$par #accept the new values and use it as the new starting values
  }
  #
  if (fit_Centotheaceae$lnLik > best_lik_Centotheaceae){
    best_lik_Centotheaceae = fit_Centotheaceae$lnLik
    best_fit_Centotheaceae = fit_Centotheaceae
    print("accept new param for Centotheaceae")
    params_Centotheaceae  = best_fit_Centotheaceae$par #accept the new values and use it as the new starting values
  }
  #
  if (fit_Chloridoideae$lnLik > best_lik_Chloridoideae){
    best_lik_Chloridoideae = fit_Chloridoideae$lnLik
    best_fit_Chloridoideae = fit_Chloridoideae
    print("accept new param for Chloridoideae")
    params_Chloridoideae  = best_fit_Chloridoideae$par #accept the new values and use it as the new starting values
  }
  #
  if (fit_Paniceae$lnLik > best_lik_Paniceae){
    best_lik_Paniceae = fit_Paniceae$lnLik
    best_fit_Paniceae = fit_Paniceae
    print("accept new param for Paniceae")
    params_Paniceae  = best_fit_Paniceae$par #accept the new values and use it as the new starting values
  }
  #
  if (fit_Paspaleae$lnLik > best_lik_Paspaleae){
    best_lik_Paspaleae = fit_Paspaleae$lnLik
    best_fit_Paspaleae = fit_Paspaleae
    print("accept new param for Paspaleae")
    params_Paspaleae  = best_fit_Paspaleae$par #accept the new values and use it as the new starting values
  }
}

# add the fitted parameter
grass_df$lambda_0 = c(params_Aristidoideae[1],params_Arundinoideae_Micrairoideae[1],
                     params_Centotheaceae[1],params_Chloridoideae[1],
                     params_Paniceae[1],params_Paspaleae[1])

grass_df$lambda_1 = c(params_Aristidoideae[2],params_Arundinoideae_Micrairoideae[2],
                      params_Centotheaceae[2],params_Chloridoideae[2],
                      params_Paniceae[2],params_Paspaleae[2])

grass_df$mu_0 = c(params_Aristidoideae[3],params_Arundinoideae_Micrairoideae[3],
                      params_Centotheaceae[3],params_Chloridoideae[3],
                      params_Paniceae[3],params_Paspaleae[3])

grass_df$mu_1 = c(params_Aristidoideae[4],params_Arundinoideae_Micrairoideae[4],
                  params_Centotheaceae[4],params_Chloridoideae[4],
                  params_Paniceae[4],params_Paspaleae[4])

grass_df$q_01 = c(params_Aristidoideae[5],params_Arundinoideae_Micrairoideae[5],
                  params_Centotheaceae[5],params_Chloridoideae[5],
                  params_Paniceae[5],params_Paspaleae[5])

grass_df$q_10 = c(params_Aristidoideae[6],params_Arundinoideae_Micrairoideae[6],
                  params_Centotheaceae[6],params_Chloridoideae[6],
                  params_Paniceae[6],params_Paspaleae[6])

# theoretical stationary dist based on best fitting parameters
grass_df$statio_0 = c(compute_statio(params_Aristidoideae),
                      compute_statio(params_Arundinoideae_Micrairoideae),
                      compute_statio(params_Centotheaceae),
                      compute_statio(params_Chloridoideae),
                      compute_statio(params_Paniceae),
                      compute_statio(params_Paspaleae))

grass_df$statio_1 = 1-grass_df$statio_0

# observed frequency
grass_df$emp_0 = c(length(which(Aristidoideae_char_dat$tip_states==0))/length(Aristidoideae_char_dat$tip_states),
                   length(which(Arundinoideae_Micrairoideae_char_dat$tip_states==0))/length(Arundinoideae_Micrairoideae_char_dat$tip_states),
                   length(which(Centotheaceae_char_dat$tip_states==0))/length(Centotheaceae_char_dat$tip_states),
                   length(which(Chloridoideae_char_dat$tip_states==0))/length(Chloridoideae_char_dat$tip_states),
                   length(which(Paniceae_char_dat$tip_states==0))/length(Paniceae_char_dat$tip_states),
                   length(which(Paspaleae_char_dat$tip_states==0))/length(Paspaleae_char_dat$tip_states))

grass_df$emp_1 = 1-grass_df$emp_0

# compute Euclidean distance between obs freq and statio freq across clades
grass_df$distance = sqrt((grass_df$emp_0-grass_df$statio_0)^2+(grass_df$emp_1-grass_df$statio_1)^2)

# save output
write.table(grass_df,paste0(out_fp,"grass_summary.csv"),sep = ";",row.names = F)

#then draw 3D plot