#### LOAD LIBRARIES ####
library(diversitree)
library(ape)
library(phangorn)
library(R.utils)

#### FILE SETTING ####
# fp     = "/Users/albertsoewongsono/Documents/Code\ Testing/rate_pattern_diffusion_SSE/"
fp     = "/storage/albert/rate_pattern_diffusion_SSE/"
in_fp  = paste0(fp,"data/Empirical/squamate_phylogeny/")
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


# load squamate clades
gekkota_mcc = read.tree(paste0(in_fp,"clade_data/","gekkota.tre"))
alethinophidia_mcc = read.tree(paste0(in_fp,"clade_data/","alethinophidia.tre"))
acrodontia_mcc = read.tree(paste0(in_fp,"clade_data/","acrodontia.tre"))
laterata_mcc = read.tree(paste0(in_fp,"clade_data/","laterata.tre"))
scolecophidia_mcc = read.tree(paste0(in_fp,"clade_data/","scolecophidia.tre"))
anguimoprha_mcc = read.tree(paste0(in_fp,"clade_data/","anguimoprha.tre"))
pleurodonta_mcc = read.tree(paste0(in_fp,"clade_data/","pleurodonta.tre"))
dibamidae_mcc = read.tree(paste0(in_fp,"clade_data/","dibamidae.tre"))
scincoidea_mcc = read.tree(paste0(in_fp,"clade_data/","scincoidea.tre"))
amphisbaenia_mcc = read.tree(paste0(in_fp,"clade_data/","amphisbaenia.tre"))

# load character data 
# diet: 0 (no_plant) & 1 (plant)
gekkota_char_dat = read.table(paste0(in_fp,"trait_data/","gekkota_diet.csv"),sep=";",header = T)
alethinophidia_char_dat = read.table(paste0(in_fp,"trait_data/","alethinophidia_diet.csv"),sep=";",header = T)
acrodontia_char_dat = read.table(paste0(in_fp,"trait_data/","acrodontia_diet.csv"),sep=";",header = T)
laterata_char_dat = read.table(paste0(in_fp,"trait_data/","laterata_diet.csv"),sep=";",header = T)
scolecophidia_char_dat = read.table(paste0(in_fp,"trait_data/","scolecophidia_diet.csv"),sep=";",header = T)
anguimoprha_char_dat = read.table(paste0(in_fp,"trait_data/","anguimoprha_diet.csv"),sep=";",header = T)
pleurodonta_char_dat = read.table(paste0(in_fp,"trait_data/","pleurodonta_diet.csv"),sep=";",header = T)
dibamidae_char_dat = read.table(paste0(in_fp,"trait_data/","dibamidae_diet.csv"),sep=";",header = T)
scincoidea_char_dat = read.table(paste0(in_fp,"trait_data/","scincoidea_diet.csv"),sep=";",header = T)
amphisbaenia_char_dat = read.table(paste0(in_fp,"trait_data/","amphisbaenia_diet.csv"),sep=";",header = T)

gekkota_char_dat$Diet[gekkota_char_dat$Diet=="no_plant"] <- 0
gekkota_char_dat$Diet[gekkota_char_dat$Diet=="plant"] <- 1

alethinophidia_char_dat$Diet[alethinophidia_char_dat$Diet=="no_plant"] <- 0
alethinophidia_char_dat$Diet[alethinophidia_char_dat$Diet=="plant"] <- 1

acrodontia_char_dat$Diet[acrodontia_char_dat$Diet=="no_plant"] <- 0
acrodontia_char_dat$Diet[acrodontia_char_dat$Diet=="plant"] <- 1

laterata_char_dat$Diet[laterata_char_dat$Diet=="no_plant"] <- 0
laterata_char_dat$Diet[laterata_char_dat$Diet=="plant"] <- 1

scolecophidia_char_dat$Diet[scolecophidia_char_dat$Diet=="no_plant"] <- 0
scolecophidia_char_dat$Diet[scolecophidia_char_dat$Diet=="plant"] <- 1

anguimoprha_char_dat$Diet[anguimoprha_char_dat$Diet=="no_plant"] <- 0
anguimoprha_char_dat$Diet[anguimoprha_char_dat$Diet=="plant"] <- 1

pleurodonta_char_dat$Diet[pleurodonta_char_dat$Diet=="no_plant"] <- 0
pleurodonta_char_dat$Diet[pleurodonta_char_dat$Diet=="plant"] <- 1

dibamidae_char_dat$Diet[dibamidae_char_dat$Diet=="no_plant"] <- 0
dibamidae_char_dat$Diet[dibamidae_char_dat$Diet=="plant"] <- 1

scincoidea_char_dat$Diet[scincoidea_char_dat$Diet=="no_plant"] <- 0
scincoidea_char_dat$Diet[scincoidea_char_dat$Diet=="plant"] <- 1

amphisbaenia_char_dat$Diet[amphisbaenia_char_dat$Diet=="no_plant"] <- 0
amphisbaenia_char_dat$Diet[amphisbaenia_char_dat$Diet=="plant"] <- 1

# remove unused row 
gekkota_char_dat = gekkota_char_dat[-2]
alethinophidia_char_dat = alethinophidia_char_dat[-2]
acrodontia_char_dat = acrodontia_char_dat[-2]
laterata_char_dat = laterata_char_dat[-2]
scolecophidia_char_dat = scolecophidia_char_dat[-2]
anguimoprha_char_dat = anguimoprha_char_dat[-2]
pleurodonta_char_dat = pleurodonta_char_dat[-2]
dibamidae_char_dat = dibamidae_char_dat[-2]
scincoidea_char_dat = scincoidea_char_dat[-2]
amphisbaenia_char_dat = amphisbaenia_char_dat[-2]

# change to numeric str for tip states
# commented ones mean no state variation observed at the tips
gekkota_char_dat$Diet = as.numeric(gekkota_char_dat$Diet)
# alethinophidia_char_dat$Diet = as.numeric(alethinophidia_char_dat$Diet)
acrodontia_char_dat$Diet = as.numeric(acrodontia_char_dat$Diet)
laterata_char_dat$Diet = as.numeric(laterata_char_dat$Diet)
# scolecophidia_char_dat$Diet = as.numeric(scolecophidia_char_dat$Diet)
anguimoprha_char_dat$Diet = as.numeric(anguimoprha_char_dat$Diet)
pleurodonta_char_dat$Diet = as.numeric(pleurodonta_char_dat$Diet)
# dibamidae_char_dat$Diet = as.numeric(dibamidae_char_dat$Diet)
scincoidea_char_dat$Diet = as.numeric(scincoidea_char_dat$Diet)
# amphisbaenia_char_dat$Diet = as.numeric(amphisbaenia_char_dat$Diet)

# name the rows using species name
rownames(gekkota_char_dat) = gekkota_char_dat$Species.name..Binomial.
rownames(acrodontia_char_dat) = acrodontia_char_dat$Species.name..Binomial.
rownames(laterata_char_dat) = laterata_char_dat$Species.name..Binomial.
rownames(anguimoprha_char_dat) = anguimoprha_char_dat$Species.name..Binomial.
rownames(pleurodonta_char_dat) = pleurodonta_char_dat$Species.name..Binomial.
rownames(scincoidea_char_dat) = scincoidea_char_dat$Species.name..Binomial.

# reorder rows in character data to match the order of species on trees 
gekkota_char_dat = gekkota_char_dat[gekkota_mcc$tip.label,]
acrodontia_char_dat = acrodontia_char_dat[acrodontia_mcc$tip.label,]
laterata_char_dat = laterata_char_dat[laterata_mcc$tip.label,]
anguimoprha_char_dat = anguimoprha_char_dat[anguimoprha_mcc$tip.label,]
pleurodonta_char_dat = pleurodonta_char_dat[pleurodonta_mcc$tip.label,]
scincoidea_char_dat = scincoidea_char_dat[scincoidea_mcc$tip.label,]

# create a vector of tip labels and their states for each character
##### CHARACTER: DIET
gekkota_states <- gekkota_char_dat$Diet
names(gekkota_states) = gekkota_char_dat$Species.name..Binomial.

acrodontia_states <- acrodontia_char_dat$Diet
names(acrodontia_states) = acrodontia_char_dat$Species.name..Binomial.

laterata_states <- laterata_char_dat$Diet
names(laterata_states) = laterata_char_dat$Species.name..Binomial.

anguimoprha_states <- anguimoprha_char_dat$Diet
names(anguimoprha_states) = anguimoprha_char_dat$Species.name..Binomial.

pleurodonta_states <- pleurodonta_char_dat$Diet
names(pleurodonta_states) = pleurodonta_char_dat$Species.name..Binomial.

scincoidea_states <- scincoidea_char_dat$Diet
names(scincoidea_states) = scincoidea_char_dat$Species.name..Binomial.

# create data frame
squamate_df = data.frame(clade = NA, tip_count = NA, root_age = NA)
#
squamate_df[1,] = c("gekkota",length(gekkota_mcc$tip.label),round(max(node.depth.edgelength(gekkota_mcc)),2))
squamate_df[2,] = c("acrodontia",length(acrodontia_mcc$tip.label),round(max(node.depth.edgelength(acrodontia_mcc)),2))
squamate_df[3,] = c("laterata",length(laterata_mcc$tip.label),round(max(node.depth.edgelength(laterata_mcc)),2))
squamate_df[4,] = c("anguimoprha",length(anguimoprha_mcc$tip.label),round(max(node.depth.edgelength(anguimoprha_mcc)),2))
squamate_df[5,] = c("pleurodonta",length(pleurodonta_mcc$tip.label),round(max(node.depth.edgelength(pleurodonta_mcc)),2))
squamate_df[6,] = c("scincoidea",length(scincoidea_mcc$tip.label),round(max(node.depth.edgelength(scincoidea_mcc)),2))

# create bisse class for each clade 
bisse_gekkota                  = make.bisse(gekkota_mcc,gekkota_states)
bisse_acrodontia               = make.bisse(acrodontia_mcc,acrodontia_states)
bisse_laterata                 = make.bisse(laterata_mcc,laterata_states)
bisse_anguimoprha              = make.bisse(anguimoprha_mcc,anguimoprha_states)
bisse_pleurodonta              = make.bisse(pleurodonta_mcc,pleurodonta_states)
bisse_scincoidea               = make.bisse(scincoidea_mcc,scincoidea_states)

# fit bisse 
# draw initial starting point
lambda_0 = runif(1,0,0.5)
lambda_1 = runif(1,0,0.5)
mu_0     = runif(1,0,lambda_0)
mu_1     = runif(1,0,lambda_1)
q_01     = runif(1,0,0.5)
q_10     = runif(1,0,0.5)

params_gekkota         = c(lambda_0,lambda_1,
                         mu_0,mu_1,
                         q_01,q_10)

params_acrodontia  <-
  params_laterata  <- params_anguimoprha <- 
  params_pleurodonta   <- params_scincoidea <- 
  params_gekkota

# create the likelihood functions
best_fit_gekkota = safe_find_mle(bisse_gekkota,params_gekkota)
best_lik_gekkota = best_fit_gekkota$lnLik
params_gekkota = best_fit_gekkota$par # start with the optimized param after the intiial draw

best_fit_acrodontia = safe_find_mle(bisse_acrodontia,params_acrodontia)
best_lik_acrodontia = best_fit_acrodontia$lnLik
params_acrodontia = best_fit_acrodontia$par # start with the optimized param after the intiial draw

best_fit_laterata = safe_find_mle(bisse_laterata,params_laterata)
best_lik_laterata = best_fit_laterata$lnLik
params_laterata = best_fit_laterata$par # start with the optimized param after the intiial draw

best_fit_anguimoprha = safe_find_mle(bisse_anguimoprha,params_anguimoprha)
best_lik_anguimoprha = best_fit_anguimoprha$lnLik
params_anguimoprha = best_fit_anguimoprha$par # start with the optimized param after the intiial draw

best_fit_pleurodonta = safe_find_mle(bisse_pleurodonta,params_pleurodonta)
best_lik_pleurodonta = best_fit_pleurodonta$lnLik
params_pleurodonta = best_fit_pleurodonta$par # start with the optimized param after the intiial draw

best_fit_scincoidea = safe_find_mle(bisse_scincoidea,params_scincoidea)
best_lik_scincoidea = best_fit_scincoidea$lnLik
params_scincoidea = best_fit_scincoidea$par # start with the optimized param after the intiial draw


for (j in 1:num_starts){
  print(paste0("running iteration ",j))
  # draw new starting point for both clades for each character
  new_params_gekkota = jitters_func(params_gekkota)
  new_params_acrodontia = jitters_func(params_acrodontia)
  new_params_laterata = jitters_func(params_laterata)
  new_params_anguimoprha = jitters_func(params_anguimoprha)
  new_params_pleurodonta = jitters_func(params_pleurodonta)
  new_params_scincoidea = jitters_func(params_scincoidea)
  # Compute the new likelihood
  fit_gekkota = safe_find_mle(bisse_gekkota,new_params_gekkota)
  fit_acrodontia = safe_find_mle(bisse_acrodontia,new_params_acrodontia)
  fit_laterata = safe_find_mle(bisse_laterata,new_params_laterata)
  fit_anguimoprha = safe_find_mle(bisse_anguimoprha,new_params_anguimoprha)
  fit_pleurodonta = safe_find_mle(bisse_pleurodonta,new_params_pleurodonta)
  fit_scincoidea = safe_find_mle(bisse_scincoidea,new_params_scincoidea)
  #
  if (fit_gekkota$lnLik > best_lik_gekkota){
    best_lik_gekkota = fit_gekkota$lnLik
    best_fit_gekkota = fit_gekkota
    print("accept new param for gekkota")
    params_gekkota  = best_fit_gekkota$par #accept the new values and use it as the new starting values
  }
  #
  if (fit_acrodontia$lnLik > best_lik_acrodontia){
    best_lik_acrodontia = fit_acrodontia$lnLik
    best_fit_acrodontia = fit_acrodontia
    print("accept new param for acrodontia")
    params_acrodontia  = best_fit_acrodontia$par #accept the new values and use it as the new starting values
  }
  #
  if (fit_laterata$lnLik > best_lik_laterata){
    best_lik_laterata = fit_laterata$lnLik
    best_fit_laterata = fit_laterata
    print("accept new param for laterata")
    params_laterata  = best_fit_laterata$par #accept the new values and use it as the new starting values
  }
  #
  if (fit_anguimoprha$lnLik > best_lik_anguimoprha){
    best_lik_anguimoprha = fit_anguimoprha$lnLik
    best_fit_anguimoprha = fit_anguimoprha
    print("accept new param for anguimoprha")
    params_anguimoprha  = best_fit_anguimoprha$par #accept the new values and use it as the new starting values
  }
  #
  if (fit_pleurodonta$lnLik > best_lik_pleurodonta){
    best_lik_pleurodonta = fit_pleurodonta$lnLik
    best_fit_pleurodonta = fit_pleurodonta
    print("accept new param for pleurodonta")
    params_pleurodonta  = best_fit_pleurodonta$par #accept the new values and use it as the new starting values
  }
  #
  if (fit_scincoidea$lnLik > best_lik_scincoidea){
    best_lik_scincoidea = fit_scincoidea$lnLik
    best_fit_scincoidea = fit_scincoidea
    print("accept new param for scincoidea")
    params_scincoidea  = best_fit_scincoidea$par #accept the new values and use it as the new starting values
  }
}

# add the fitted parameter
squamate_df$lambda_0 = c(params_gekkota[1],params_acrodontia[1],
                      params_laterata[1],params_anguimoprha[1],
                      params_pleurodonta[1],params_scincoidea[1])

squamate_df$lambda_1 = c(params_gekkota[2],params_acrodontia[2],
                         params_laterata[2],params_anguimoprha[2],
                         params_pleurodonta[2],params_scincoidea[2])

squamate_df$mu_0 = c(params_gekkota[3],params_acrodontia[3],
                         params_laterata[3],params_anguimoprha[3],
                         params_pleurodonta[3],params_scincoidea[3])

squamate_df$mu_1 = c(params_gekkota[4],params_acrodontia[4],
                     params_laterata[4],params_anguimoprha[4],
                     params_pleurodonta[4],params_scincoidea[4])

squamate_df$q_01 = c(params_gekkota[5],params_acrodontia[5],
                     params_laterata[5],params_anguimoprha[5],
                     params_pleurodonta[5],params_scincoidea[5])

squamate_df$q_10 = c(params_gekkota[6],params_acrodontia[6],
                     params_laterata[6],params_anguimoprha[6],
                     params_pleurodonta[6],params_scincoidea[6])

# theoretical stationary dist based on best fitting parameters
squamate_df$statio_0 = c(compute_statio(params_gekkota),
                      compute_statio(params_acrodontia),
                      compute_statio(params_laterata),
                      compute_statio(params_anguimoprha),
                      compute_statio(params_pleurodonta),
                      compute_statio(params_scincoidea))

squamate_df$statio_1 = 1-squamate_df$statio_0

# observed frequency
squamate_df$emp_0 = c(length(which(gekkota_char_dat$Diet==0))/length(gekkota_char_dat$Diet),
                   length(which(acrodontia_char_dat$Diet==0))/length(acrodontia_char_dat$Diet),
                   length(which(laterata_char_dat$Diet==0))/length(laterata_char_dat$Diet),
                   length(which(anguimoprha_char_dat$Diet==0))/length(anguimoprha_char_dat$Diet),
                   length(which(pleurodonta_char_dat$Diet==0))/length(pleurodonta_char_dat$Diet),
                   length(which(scincoidea_char_dat$Diet==0))/length(scincoidea_char_dat$Diet))

squamate_df$emp_1 = 1-squamate_df$emp_0

# compute Euclidean distance between obs freq and statio freq across clades
squamate_df$distance = sqrt((squamate_df$emp_0-squamate_df$statio_0)^2+(squamate_df$emp_1-squamate_df$statio_1)^2)

# save output
write.table(squamate_df,paste0(out_fp,"squamate_summary.csv"),sep = ";",row.names = F)

#then draw 3D plot