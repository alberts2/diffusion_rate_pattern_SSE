#### LOAD LIBRARIES ####
library(diversitree)
library(ape)
library(phangorn)
library(R.utils)

#### FILE SETTING ####
#fp     = "/Users/albertsoewongsono/Documents/Code\ Testing/rate_pattern_diffusion_SSE/"
fp     = "/storage/albert/rate_pattern_diffusion_SSE/"
in_fp  = paste0(fp,"data/Simulation/slow_rates_40percent_miss/combined_slow/")
out_fp = paste0(fp,"plot/")

setwd(fp)

num_starts = 10  #number of likelihood searches
# num_starts = 2  #number of likelihood searches
tree_sizes = c(25,50,100,200,400,800)  #different tree sizes in the dataset
# tree_sizes = c(25)  #different tree sizes in the dataset
#num_sim    = 10 #number of simulated data
num_sim    = 100 #number of simulated data
total_sim  = num_sim*length(tree_sizes)


# pars (lambda0,lambda1, mu0, mu1, q01, q10)
jitters_func <- function(pars){
  # rand <- runif(length(pars), 0, 1) + 0.1  # one random factor per element
  # pars * rand
  new_pars = pars + rnorm(length(pars), mean = 0, sd = 0.1 * pars) #small 10% deviation from current values
  new_pars = pmax(new_pars, 1e-6)
  new_pars = pmin(new_pars, 10)
  return(new_pars)
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
safe_find_mle <- function(model, params,
                          lower = rep(0,6), upper = rep(10,6),
                          timeout = 10000, max_attempts = 1000) {
  
  for (i in seq_len(max_attempts)) {
    cat("Attempt", i, "\n")
    
    result <- tryCatch({
      withTimeout({
        nan_flag <- FALSE
        #
        fit <- withCallingHandlers({
          diversitree::find.mle(model, params, method = "subplex",
                                lower = lower, upper = upper,
                                control = list(reltol = 1e-5, maxit = 10000))
        }, warning = function(w) {
          if (grepl("NaNs produced", w$message)) nan_flag <<- TRUE
          invokeRestart("muffleWarning")
        })
        
        # If NaN occurred, treat as failure
        if (nan_flag || any(!is.finite(fit$lnLik))) stop("NaN or invalid likelihood") #likelihood fails to produce valid result
        # If not reached convergence, treat as failure
        if (!is.null(fit$convergence) && fit$convergence != 0) stop("Optimizer did not converge") #ensure to reach convergence 
        fit
      }, timeout = timeout, onTimeout = "error")
    }, error = function(e) {
      cat("Attempt", i, "failed, timed out, has not reached convergence or NaN produced:", e$message, "\n")
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

# create data frame for storing information from each tree in this slow_rate dataset
treeset_slow_df = data.frame(tree = rep(NA,total_sim), tip_count = rep(NA,total_sim), 
                             root_age = rep(NA,total_sim), sim_type = rep(NA,total_sim),
                             lambda0_truth = rep(NA,total_sim),lambda0_est = rep(NA,total_sim),
                             lambda1_truth = rep(NA,total_sim),lambda1_est = rep(NA,total_sim),
                             mu0_truth = rep(NA,total_sim), mu0_est = rep(NA,total_sim),
                             mu1_truth = rep(NA,total_sim), mu1_est = rep(NA,total_sim),
                             q01_truth = rep(NA,total_sim), q01_est = rep(NA,total_sim),
                             q10_truth = rep(NA,total_sim), q10_est = rep(NA,total_sim),
                             freq0_obs = rep(NA,total_sim), freq0_statio = rep(NA,total_sim),
                             freq1_obs = rep(NA,total_sim), freq1_statio = rep(NA,total_sim),
                             distance_to_statio_0 = rep(NA,total_sim))


for (t in 1:length(tree_sizes)){ #loop over each tree size 
  for (i in 1:num_sim){ #loop over each tree with the given tree size 
    set.seed(1) # set seed so each search will be consistent per sampled tree
    # load the tree 
    tree_dat   = read.tree(paste0(in_fp,"tree_slow_",tree_sizes[t],"taxa_",i,".tre"))
    # add information about the tree to data_frame
    treeset_slow_df$tree[num_sim*(t-1)+i]      = paste0("tree ", i, " with size ", tree_sizes[t])
    treeset_slow_df$tip_count[num_sim*(t-1)+i] = length(tree_dat$tip.label)
    treeset_slow_df$root_age[num_sim*(t-1)+i]  = max(node.depth.edgelength(tree_dat))
    treeset_slow_df$sim_type[num_sim*(t-1)+i]  = "slow"
    # load character data at tips after pruning
    char_dat   = read.table(paste0(in_fp,"tipstates_slow_",tree_sizes[t],"taxa_",i,".csv"),sep=";",header = T)
    # load the character states before pruning
    full_char  = read.table(paste0(in_fp,"all_tipstates_slow_",tree_sizes[t],"taxa_",i,".csv"),sep=";",header = T)
    # compute how many species in 0 and 1 are included in the phylogeny
    prop_0 = nrow(char_dat[char_dat$state==0,])/length(full_char[full_char$x==0,])
    prop_1 = nrow(char_dat[char_dat$state==1,])/length(full_char[full_char$x==1,])
    # add observed frequencies from the tree
    treeset_slow_df$freq0_obs[num_sim*(t-1)+i] = sum(char_dat$state==0)/nrow(char_dat)
    treeset_slow_df$freq1_obs[num_sim*(t-1)+i] = 1-treeset_slow_df$freq0_obs[num_sim*(t-1)+i]
    # convert to a vector
    char_dat_mod        = char_dat$state
    names(char_dat_mod) = char_dat$tip
    char_dat_mod <- char_dat_mod[tree_dat$tip.label] #reorder according to tip name order
    # load the true generating params for that simulation
    true_param = read.table(paste0(in_fp,"pars_slow_",tree_sizes[t],"taxa_",i,".csv"),sep=";",header = T)
    # add true generating param to dataframe
    treeset_slow_df$lambda0_truth[num_sim*(t-1)+i] = true_param$lambda0
    treeset_slow_df$lambda1_truth[num_sim*(t-1)+i] = true_param$lambda1
    treeset_slow_df$mu0_truth[num_sim*(t-1)+i]     = true_param$mu0
    treeset_slow_df$mu1_truth[num_sim*(t-1)+i]     = true_param$mu1
    treeset_slow_df$q01_truth[num_sim*(t-1)+i]     = true_param$q01
    treeset_slow_df$q10_truth[num_sim*(t-1)+i]     = true_param$q10
    # create bisse object for that tree
    bisse_class = make.bisse(tree_dat,char_dat_mod,sampling.f = c(prop_0,prop_1))
    # fit bisse 
    # draw initial starting point
    lambda_0 = runif(1,0,0.1)
    lambda_1 = runif(1,0,0.1)
    mu_0     = runif(1,0,0.1)
    mu_1     = runif(1,0,0.1)
    q_01     = runif(1,0,0.1)
    q_10     = runif(1,0,0.1)
    # create params object
    params= c(lambda_0,lambda_1,
              mu_0,mu_1,
              q_01,q_10)
    # create the likelihood function
    best_fit = safe_find_mle(bisse_class,params)
    best_lik = best_fit$lnLik
    params   = best_fit$par # start with the optimized param after the intial draw
    # search the max likelihood for that tree 
    for (j in 2:num_starts){
      print(paste0("running iteration ",j))
      # draw new starting point for both clades for each character
      lambda_0 = runif(1,0,0.1)
      lambda_1 = runif(1,0,0.1)
      mu_0     = runif(1,0,0.1)
      mu_1     = runif(1,0,0.1)
      q_01     = runif(1,0,0.1)
      q_10     = runif(1,0,0.1)
      # create params object
      new_params= c(lambda_0,lambda_1,
                mu_0,mu_1,
                q_01,q_10)
      # Compute the new likelihood
      fit = safe_find_mle(bisse_class,new_params)
      #
      if (fit$lnLik > best_lik){ #update the likelihood and its params
        best_lik = fit$lnLik
        best_fit = fit
        print(paste0("accept new param for tree ",i, " with size ", tree_sizes[t], " from search ",j))
        params  = best_fit$par #accept the new values and use it as the new starting values
      }
      print(paste0("accept current param for tree ",i, " with size ", tree_sizes[t]," from search ",j))
    }
    # added the optima parameters to treeset_slow_df
    treeset_slow_df$lambda0_est[num_sim*(t-1)+i] = params[1]
    treeset_slow_df$lambda1_est[num_sim*(t-1)+i] = params[2]
    treeset_slow_df$mu0_est[num_sim*(t-1)+i]     = params[3]
    treeset_slow_df$mu1_est[num_sim*(t-1)+i]     = params[4]
    treeset_slow_df$q01_est[num_sim*(t-1)+i]     = params[5] 
    treeset_slow_df$q10_est[num_sim*(t-1)+i]     = params[6]
    # added theoretical stationary frequency using the optimal params 
    treeset_slow_df$freq0_statio[num_sim*(t-1)+i] = compute_statio(params)
    treeset_slow_df$freq1_statio[num_sim*(t-1)+i] = 1-treeset_slow_df$freq0_statio[num_sim*(t-1)+i]
    # compute abs distance between obs_0 and statio_0 for this tree
    treeset_slow_df$distance_to_statio_0[num_sim*(t-1)+i] = abs(treeset_slow_df$freq0_obs[num_sim*(t-1)+i]-treeset_slow_df$freq0_statio[num_sim*(t-1)+i])
  }
}

# save output
write.table(treeset_slow_df,paste0(out_fp,"bisse_slow_miss40_summary_sampling.csv"),sep = ";",row.names = F)
