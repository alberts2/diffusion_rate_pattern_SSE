# This is sim_bisse.R but without rejection sampling on proportion of tip states 

#### LOAD LIBRARIES ####
library(diversitree)
library(ape)

#### FILE SETTING ####
#fp     = "/Users/albertsoewongsono/Documents/Code\ Testing/rate_pattern_diffusion_SSE/"
fp     = "/storage/albert/rate_pattern_diffusion_SSE/"
out_fp = paste0(fp,"data/Simulation/no_rejection_states/")

#### PARALLELIZE ####
args <- commandArgs(trailingOnly = TRUE)
batch_size <- as.numeric(args[1]) # number of trees per batch
batch_id   <- as.numeric(args[2]) # number of batches

redraw_function <- function(lower,upper){
  lambda_0 = runif(1,lower,upper)
  lambda_1 = runif(1,lower,upper)
  mu_0     = runif(1,lower,upper)
  mu_1     = runif(1,lower,upper)
  q_01     = runif(1,lower,upper)
  q_10     = runif(1,lower,upper)
  pars = c(lambda0 = lambda_0, 
           lambda1 = lambda_1,
           mu0     = mu_0,
           mu1     = mu_1,
           q01     = q_01,
           q10     = q_10)
  return(pars)
}

# simulate old trees with varying sizes using bisse 
sim_bisse_slow <- function(batch_size,lower_bound = 0, upper_bound = 0.01,
                           tree_sizes = c(25,50,100,200,400,800)){
  for (i in 1:batch_size){
    # draw parameter values 
    set.seed(batch_id*1e6+i)
    lambda_0 = runif(1,lower_bound,upper_bound)
    lambda_1 = runif(1,lower_bound,upper_bound)
    mu_0     = runif(1,lower_bound,upper_bound)
    mu_1     = runif(1,lower_bound,upper_bound)
    q_01     = runif(1,lower_bound,upper_bound)
    q_10     = runif(1,lower_bound,upper_bound)
    #
    pars = c(lambda0 = lambda_0, 
             lambda1 = lambda_1,
             mu0     = mu_0,
             mu1     = mu_1,
             q01     = q_01,
             q10     = q_10)
    # simulate trees with different sizes
    # all conditioned on starting from the same root state A
    for (size in tree_sizes){
      tree_slow = tree.bisse(pars,max.taxa = size,x0=0)
      # redraw parameters if empty trees or if the distribution of tip states are very imbalanced
      while(is.null(tree_slow)){
        # dont'set seed here, so does not get stuck
        pars      = redraw_function(lower_bound,upper_bound)
        tree_slow = tree.bisse(pars,max.taxa = size,x0=0)
      }
      # SAVE GENERATING PARAMETERS
      # convert to dataframe
      pars_convert = as.data.frame(t(pars))
      #
      write.table(pars_convert,paste0(out_fp,"slow_rates/",size,"_taxa/","pars_slow_",size,"taxa_",batch_size*(batch_id-1)+i,".csv"),sep = ";",row.names = F)
      # SAVE TREES
      write.tree(tree_slow,paste0(out_fp,"slow_rates/",size,"_taxa/","tree_slow_",size,"taxa_",batch_size*(batch_id-1)+i,".tre"))
      # SAVE CHARACTER DATA AT TIPS
      tip_states <- data.frame(
        tip   = tree_slow$tip.label,
        state = tree_slow$tip.state
      )
      write.table(tip_states,paste0(out_fp, "slow_rates/", size, "_taxa/","tipstates_slow_", size, "taxa_", batch_size*(batch_id-1)+i, ".csv"),sep = ";",row.names = FALSE)
      #
      print(paste0("simulation ", i ," from batch id ",batch_id, " for slow tree with ",size," taxa finished"))
    }
  }
}

sim_bisse_fast <- function(batch_size,lower_bound = 0.1, upper_bound = 1,
                           tree_sizes = c(25,50,100,200,400,800)){
  for (i in 1:batch_size){
    # draw parameter values 
    set.seed(batch_id*1e6+i)
    lambda_0 = runif(1,lower_bound,upper_bound)
    lambda_1 = runif(1,lower_bound,upper_bound)
    mu_0     = runif(1,lower_bound,upper_bound)
    mu_1     = runif(1,lower_bound,upper_bound)
    q_01     = runif(1,lower_bound,upper_bound)
    q_10     = runif(1,lower_bound,upper_bound)
    #
    pars = c(lambda0 = lambda_0, 
             lambda1 = lambda_1,
             mu0     = mu_0,
             mu1     = mu_1,
             q01     = q_01,
             q10     = q_10)
    # simulate trees with different sizes
    # all conditioned on starting from the same root state A
    for (size in tree_sizes){
      tree_fast = tree.bisse(pars,max.taxa = size,x0=0)
      # redraw parameters if empty trees or if the distribution of tip states are very imbalanced
      while(is.null(tree_fast)){
        # dont'set seed here, so does not get stuck
        pars      = redraw_function(lower_bound,upper_bound)
        tree_fast = tree.bisse(pars,max.taxa = size,x0=0)
      }
      # SAVE GENERATING PARAMETERS
      # convert to dataframe
      pars_convert = as.data.frame(t(pars))
      #
      write.table(pars_convert,paste0(out_fp,"fast_rates/",size,"_taxa/","pars_fast_",size,"taxa_",batch_size*(batch_id-1)+i,".csv"),sep = ";",row.names = F)
      # SAVE TREES
      write.tree(tree_fast,paste0(out_fp,"fast_rates/",size,"_taxa/","tree_fast_",size,"taxa_",batch_size*(batch_id-1)+i,".tre"))
      # SAVE CHARACTER DATA AT TIPS
      tip_states <- data.frame(
        tip   = tree_fast$tip.label,
        state = tree_fast$tip.state
      )
      write.table(tip_states,paste0(out_fp, "fast_rates/", size, "_taxa/","tipstates_fast_", size, "taxa_", batch_size*(batch_id-1)+i, ".csv"),sep = ";",row.names = FALSE)
      #
      print(paste0("simulation ", i ," from batch id ",batch_id, " for fast tree with ",size," taxa finished"))
    }
  }
}

### RUNNING SIMULATIONS #####
sim_bisse_slow(batch_size)
sim_bisse_fast(batch_size)
