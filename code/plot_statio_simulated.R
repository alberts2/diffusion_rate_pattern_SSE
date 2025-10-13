# This script plots relative clade age, clade size, and distance to stationary across phylogenies 
# using inference result (for parameter estimates) from BiSSE analysis 

# NOTE: 26 Sep
# Change Euclidean distance to absolute distance in one of the states
# added color-code according to evolutionary scenarios

#### LOAD LIBRARIES ####
library(ggplot2)
library(gridExtra)
library(ape)
library(scatterplot3d)

#### FILE SETTING ####
fp     = "/Users/albertsoewongsono/Documents/Code\ Testing/rate_pattern_diffusion_SSE/"
in_fp  = paste0(fp,"data/Inference_Simulation/")
out_fp = paste0(fp,"plot/")

#### LOAD INFERENCE FILES (case where simulation setting = inference setting) ####
# FAST BISSE
fast_df = read.table(paste0(in_fp,"bisse_fast_summary.csv"),sep=";",header = T)
# SLOW BISSE
slow_df = read.table(paste0(in_fp,"bisse_slow_summary.csv"),sep=";",header = T)
# COMBINED 
comb_df = rbind(fast_df,slow_df)

#### LOAD INFERENCE FILES (case where simulation setting != inference setting) ####
# 20% missing taxa
fast_df = read.table(paste0(in_fp,"bisse_fast_miss20_summary.csv"),sep=";",header = T)
slow_df = read.table(paste0(in_fp,"bisse_slow_miss20_summary.csv"),sep=";",header = T)
comb_df = rbind(fast_df,slow_df)
# 40% missing taxa
fast_df = read.table(paste0(in_fp,"bisse_fast_miss40_summary.csv"),sep=";",header = T)
slow_df = read.table(paste0(in_fp,"bisse_slow_miss40_summary.csv"),sep=";",header = T)
comb_df = rbind(fast_df,slow_df)
# 80% missing taxa
fast_df = read.table(paste0(in_fp,"bisse_fast_miss80_summary.csv"),sep=";",header = T)
slow_df = read.table(paste0(in_fp,"bisse_slow_miss80_summary.csv"),sep=";",header = T)
comb_df = rbind(fast_df,slow_df)
# # Combined (20%, 40%, and 80%)
fast_df_80 = read.table(paste0(in_fp,"bisse_fast_miss80_summary.csv"),sep=";",header = T)
slow_df_80 = read.table(paste0(in_fp,"bisse_slow_miss80_summary.csv"),sep=";",header = T)
#
fast_df_40 = read.table(paste0(in_fp,"bisse_fast_miss40_summary.csv"),sep=";",header = T)
slow_df_40 = read.table(paste0(in_fp,"bisse_slow_miss40_summary.csv"),sep=";",header = T)
#
fast_df_20 = read.table(paste0(in_fp,"bisse_fast_miss20_summary.csv"),sep=";",header = T)
slow_df_20 = read.table(paste0(in_fp,"bisse_slow_miss20_summary.csv"),sep=";",header = T)

comb_df = rbind(fast_df_80, fast_df_40, fast_df_20,
                slow_df_80, slow_df_40, slow_df_20)


#### LOAD INFERENCE FILES (case where simulation setting = inference setting and no rejection on tip states) ####
# FAST BISSE
fast_df = read.table(paste0(in_fp,"bisse_fast_summary_no_rejection.csv"),sep=";",header = T)
# SLOW BISSE
slow_df = read.table(paste0(in_fp,"bisse_slow_summary_no_rejection.csv"),sep=";",header = T)
# COMBINED 
comb_df = rbind(fast_df,slow_df)

#### LOAD INFERENCE FILES (case where simulation setting != inference setting, no rejection on tip states) ####
# # Combined (20%, 40%, and 80%)
fast_df_80 = read.table(paste0(in_fp,"bisse_fast_miss80_summary_no_rejection.csv"),sep=";",header = T)
slow_df_80 = read.table(paste0(in_fp,"bisse_slow_miss80_summary_no_rejection.csv"),sep=";",header = T)
#
fast_df_40 = read.table(paste0(in_fp,"bisse_fast_miss40_summary_no_rejection.csv"),sep=";",header = T)
slow_df_40 = read.table(paste0(in_fp,"bisse_slow_miss40_summary_no_rejection.csv"),sep=";",header = T)
#
fast_df_20 = read.table(paste0(in_fp,"bisse_fast_miss20_summary_no_rejection.csv"),sep=";",header = T)
slow_df_20 = read.table(paste0(in_fp,"bisse_slow_miss20_summary_no_rejection.csv"),sep=";",header = T)

comb_df = rbind(fast_df_80, fast_df_40, fast_df_20,
                slow_df_80, slow_df_40, slow_df_20)


#### COMPUTE STATIO FREQS ####
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

comb_df$true_statio_0 = rep(NA,nrow(comb_df)) # true statio freqs 0
comb_df$true_statio_1 = rep(NA,nrow(comb_df)) # true statio freqs 0

comb_df$est_evo        = rep(NA,nrow(comb_df)) # estimate evolutionary scenario based on est parameters and est statio freqs
comb_df$true_evo       = rep(NA,nrow(comb_df)) # true evolutionary scenario based on true parameters and true statio freqs 
comb_df$mismatch_class = rep(NA,nrow(comb_df)) # mismatch class based on expected stationary freqs and true stationary freqs
comb_df$outlier_direct = rep(NA,nrow(comb_df)) # outlier based on sign between \Pi_0 - \Pi_1 with obs_Pi_0 - obs_Pi_1
comb_df$outlier_dist   = rep(NA,nrow(comb_df)) # outlier based on whether not f_A is outside of 95% CI of binom(n=num_tips, p=\hat{\Pi}_A)
comb_df$outlier_both   = rep(NA,nrow(comb_df)) # outlier if it meets outside 95% CI and directionality differs
# relative Mean Abs Error (MAE) and Root Mean Squared Error (RMSE) for the parameters
# i.e. normalize with the average across true parameters values 
comb_df$rel_MAE        = rep(NA,nrow(comb_df))
comb_df$rel_RMSE       = rep(NA,nrow(comb_df))
# classify bad or good estimates according to rel_MAE and rel_RMSE
# if the errors <= 0.3, considered acceptable 
comb_df$class_MAE      = rep(NA,nrow(comb_df))
comb_df$class_RMSE     = rep(NA,nrow(comb_df))

#
for (i in 1:nrow(comb_df)){ #loop over each clade
  # Compute true statio freqs
  true_lambda_0 = comb_df$lambda0_truth[i]
  true_lambda_1 = comb_df$lambda1_truth[i]
  true_mu_0     = comb_df$mu0_truth[i]
  true_mu_1     = comb_df$mu1_truth[i]
  true_q01      = comb_df$q01_truth[i]
  true_q10      = comb_df$q10_truth[i]
  comb_df$true_statio_0[i] = compute_statio(c(true_lambda_0,true_lambda_1,true_mu_0,true_mu_1,
                                              true_q01,true_q10))
  comb_df$true_statio_1[i] = 1-comb_df$true_statio_0[i]
  #
  sign_true   = sign(comb_df$true_statio_0[i] - comb_df$true_statio_1[i])#sign for true statio freqs
  sign_obs    = sign(comb_df$freq0_obs[i] - comb_df$freq1_obs[i]) #sign for observed frequency (1 means \Pi_A > \Pi_B)
  sign_statio = sign(comb_df$freq0_statio[i] - comb_df$freq1_statio[i]) #sign for statio frequency (\Pi_hat_A - \Pi_hat_B)
  # net diversification and transition direction based on the estimated parameters 
  dir_net_diver_A = comb_df$lambda0_est[i] - comb_df$mu0_est[i] # species loss or gain in A
  dir_net_diver_B = comb_df$lambda1_est[i] - comb_df$mu1_est[i] # species loss or gain in B
  dir_transition  = comb_df$q01_est[i] - comb_df$q10_est[i] # q_01 > q_10 or q_01 < q_10
  # net diversification and transition direction based on the true parameters 
  dir_net_diver_A_truth = comb_df$lambda0_truth[i] - comb_df$mu0_truth[i] # species loss or gain in A
  dir_net_diver_B_truth = comb_df$lambda1_truth[i] - comb_df$mu1_truth[i] # species loss or gain in B
  dir_transition_truth  = comb_df$q01_truth[i] - comb_df$q10_truth[i] # q_01 > q_10 or q_01 < q_10
  # determine estimate evolutionary scenario
  if (dir_net_diver_B < 0 & dir_net_diver_A > 0 & dir_transition > 0 & sign_statio == -1){
    comb_df$est_evo[i] = "scenario 1a"
  } else if (dir_net_diver_B > 0 & dir_net_diver_A < 0 & sign_statio == -1){
    comb_df$est_evo[i] = "scenario 2a"
  } else if (dir_net_diver_A == dir_net_diver_B & dir_transition > 0 & sign_statio == -1){
    comb_df$est_evo[i] = "cenario 3a"
  } else if (dir_net_diver_B > 0 & dir_net_diver_A > 0 & dir_transition > 0 & sign_statio == -1){
    comb_df$est_evo[i] = "scenario 4a"
  } else if (dir_net_diver_B < 0 & dir_net_diver_A < 0 & dir_transition > 0 & sign_statio == -1){
    comb_df$est_evo[i] = "scenario 5a"
  } else if (dir_net_diver_B > 0 & dir_net_diver_A > 0 & dir_transition < 0 & sign_statio == -1){
    comb_df$est_evo[i] = "scenario 6a"
  } else if (dir_net_diver_B < 0 & dir_net_diver_A < 0 & dir_transition < 0 & sign_statio == -1){
    comb_df$est_evo[i] = "scenario 7a"
  } else if (dir_net_diver_B > 0 & dir_net_diver_A < 0 & dir_transition < 0 & sign_statio == 1){
    comb_df$est_evo[i] = "scenario 1b"
  } else if (dir_net_diver_B < 0 & dir_net_diver_A > 0 & sign_statio == 1){
    comb_df$est_evo[i] = "scenario 2b"
  } else if (dir_net_diver_A == dir_net_diver_B & dir_transition < 0 & sign_statio == 1){
    comb_df$est_evo[i] = "scenario 3b"
  } else if (dir_net_diver_B > 0 & dir_net_diver_A > 0 & dir_transition < 0 & sign_statio == 1){
    comb_df$est_evo[i] = "scenario 4b"
  } else if (dir_net_diver_B < 0 & dir_net_diver_A < 0 & dir_transition < 0 & sign_statio == 1){
    comb_df$est_evo[i] = "scenario 5b"
  } else if (dir_net_diver_B > 0 & dir_net_diver_A > 0 & dir_transition > 0 & sign_statio == 1){
    comb_df$est_evo[i] = "scenario 6b"
  } else if (dir_net_diver_B < 0 & dir_net_diver_A < 0 & dir_transition > 0 & sign_statio == 1){
    comb_df$est_evo[i] = "scenario 7b"
  } else if (dir_net_diver_B > 0 & dir_net_diver_A < 0 & dir_transition < 0 & sign_statio == 0){
    comb_df$est_evo[i] = "scenario 1c"
  } else if (dir_net_diver_B < 0 & dir_net_diver_A > 0 & dir_transition > 0 & sign_statio == 0){
    comb_df$est_evo[i] = "scenario 2c"
  } else if (dir_net_diver_B == dir_net_diver_A & dir_transition == 0 &sign_statio == 0){
    comb_df$est_evo[i] = "scenario 3c"
  }
  # determine true evolutionary scenario (using true statio freqs)
  if (dir_net_diver_B_truth < 0 & dir_net_diver_A_truth > 0 & dir_transition_truth > 0 & sign_true == -1){
    comb_df$true_evo[i] = "scenario 1a"
  } else if (dir_net_diver_B_truth > 0 & dir_net_diver_A_truth < 0 & sign_true == -1){
    comb_df$true_evo[i] = "scenario 2a"
  } else if (dir_net_diver_A_truth == dir_net_diver_B_truth & dir_transition_truth > 0 & sign_true == -1){
    comb_df$true_evo[i] = "cenario 3a"
  } else if (dir_net_diver_B_truth > 0 & dir_net_diver_A_truth > 0 & dir_transition_truth > 0 & sign_true == -1){
    comb_df$true_evo[i] = "scenario 4a"
  } else if (dir_net_diver_B_truth < 0 & dir_net_diver_A_truth < 0 & dir_transition_truth > 0 & sign_true == -1){
    comb_df$true_evo[i] = "scenario 5a"
  } else if (dir_net_diver_B_truth > 0 & dir_net_diver_A_truth > 0 & dir_transition_truth < 0 & sign_true == -1){
    comb_df$true_evo[i] = "scenario 6a"
  } else if (dir_net_diver_B_truth < 0 & dir_net_diver_A_truth < 0 & dir_transition_truth < 0 & sign_true == -1){
    comb_df$true_evo[i] = "scenario 7a"
  } else if (dir_net_diver_B_truth > 0 & dir_net_diver_A_truth < 0 & dir_transition_truth < 0 & sign_true == 1){
    comb_df$true_evo[i] = "scenario 1b"
  } else if (dir_net_diver_B_truth < 0 & dir_net_diver_A_truth > 0 & sign_true == 1){
    comb_df$true_evo[i] = "scenario 2b"
  } else if (dir_net_diver_A_truth == dir_net_diver_B_truth & dir_transition_truth < 0 & sign_true == 1){
    comb_df$true_evo[i] = "scenario 3b"
  } else if (dir_net_diver_B_truth > 0 & dir_net_diver_A_truth > 0 & dir_transition_truth < 0 & sign_true == 1){
    comb_df$true_evo[i] = "scenario 4b"
  } else if (dir_net_diver_B_truth < 0 & dir_net_diver_A_truth < 0 & dir_transition_truth < 0 & sign_true == 1){
    comb_df$true_evo[i] = "scenario 5b"
  } else if (dir_net_diver_B_truth > 0 & dir_net_diver_A_truth > 0 & dir_transition_truth > 0 & sign_true == 1){
    comb_df$true_evo[i] = "scenario 6b"
  } else if (dir_net_diver_B_truth < 0 & dir_net_diver_A_truth < 0 & dir_transition_truth > 0 & sign_true == 1){
    comb_df$true_evo[i] = "scenario 7b"
  } else if (dir_net_diver_B_truth > 0 & dir_net_diver_A_truth < 0 & dir_transition_truth < 0 & sign_true == 0){
    comb_df$true_evo[i] = "scenario 1c"
  } else if (dir_net_diver_B_truth < 0 & dir_net_diver_A_truth > 0 & dir_transition_truth > 0 & sign_true == 0){
    comb_df$true_evo[i] = "scenario 2c"
  } else if (dir_net_diver_B_truth == dir_net_diver_A_truth & dir_transition_truth == 0 &sign_true == 0){
    comb_df$true_evo[i] = "scenario 3c"
  } else {
    comb_df$true_evo[i] = "none"
  }
  # determine whether there is a "mismatch" between expected evo scenario under estimated parameters and true evo scenario
  # under true generating parameters. There are three types of mismatches:
  # 1) mismatch only in directionality of statio freqs (determined by different letters with the same number) (miss_class)
  # 2) mismatch only in evo scenario within the same evo scenario (determined by different numbers with the same letter) (miss_evo)
  # 3) mismatch in both (determined by different numbers and letters) (miss_both)
  # there are 1 type of match:
  # 1) match in both directionality class and evo scenarios (same number and letter) (match_both)
  #
  # Extract the directionality class in both est_evo and true_evo (the letter)
  direct_est  = sub(".*[0-9]([a-zA-Z])", "\\1", comb_df$est_evo[i])
  direct_true = sub(".*[0-9]([a-zA-Z])", "\\1", comb_df$true_evo[i])
  # Extract the evolutionary scenario in both est_evo and true_evo (the number)
  evo_est   = as.numeric(sub(".* ([0-9]+)[a-zA-Z]", "\\1", comb_df$est_evo[i]))
  evo_true  = as.numeric(sub(".* ([0-9]+)[a-zA-Z]", "\\1", comb_df$true_evo[i]))
  # Determine which mismatch class 
  if ((direct_est != direct_true) && (evo_est == evo_true)){
    comb_df$mismatch_class[i] = "miss_class"
  } else if ((direct_est == direct_true) && (evo_est != evo_true)){
    comb_df$mismatch_class[i] = "miss_evo"
  } else if ((direct_est != direct_true) && (evo_est != evo_true)){
    comb_df$mismatch_class[i] = "miss_both"
  } else if ((direct_est == direct_true) && (evo_est == evo_true)){
    comb_df$mismatch_class[i] = "match_both"
  }
  # criterion 1: determine if the direction of state frequencies at tip match with its expected stationary using est parameters
  if (sign_obs == sign_statio){ #if the direction of difference in observed frequencies matches direction of difference in statio freq
    comb_df$outlier_direct[i] = "no outlier"
  } else {
    comb_df$outlier_direct[i] = "outlier"
  }
  # criterion 2: determine if f_A is outside 95% CI of binom experiment with p = \Pi_A from estimated parameters and n = num_tips
  # we use normal approximation to get the CI
  prob_success = comb_df$freq0_statio[i]
  n_trials     = comb_df$tip_count[i]
  stand_err    = sqrt(prob_success * (1 - prob_success) / n_trials)
  z            = qnorm(0.975)
  lower_bnd    = prob_success - z*stand_err
  upper_bnd    = prob_success + z*stand_err
  #
  if ((comb_df$freq0_obs[i] < lower_bnd) || (comb_df$freq0_obs[i] > upper_bnd)){
    comb_df$outlier_dist[i] = "outlier"
  } else {
    comb_df$outlier_dist[i] = "no outlier"
  }
  # criterion 3: outlier if f_A is outside 95% CI and f_A-f_B has different sign that est_pi_A - est_pi_B
  if (((comb_df$freq0_obs[i] < lower_bnd) || (comb_df$freq0_obs[i] > upper_bnd)) && sign_obs != sign_statio){
    comb_df$outlier_both[i] = "outlier"
  } else {
    comb_df$outlier_both[i] = "no outlier"
  }
  # compute rel. MAE and rel. RMSE
  # scaled to the magnitude of the true parameters in each sample (dimensionless)
  mean_scale          = mean(c(true_lambda_0,true_lambda_1,
                               true_mu_0,true_mu_1,
                               true_q01,true_q10))
  #
  MAE                 = (abs(comb_df$lambda0_truth[i]-comb_df$lambda0_est[i]) + abs(comb_df$lambda1_truth[i]-comb_df$lambda1_est[i]) +
                           abs(comb_df$mu0_truth[i]-comb_df$mu0_est[i]) + abs(comb_df$mu1_truth[i]-comb_df$mu1_est[i]) + 
                           abs(comb_df$q01_truth[i]-comb_df$q01_est[i]) + abs(comb_df$q10_truth[i]-comb_df$q10_est[i]))/6
  comb_df$rel_MAE[i]  = MAE/mean_scale
  #
  RMSE                = sqrt(((comb_df$lambda0_truth[i]-comb_df$lambda0_est[i])^2 + (comb_df$lambda1_truth[i]-comb_df$lambda1_est[i])^2 +
                                (comb_df$mu0_truth[i]-comb_df$mu0_est[i])^2 + (comb_df$mu1_truth[i]-comb_df$mu1_est[i])^2 + 
                                (comb_df$q01_truth[i]-comb_df$q01_est[i])^2 + (comb_df$q10_truth[i]-comb_df$q10_est[i])^2)/6)
  comb_df$rel_RMSE[i] = RMSE/mean_scale
  # classify whether the estimates are acceptable or not based on relMAE and and relRMSE 
  if (comb_df$rel_MAE[i] <= 0.3){
    comb_df$class_MAE[i] = "acceptable"
  } else {
    comb_df$class_MAE[i] = "bad"
  }
  #
  if (comb_df$rel_RMSE[i] <= 0.3){
    comb_df$class_RMSE[i] = "acceptable"
  } else {
    comb_df$class_RMSE[i] = "bad"
  }
}

#### GET RELEVANT INFORMATION ####
summary_df = comb_df[,c("tree","tip_count","root_age",
                        "sim_type","distance_to_statio_0","est_evo",
                        "true_evo","mismatch_class","freq0_obs","freq0_statio","true_statio_0",
                        "freq1_obs","freq1_statio","true_statio_1",
                        "outlier_direct","outlier_dist","outlier_both","rel_MAE","rel_RMSE",
                        "class_MAE","class_RMSE")]

outlier_1 = summary_df[summary_df$outlier_direct=="outlier",] #all the outliers based on criterion 1
outlier_2 = summary_df[summary_df$outlier_dist=="outlier",]   #all the outliers based on criterion 2
outlier_3 = summary_df[summary_df$outlier_both=="outlier",]   #all the outliers based on criterion 3

#### CREATE COUNT HEATMAP FOR TRUE EVO VS EXPECTED SCENARIO ACROSS SAMPLES 
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
all_scenarios = c("scenario 1a", "scenario 2a", "scenario 3a", "scenario 4a",
                  "scenario 5a", "scenario 6a", "scenario 7a", "scenario 1b",
                  "scenario 2b", "scenario 3b", "scenario 4b", "scenario 5b",
                  "scenario 6b", "scenario 7b", "scenario 1c", "scenario 2c",
                  "scenario 3c") # list of all scenarios 

heatmap_df        = expand.grid(true_scenario_list = all_scenarios,
                                est_scenario_list = all_scenarios)
heatmap_df$count  = rep(0,nrow(heatmap_df))

# fill in the count
for (i in 1:nrow(comb_df)){
  true_scenario = comb_df$true_evo[i]
  est_scenario  = comb_df$est_evo[i]
  #
  which_ind = which(heatmap_df$true_scenario_list== true_scenario & 
                    heatmap_df$est_scenario_list == est_scenario)
  #
  heatmap_df$count[which_ind] = heatmap_df$count[which_ind]+1
}

# get the proportion across each true scenario counts 
row_sums        = tapply(heatmap_df$count, heatmap_df$true_scenario_list, sum) # compute row sums for true scenarios
# lookup the corresponding rowsums for each true scenario in heatmap_df
look_up_rowsums = row_sums[heatmap_df$true_scenario_list]
# get the proportion for each combination of true vs estimated divided by total trees belong to the true scenario
heatmap_df$proportion = heatmap_df$count/look_up_rowsums
# assign 0 for division by zero problem (those true evo scenarios that do not have samples)
heatmap_df$proportion[is.na(heatmap_df$proportion)] = 0
# join count and proportion for plotting
heatmap_df$joint_count_prop <- paste0(
  heatmap_df$count,
  "\n(",
  sprintf("%.2f", heatmap_df$proportion),
  ")"
)

#### PLOT ####
# Note: currently there's a problem with the dataset, trees across different batches are the same
# e.g., tree 1 = tree 11 = tree 21 etc. 
# current redoing the fitting step. (fixed) 
#

summary_df$outlier_direct = factor(summary_df$outlier_direct)
summary_df$outlier_dist   = factor(summary_df$outlier_dist)
summary_df$outlier_both   = factor(summary_df$outlier_both)

#### PLOTTING ####

# Plot 1: outlier + errors in parameter estimation
# basically this says: outlier -> bad parameter estimation (good estimation -> no outlier)
# this direction should stay true whenever we identify an outlier according to one or more criterions defined above 
# but bad parameter estimation -\> outlier (this is not a good direction since we never know true generating parameters
# in empirical setting)

# based on rel MAE
count_outlier_direct_mae = ggplot(summary_df, aes(x = outlier_direct, fill = class_MAE)) +
                                   geom_bar(position = "fill") +  # counts stacked
                                   labs(x = "outlier based on stationary direction", y = "Proportion", fill = "MAE param estimates") +
                                   scale_fill_manual(values = c("acceptable" = "turquoise", "bad" = "#DD7070")) +
                                   theme_minimal()

count_outlier_dist_mae = ggplot(summary_df, aes(x = outlier_dist, fill = class_MAE)) +
                                geom_bar(position = "fill") +  # counts stacked
                                labs(x = "outlier based on binomial CI", y = "Proportion", fill = "MAE param estimates") +
                                scale_fill_manual(values = c("acceptable" = "turquoise", "bad" = "#DD7070")) +
                                theme_minimal()

count_outlier_both_mae = ggplot(summary_df, aes(x = outlier_both, fill = class_MAE)) +
                                geom_bar(position = "fill") +  # counts stacked
                                labs(x = "outlier based on both criterions", y = "Proportion", fill = "MAE param estimates") +
                                scale_fill_manual(values = c("acceptable" = "turquoise", "bad" = "#DD7070")) +
                                theme_minimal()
# based on rel RMSE 
count_outlier_direct_rmse = ggplot(summary_df, aes(x = outlier_direct, fill = class_RMSE)) +
                                  geom_bar(position = "fill") +  # counts stacked
                                  labs(x = "outlier based on stationary direction", y = "Proportion", fill = "RMSE param estimates") +
                                  scale_fill_manual(values = c("acceptable" = "turquoise", "bad" = "#DD7070")) +
                                  theme_minimal()

count_outlier_dist_rmse = ggplot(summary_df, aes(x = outlier_dist, fill = class_RMSE)) +
                                  geom_bar(position = "fill") +  # counts stacked
                                  labs(x = "outlier based on binomial CI", y = "Proportion", fill = "RMSE param estimates") +
                                  scale_fill_manual(values = c("acceptable" = "turquoise", "bad" = "#DD7070")) +
                                  theme_minimal()

count_outlier_both_rmse = ggplot(summary_df, aes(x = outlier_both, fill = class_RMSE)) +
                                  geom_bar(position = "fill") +  # counts stacked
                                  labs(x = "outlier based on both criterions", y = "Proportion", fill = "RMSE param estimates") +
                                  scale_fill_manual(values = c("acceptable" = "turquoise", "bad" = "#DD7070")) +
                                  theme_minimal()

# Plot 2: outlier + type of mismatch in parameter estimation

# based on criterion 1: direction
count_direct_miss = ggplot(summary_df, aes(x = outlier_direct, fill = mismatch_class)) +
                           geom_bar(position = "fill") +  # counts stacked
                           labs(x = "outlier based on stationary direction", y = "Proportion", fill = "Types of mismatch") +
                           scale_fill_manual(values = c("match_both" = "turquoise", "miss_both" = "#DD7070",
                                                        "miss_class" = "#DECBE4", "miss_evo" = "#E5D8BD")) +
                           theme_minimal()

# based on criterion 2: distribution
count_dist_miss = ggplot(summary_df, aes(x = outlier_dist, fill = mismatch_class)) +
                          geom_bar(position = "fill") +  # counts stacked
                          labs(x = "outlier based on binomial CI", y = "Proportion", fill = "Types of mismatch") +
                          scale_fill_manual(values = c("match_both" = "turquoise", "miss_both" = "#DD7070",
                                                       "miss_class" = "#DECBE4", "miss_evo" = "#E5D8BD")) +
                          theme_minimal()

# based on criterion 3: both
count_both_miss = ggplot(summary_df, aes(x = outlier_both, fill = mismatch_class)) +
                          geom_bar(position = "fill") +  # counts stacked
                          labs(x = "outlier based on both criterions", y = "Proportion", fill = "Types of mismatch") +
                          scale_fill_manual(values = c("match_both" = "turquoise", "miss_both" = "#DD7070",
                                                       "miss_class" = "#DECBE4", "miss_evo" = "#E5D8BD")) +
                          theme_minimal()

# heatmap for true scenario vs est scenario (by count)
p_heatmap_count = ggplot(heatmap_df, aes(x = est_scenario_list, y = true_scenario_list, fill = count)) +
                    geom_tile(color = "grey90") +
                    geom_text(aes(label = joint_count_prop), color = "black", size = 3) +
                    scale_fill_gradient(low = "white", high = "steelblue") +
                    labs(title = "True vs Estimated Evolutionary Scenario (by count)",
                         x = "Estimated scenario",
                         y = "True scenario") +
                    theme_minimal(base_size = 12) +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1),
                          plot.title = element_text(hjust = 0.5))
# heatmap for true scenario vs est scenario (by proportion)
p_heatmap_prop = ggplot(heatmap_df, aes(x = est_scenario_list, y = true_scenario_list, fill = proportion)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = joint_count_prop), color = "black", size = 3) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "True vs Estimated Evolutionary Scenario (by proportion)",
       x = "Estimated scenario",
       y = "True scenario") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))


### NOTE: My trees are simulated condition on tree sizes, not time 
# color = "outlier" or "not outlier" (according to three different criterions)
# shape = est_evo 
# shape = true_evo

# define global mapping for shapes to be consistent across est_evo and true_evo
all_levels <- sort(unique(c("scenario 1a", "scenario 2a", "scenario 3a",
                            "scenario 4a", "scenario 5a", "scenario 6a",
                            "scenario 7a", "scenario 1b", "scenario 2b",
                            "scenario 3b", "scenario 4b", "scenario 5b",
                            "scenario 6b", "scenario 7b", "scenario 1c",
                            "scenario 2c", "scenario 3c")))

shape_mapping <- setNames(seq(1,17), all_levels)

# Plot 2: relMAE vs Tree Size 
# relMAE vs Tree Size by est_evo (shape) and outlier_direct (color)
p_MAE_size_est_direct =  ggplot(summary_df, aes(x = tip_count, y = rel_MAE, color = outlier_direct, shape = est_evo)) +
                          geom_point(size = 3) +
                          scale_x_continuous(breaks = c(0, 25, 50, 100, 200, 400, 800)) +
                          scale_color_manual(values = c("no outlier" = "blue", "outlier" = "red")) +
                          coord_cartesian(ylim = c(0, 5)) +
                          scale_shape_manual(values = shape_mapping)+
                          theme_minimal()
# relMAE vs Tree Size by est_evo (shape) and outlier_dist (color)
p_MAE_size_est_dist =  ggplot(summary_df, aes(x = tip_count, y = rel_MAE, color = outlier_dist, shape = est_evo)) +
                              geom_point(size = 3) +
                              scale_x_continuous(breaks = c(0, 25, 50, 100, 200, 400, 800)) +
                              scale_color_manual(values = c("no outlier" = "blue", "outlier" = "red")) +
                              coord_cartesian(ylim = c(0, 5)) +
                              scale_shape_manual(values = shape_mapping)+
                              theme_minimal()
# relMAE vs Tree Size by est_evo (shape) and outlier_both (color)
p_MAE_size_est_both =  ggplot(summary_df, aes(x = tip_count, y = rel_MAE, color = outlier_both, shape = est_evo)) +
                              geom_point(size = 3) +
                              scale_x_continuous(breaks = c(0, 25, 50, 100, 200, 400, 800)) +
                              scale_color_manual(values = c("no outlier" = "blue", "outlier" = "red")) +
                              coord_cartesian(ylim = c(0, 5)) +
                              scale_shape_manual(values = shape_mapping)+
                              theme_minimal()


# relMAE vs Tree Size by true_evo (shape) and outlier_direct (color)
p_MAE_size_true_direct =  ggplot(summary_df, aes(x = tip_count, y = rel_MAE, color = outlier_direct, shape = true_evo)) +
                                  geom_point(size = 3) +
                                  scale_x_continuous(breaks = c(0, 25, 50, 100, 200, 400, 800)) +
                                  scale_color_manual(values = c("no outlier" = "blue", "outlier" = "red")) +
                                  coord_cartesian(ylim = c(0, 5)) +
                                  scale_shape_manual(values = shape_mapping) +
                                  theme_minimal()
# relMAE vs Tree Size by true_evo (shape) and outlier_dist (color)
p_MAE_size_true_dist =  ggplot(summary_df, aes(x = tip_count, y = rel_MAE, color = outlier_dist, shape = true_evo)) +
                                geom_point(size = 3) +
                                scale_x_continuous(breaks = c(0, 25, 50, 100, 200, 400, 800)) +
                                scale_color_manual(values = c("no outlier" = "blue", "outlier" = "red")) +
                                coord_cartesian(ylim = c(0, 5)) +
                                scale_shape_manual(values = shape_mapping) +
                                theme_minimal()
# relMAE vs Tree Size by true_evo (shape) and outlier_both (color)
p_MAE_size_true_both =  ggplot(summary_df, aes(x = tip_count, y = rel_MAE, color = outlier_both, shape = true_evo)) +
                                geom_point(size = 3) +
                                scale_x_continuous(breaks = c(0, 25, 50, 100, 200, 400, 800)) +
                                scale_color_manual(values = c("no outlier" = "blue", "outlier" = "red")) +
                                coord_cartesian(ylim = c(0, 5)) +
                                scale_shape_manual(values = shape_mapping) +
                                theme_minimal()

# Plot 3: relRMSE vs Tree Size 
# relRMSE vs Tree Size by est_evo (shape) and outlier_direct (color)
p_RMSE_size_est_direct = ggplot(summary_df, aes(x = tip_count, y = rel_RMSE, color = outlier_direct, shape = est_evo)) +
                                geom_point(size = 3) +
                                scale_x_continuous(breaks = c(0, 25, 50, 100, 200, 400, 800)) +
                                scale_color_manual(values = c("no outlier" = "blue", "outlier" = "red")) +
                                coord_cartesian(ylim = c(0, 5)) +
                                scale_shape_manual(values = shape_mapping) +
                                theme_minimal()
# relRMSE vs Tree Size by est_evo (shape) and outlier_dist (color)
p_RMSE_size_est_dist = ggplot(summary_df, aes(x = tip_count, y = rel_RMSE, color = outlier_dist, shape = est_evo)) +
                              geom_point(size = 3) +
                              scale_x_continuous(breaks = c(0, 25, 50, 100, 200, 400, 800)) +
                              scale_color_manual(values = c("no outlier" = "blue", "outlier" = "red")) +
                              coord_cartesian(ylim = c(0, 5)) +
                              scale_shape_manual(values = shape_mapping) +
                              theme_minimal()
# relRMSE vs Tree Size by est_evo (shape) and outlier_both (color)
p_RMSE_size_est_both = ggplot(summary_df, aes(x = tip_count, y = rel_RMSE, color = outlier_both, shape = est_evo)) +
                              geom_point(size = 3) +
                              scale_x_continuous(breaks = c(0, 25, 50, 100, 200, 400, 800)) +
                              scale_color_manual(values = c("no outlier" = "blue", "outlier" = "red")) +
                              coord_cartesian(ylim = c(0, 5)) +
                              scale_shape_manual(values = shape_mapping) +
                              theme_minimal()

# relRMSE vs Tree Size by true_evo (shape) and outlier_direct (color)
p_RMSE_size_true_direct = ggplot(summary_df, aes(x = tip_count, y = rel_RMSE, color = outlier_direct, shape = true_evo)) +
                              geom_point(size = 3) +
                              scale_x_continuous(breaks = c(0, 25, 50, 100, 200, 400, 800)) +
                              scale_color_manual(values = c("no outlier" = "blue", "outlier" = "red")) +
                              coord_cartesian(ylim = c(0, 5)) +
                              scale_shape_manual(values = shape_mapping) +
                              theme_minimal()

# relRMSE vs Tree Size by true_evo (shape) and outlier_dist (color)
p_RMSE_size_true_dist = ggplot(summary_df, aes(x = tip_count, y = rel_RMSE, color = outlier_dist, shape = true_evo)) +
                              geom_point(size = 3) +
                              scale_x_continuous(breaks = c(0, 25, 50, 100, 200, 400, 800)) +
                              scale_color_manual(values = c("no outlier" = "blue", "outlier" = "red")) +
                              coord_cartesian(ylim = c(0, 5)) +
                              scale_shape_manual(values = shape_mapping) +
                              theme_minimal()

# relRMSE vs Tree Size by true_evo (shape) and outlier_both (color)
p_RMSE_size_true_both = ggplot(summary_df, aes(x = tip_count, y = rel_RMSE, color = outlier_both, shape = true_evo)) +
                                geom_point(size = 3) +
                                scale_x_continuous(breaks = c(0, 25, 50, 100, 200, 400, 800)) +
                                scale_color_manual(values = c("no outlier" = "blue", "outlier" = "red")) +
                                coord_cartesian(ylim = c(0, 5)) +
                                scale_shape_manual(values = shape_mapping) +
                                theme_minimal()

# PLOT 4: \est_pi_A - obs_pi_A vs \est_pi_A - \true_pi_A  across clade sizes
rng       = range(c(summary_df$freq0_statio - summary_df$freq0_obs,
                    summary_df$freq0_statio - summary_df$true_statio_0)) # so axes limits are the same across all plots
#
p_diff_25 = ggplot(summary_df[summary_df$tip_count == 25,], aes(x = freq0_statio - freq0_obs, y = freq0_statio-true_statio_0)) +
                    geom_point(size = 3) + 
                    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
                    coord_equal(xlim = rng, ylim = rng)+
                    ggtitle("Tip count = 25") +
                    theme_minimal()
#
p_diff_50 = ggplot(summary_df[summary_df$tip_count == 50,], aes(x = freq0_statio - freq0_obs, y = freq0_statio-true_statio_0)) +
                    geom_point(size = 3) + 
                    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
                    coord_equal(xlim = rng, ylim = rng)+
                    ggtitle("Tip count = 50") +
                    theme_minimal()
#
p_diff_100 = ggplot(summary_df[summary_df$tip_count == 100,], aes(x = freq0_statio - freq0_obs, y = freq0_statio-true_statio_0)) +
                    geom_point(size = 3) + 
                    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
                    coord_equal(xlim = rng, ylim = rng)+
                    ggtitle("Tip count = 100") +
                    theme_minimal()
#
p_diff_200 = ggplot(summary_df[summary_df$tip_count == 200,], aes(x = freq0_statio - freq0_obs, y = freq0_statio-true_statio_0)) +
                    geom_point(size = 3) + 
                    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
                    coord_equal(xlim = rng, ylim = rng)+
                    ggtitle("Tip count = 200") +
                    theme_minimal()
#
p_diff_400 = ggplot(summary_df[summary_df$tip_count == 400,], aes(x = freq0_statio - freq0_obs, y = freq0_statio-true_statio_0)) +
                    geom_point(size = 3) + 
                    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
                    coord_equal(xlim = rng, ylim = rng)+
                    ggtitle("Tip count = 400") +
                    theme_minimal()
#
p_diff_800 = ggplot(summary_df[summary_df$tip_count == 800,], aes(x = freq0_statio - freq0_obs, y = freq0_statio-true_statio_0)) +
                    geom_point(size = 3) + 
                    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
                    coord_equal(xlim = rng, ylim = rng)+
                    ggtitle("Tip count = 800") +
                    theme_minimal()
#
p_diff_comb = grid.arrange(p_diff_25, p_diff_50, p_diff_100,
                           p_diff_200, p_diff_400, p_diff_800,
                           nrow = 2, ncol = 3)

pdf(paste0(out_fp,"combined_sim_bisse_plots.pdf"), width = 10, height = 8)  # size in inches

print(count_outlier_direct_mae)
print(count_outlier_dist_mae)
print(count_outlier_both_mae)
print(count_outlier_direct_rmse)
print(count_outlier_dist_rmse)
print(count_outlier_both_rmse)
print(count_direct_miss)
print(count_dist_miss)
print(count_both_miss)
print(p_heatmap_count)
print(p_heatmap_prop)
print(p_MAE_size_est_direct)
print(p_MAE_size_est_dist)
print(p_MAE_size_est_both)
print(p_MAE_size_true_direct)
print(p_MAE_size_true_dist)
print(p_MAE_size_true_both)
print(p_RMSE_size_est_direct)
print(p_RMSE_size_est_dist)
print(p_RMSE_size_est_both)
print(p_RMSE_size_true_direct)
print(p_RMSE_size_true_dist)
print(p_RMSE_size_true_both)
grid.arrange(p_diff_25, p_diff_50, p_diff_100,
             p_diff_200, p_diff_400, p_diff_800,
             nrow = 2, ncol = 3)

dev.off()
#
# PLOT 5
outlier_direct <- plot_ly(summary_df,
                      x = ~root_age, 
                      y = ~tip_count, 
                      z = ~distance_to_statio_0,
                      type   = "scatter3d", 
                      mode   = "markers",
                      color  = ~outlier_direct,
                      colors = c("blue", "red"),
                      symbol = ~sim_type,
                      symbols = c("circle", "square")
              ) %>%    # color by group
                layout(
                  title = "Outlier based on direction in frequencies criterion ",
                  scene = list(
                    xaxis = list(title = "root age", showbackground = FALSE),
                    yaxis = list(title = "tree size", showbackground = FALSE),
                    zaxis = list(title = "distance from true 0", showbackground = FALSE)
                  ),
                  legend = list(traceorder = "normal")  # preserve the factor level order
                )

outlier_dist <- plot_ly(summary_df,
                  x = ~root_age, 
                  y = ~tip_count, 
                  z = ~distance_to_statio_0,
                  type   = "scatter3d", 
                  mode   = "markers",
                  color  = ~outlier_dist,
                  colors = c("blue", "red"),
                  symbol = ~sim_type,
                  symbols = c("circle", "square")
          ) %>%    # color by group
            layout(
              title = "Outlier based on binomial CI criterion",
              scene = list(
                xaxis = list(title = "root age", showbackground = FALSE),
                yaxis = list(title = "tree size", showbackground = FALSE),
                zaxis = list(title = "distance from true 0", showbackground = FALSE)
              ),
              legend = list(traceorder = "normal")  # preserve the factor level order
            )
#
outlier_both <- plot_ly(summary_df,
                    x = ~root_age, 
                    y = ~tip_count, 
                    z = ~distance_to_statio_0,
                    type   = "scatter3d", 
                    mode   = "markers",
                    color  = ~outlier_both,
                    colors = c("blue", "red"),
                    symbol = ~sim_type,
                    symbols = c("circle", "square")
            ) %>%    # color by group
              layout(
                title = "Outlier based on both criterions",
                scene = list(
                  xaxis = list(title = "root age", showbackground = FALSE),
                  yaxis = list(title = "tree size", showbackground = FALSE),
                  zaxis = list(title = "distance from true 0", showbackground = FALSE)
                ),
                legend = list(traceorder = "normal")  # preserve the factor level order
            )

htmlwidgets::saveWidget(outlier_direct, paste0(out_fp,"3D_outlier_direct.html"))
htmlwidgets::saveWidget(outlier_dist, paste0(out_fp,"3D_outlier_dist.html"))
htmlwidgets::saveWidget(outlier_both, paste0(out_fp,"3D_outlier_both.html"))
