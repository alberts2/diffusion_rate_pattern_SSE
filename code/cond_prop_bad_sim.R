# this script computes conditional probability of getting bad estimates across tree samples given
# outliers or not and other arguments defined in the function 

# Oct 17: changed from using relMAE and relRMSE -> MAE and RMSE, 
# because using the former, the opt threshold tends to be quite big e.g. 50% difference between true and est parameter values
# to get the biggest difference between P(bad|outlier) and p(bad|no outlier).
# note: we want opt_threshold to be as small as possible because that means we have more accurate parameter estimate in which 
# P(bad|outlier) > P(bad|no outlier) is at max. 

#### LOAD LIBRARIES ####
library(ggplot2)
library(gridExtra)
library(grid)

#### FILE SETTING ####
fp     = "/Users/albertsoewongsono/Documents/Code\ Testing/rate_pattern_diffusion_SSE/"
in_fp  = paste0(fp,"data/Inference_Simulation/")
out_fp = paste0(fp,"plot/")

#### INPUT DATA ####
### FAST VS SLOW RATES - NO MISSING TAXA - REJECTION ON DIST. STATES
slow_df = read.table(paste0(in_fp,"bisse_slow_summary.csv"),sep=";",header = T)
fast_df = read.table(paste0(in_fp,"bisse_fast_summary.csv"),sep=";",header = T)
comb_df = rbind(fast_df,slow_df)
#
### FAST VS SLOW RATES - MISSING TAXA - REJECTION ON DIST. STATES
# 20% missing taxa
# fast_df = read.table(paste0(in_fp,"bisse_fast_miss20_summary.csv"),sep=";",header = T)
# slow_df = read.table(paste0(in_fp,"bisse_slow_miss20_summary.csv"),sep=";",header = T)
# comb_df = rbind(fast_df,slow_df)
# # 40% missing taxa
# fast_df = read.table(paste0(in_fp,"bisse_fast_miss40_summary.csv"),sep=";",header = T)
# slow_df = read.table(paste0(in_fp,"bisse_slow_miss40_summary.csv"),sep=";",header = T)
# comb_df = rbind(fast_df,slow_df)
# # 80% missing taxa
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
#
### FAST VS SLOW RATES - NO MISSING TAXA - NO REJECTION ON DIST. STATES
# FAST BISSE
fast_df = read.table(paste0(in_fp,"bisse_fast_summary_no_rejection.csv"),sep=";",header = T)
# SLOW BISSE
slow_df = read.table(paste0(in_fp,"bisse_slow_summary_no_rejection.csv"),sep=";",header = T)
# COMBINED 
comb_df = rbind(fast_df,slow_df)
#
### FAST VS SLOW RATES - MISSING TAXA - NO REJECTION ON DIST. STATES
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


#### 
# dataset = dataset storing true and estimated values, tip state information, 
#           evolutionary scenario, etc 
# outlier_criterion = "criterion_1" (directionality), "criterion_2" (95% BI),
#                     "criterion_3" (both)
# treshold = a single value to classify "bad" vs "acceptable" estimates 
# estim_criterion = error function to classify "bad" vs "acceptable" estimates
#                   ("MAE" or "RMSE")
# outlier_type = "outlier" or "no_outlier" given an outlier_criterion

cond_prop_bad_sim <- function(dataset,outlier_criterion,treshold,
                             estim_criterion,outlier_type){
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
  # Input dataset
  comb_df = dataset
  #### CLASSIFYING EVO SCENARIOS ####
  comb_df$true_statio_0 = rep(NA,nrow(comb_df)) # true statio freqs 0
  comb_df$true_statio_1 = rep(NA,nrow(comb_df)) # true statio freqs 0
  # create an empty column to store outlier or not across samples
  comb_df$outlier       = rep(NA,nrow(comb_df)) 
  # create an empty column to store errors between estimate and true params across samples based on MAE or RMSE
  comb_df$rel_error     = rep(NA,nrow(comb_df))
  # classify bad or good estimates according to rel_MAE or rel_RMSE
  # if the errors <= treshold, considered acceptable 
  comb_df$class_error   = rep(NA,nrow(comb_df))
  for (i in 1:nrow(comb_df)){
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
    # precompute for criterion 1
    sign_true   = sign(comb_df$true_statio_0[i] - comb_df$true_statio_1[i])#sign for true statio freqs
    sign_obs    = sign(comb_df$freq0_obs[i] - comb_df$freq1_obs[i]) #sign for observed frequency (1 means \Pi_A > \Pi_B)
    sign_statio = sign(comb_df$freq0_statio[i] - comb_df$freq1_statio[i]) #sign for statio frequency (\Pi_hat_A - \Pi_hat_B)
    # precompute for criterion 2 and 3
    prob_success = comb_df$freq0_statio[i]
    n_trials     = comb_df$tip_count[i]
    stand_err    = sqrt(prob_success * (1 - prob_success) / n_trials)
    z            = qnorm(0.975)
    lower_bnd    = prob_success - z*stand_err
    upper_bnd    = prob_success + z*stand_err
    #
    if (outlier_criterion == "criterion_1"){
      # criterion 1: determine if the direction of state frequencies at tip match with its expected stationary using est parameters
      if (sign_obs == sign_statio){ #if the direction of difference in observed frequencies matches direction of difference in statio freq
        comb_df$outlier[i] = "no outlier"
      } else {
        comb_df$outlier[i] = "outlier"
      }
    } else if (outlier_criterion == "criterion_2"){
      # criterion 2: determine if f_A is outside 95% CI of binom experiment with p = \Pi_A from estimated parameters and n = num_tips
      # we use normal approximation to get the CI
      if ((comb_df$freq0_obs[i] < lower_bnd) || (comb_df$freq0_obs[i] > upper_bnd)){
        comb_df$outlier[i] = "outlier"
      } else {
        comb_df$outlier[i] = "no outlier"
      }
    } else if (outlier_criterion == "criterion_3"){
      if (((comb_df$freq0_obs[i] < lower_bnd) || (comb_df$freq0_obs[i] > upper_bnd)) && sign_obs != sign_statio){
        comb_df$outlier[i] = "outlier"
      } else {
        comb_df$outlier[i] = "no outlier"
      }
    }
    # # scale using the magnitude of the true parameters in each sample (dimensionless)
    # mean_scale          = mean(c(true_lambda_0,true_lambda_1,
    #                              true_mu_0,true_mu_1,
    #                              true_q01,true_q10))
    if (estim_criterion == "MAE"){
      MAE                   = (abs(comb_df$lambda0_truth[i]-comb_df$lambda0_est[i]) + abs(comb_df$lambda1_truth[i]-comb_df$lambda1_est[i]) +
                               abs(comb_df$mu0_truth[i]-comb_df$mu0_est[i]) + abs(comb_df$mu1_truth[i]-comb_df$mu1_est[i]) + 
                               abs(comb_df$q01_truth[i]-comb_df$q01_est[i]) + abs(comb_df$q10_truth[i]-comb_df$q10_est[i]))/6
      # comb_df$rel_error[i]  = MAE/mean_scale
      comb_df$rel_error[i]  = MAE
    } else if (estim_criterion == "RMSE"){
      RMSE                = sqrt(((comb_df$lambda0_truth[i]-comb_df$lambda0_est[i])^2 + (comb_df$lambda1_truth[i]-comb_df$lambda1_est[i])^2 +
                                    (comb_df$mu0_truth[i]-comb_df$mu0_est[i])^2 + (comb_df$mu1_truth[i]-comb_df$mu1_est[i])^2 + 
                                    (comb_df$q01_truth[i]-comb_df$q01_est[i])^2 + (comb_df$q10_truth[i]-comb_df$q10_est[i])^2)/6)
      # comb_df$rel_error[i] = RMSE/mean_scale
      comb_df$rel_error[i] = RMSE
    }
    # classify whether the estimates are acceptable or not
    if (comb_df$rel_error[i] <= treshold){
      comb_df$class_error[i] = "acceptable"
    } else {
      comb_df$class_error[i] = "bad"
    }
  }
  # compute percentage of sample points classified as "bad" given it is outlier or not outlier based on outlier_criterion
  # here we compute 
  if (outlier_type == "outlier"){
    # number of samples that are bad and also outlier
    num_bad_outlier = length(which(comb_df$class_error == "bad" & comb_df$outlier == "outlier"))
    # number of samples that are outlier, but not necessarily bad 
    num_outlier     = length(which(comb_df$outlier == "outlier"))
    #
    percent_bad     = (num_bad_outlier/num_outlier)*100
  } else if (outlier_type == "no_outlier"){
    # number of samples that are bad and not outlier 
    num_bad_no_outlier = length(which(comb_df$class_error == "bad" & comb_df$outlier == "no outlier"))
    # number of samples that are not outlier, but not necessarily bad 
    num_no_outlier     = length(which(comb_df$outlier == "no outlier"))
    #
    percent_bad        = (num_bad_no_outlier/num_no_outlier)*100
  }
  return(percent_bad)
}

# Do p(bad|outlier) & p(bad|no_outlier) for different treshold levels using both MAE and RMSE
# then do the same plot using different outlier criterions 

# num_tres     = 100
num_tres     = 1000
# treshold_vec = seq(0.1,1,length=num_tres)
treshold_vec = seq(0.001,1,length=num_tres)

plot_dat     = data.frame(treshold_id = rep(NA,num_tres), treshold_val = rep(NA,num_tres),
                          bad_MAE_outlier_crit_1 = rep(NA,num_tres), bad_MAE_no_outlier_crit_1 = rep(NA,num_tres),
                          bad_RMSE_outlier_crit_1 = rep(NA,num_tres), bad_RMSE_no_outlier_crit_1 = rep(NA,num_tres),
                          bad_MAE_outlier_crit_2 = rep(NA,num_tres), bad_MAE_no_outlier_crit_2 = rep(NA,num_tres),
                          bad_RMSE_outlier_crit_2 = rep(NA,num_tres), bad_RMSE_no_outlier_crit_2 = rep(NA,num_tres),
                          bad_MAE_outlier_crit_3 = rep(NA,num_tres), bad_MAE_no_outlier_crit_3 = rep(NA,num_tres),
                          bad_RMSE_outlier_crit_3 = rep(NA,num_tres), bad_RMSE_no_outlier_crit_3 = rep(NA,num_tres)
                          ) 

for (i in 1:num_tres){
  plot_dat$treshold_id[i]  = i
  plot_dat$treshold_val[i] = treshold_vec[i]
  #
  plot_dat$bad_MAE_outlier_crit_1[i]     = cond_prop_bad_sim(comb_df,"criterion_1",treshold_vec[i],"MAE","outlier")
  plot_dat$bad_MAE_no_outlier_crit_1[i]  = cond_prop_bad_sim(comb_df,"criterion_1",treshold_vec[i],"MAE","no_outlier")
  #
  plot_dat$bad_RMSE_outlier_crit_1[i]    = cond_prop_bad_sim(comb_df,"criterion_1",treshold_vec[i],"RMSE","outlier")
  plot_dat$bad_RMSE_no_outlier_crit_1[i] = cond_prop_bad_sim(comb_df,"criterion_1",treshold_vec[i],"RMSE","no_outlier")
  #
  plot_dat$bad_MAE_outlier_crit_2[i]     = cond_prop_bad_sim(comb_df,"criterion_2",treshold_vec[i],"MAE","outlier")
  plot_dat$bad_MAE_no_outlier_crit_2[i]  = cond_prop_bad_sim(comb_df,"criterion_2",treshold_vec[i],"MAE","no_outlier")
  #
  plot_dat$bad_RMSE_outlier_crit_2[i]    = cond_prop_bad_sim(comb_df,"criterion_2",treshold_vec[i],"RMSE","outlier")
  plot_dat$bad_RMSE_no_outlier_crit_2[i] = cond_prop_bad_sim(comb_df,"criterion_2",treshold_vec[i],"RMSE","no_outlier")
  #
  plot_dat$bad_MAE_outlier_crit_3[i]     = cond_prop_bad_sim(comb_df,"criterion_3",treshold_vec[i],"MAE","outlier")
  plot_dat$bad_MAE_no_outlier_crit_3[i]  = cond_prop_bad_sim(comb_df,"criterion_3",treshold_vec[i],"MAE","no_outlier")
  #
  plot_dat$bad_RMSE_outlier_crit_3[i]    = cond_prop_bad_sim(comb_df,"criterion_3",treshold_vec[i],"RMSE","outlier")
  plot_dat$bad_RMSE_no_outlier_crit_3[i] = cond_prop_bad_sim(comb_df,"criterion_3",treshold_vec[i],"RMSE","no_outlier")
  print(paste0("done for ",i))
}

#### PLOTTING ####
#### USING MAE ####
# Prob(bad|outlier,crit_1,MAE) vs Prob(bad|no outlier,crit_1,MAE)

# compute "total difference across thresholds"
# I.e. int_{threshold}|P(bad|outlier) - P(bad|not outlier|)
total_diff_crit1_MAE =  trapz(plot_dat$treshold_val,abs(plot_dat$bad_MAE_outlier_crit_1-plot_dat$bad_MAE_no_outlier_crit_1))

# compute the optima threshold where the difference between P(bad|outlier) - P(bad|no outlier) is at max
# Note I do not want abs, because I want to only choose threshold where P(bad|outlier) > P(bad|no outlier)
best_threshold_crit1_MAE = which.max(plot_dat$bad_MAE_outlier_crit_1-plot_dat$bad_MAE_no_outlier_crit_1)
if (length(best_threshold_crit1_MAE) == 0 || all(plot_dat$bad_MAE_outlier_crit_1-plot_dat$bad_MAE_no_outlier_crit_1<0)){
  best_threshold_crit1_MAE = NaN
} else {
  best_threshold_crit1_MAE = plot_dat$treshold_val[best_threshold_crit1_MAE]
}

p_bad_outlier_crit1_MAE = ggplot(plot_dat, aes(x = treshold_val)) +
                                geom_line(aes(y = bad_MAE_outlier_crit_1, color = "P(bad|outlier,crit_1)"), linewidth = 1.2) +
                                geom_line(aes(y = bad_MAE_no_outlier_crit_1, color = "P(bad|no outlier,crit_1)"), linewidth = 1.2) +
                                scale_color_manual(
                                  name = "Probability",
                                  values = c(
                                    "P(bad|outlier,crit_1)" = "steelblue",
                                    "P(bad|no outlier,crit_1)" = "tomato"
                                  )
                                ) +
                                labs(x = "Threshold value", y = "% bad") +
                                theme_minimal(base_size = 14) +
                                annotate(
                                  "text",
                                  x = Inf, y = Inf,                # top-right position
                                  label = paste0("total diff. = ", round(total_diff_crit1_MAE, 3)),
                                  hjust = 1.1, vjust = 1.5,
                                  size = 5, color = "black", fontface = "italic"
                              ) +
                              annotate(
                                "text",
                                x = Inf, y = Inf,                # same x, slightly lower y via vjust
                                label = paste0("opt threshold = ", round(best_threshold_crit1_MAE, 3)),
                                hjust = 1.1, vjust = 3.0,        # increase vjust to move it below
                                size = 5, color = "black", fontface = "italic"
                              )

# Prob(bad|outlier,crit_2,MAE) vs Prob(bad|no outlier,crit_2,MAE)

# compute "total difference across thresholds"
# I.e. int_{threshold}|P(bad|outlier) - P(bad|not outlier|)
total_diff_crit2_MAE =  trapz(plot_dat$treshold_val,abs(plot_dat$bad_MAE_outlier_crit_2-plot_dat$bad_MAE_no_outlier_crit_2))

# compute the optima threshold where the difference between P(bad|outlier) - P(bad|no outlier) is at max
# Note I do not want abs, because I want to only choose threshold where P(bad|outlier) > P(bad|no outlier)
best_threshold_crit2_MAE = which.max(plot_dat$bad_MAE_outlier_crit_2-plot_dat$bad_MAE_no_outlier_crit_2)
if (length(best_threshold_crit2_MAE) == 0 || all(plot_dat$bad_MAE_outlier_crit_2-plot_dat$bad_MAE_no_outlier_crit_2<0)){ #second condition since possibly all differences are -
  best_threshold_crit2_MAE = NaN
} else {
  best_threshold_crit2_MAE = plot_dat$treshold_val[best_threshold_crit2_MAE]
}

p_bad_outlier_crit2_MAE = ggplot(plot_dat, aes(x = treshold_val)) +
  geom_line(aes(y = bad_MAE_outlier_crit_2, color = "P(bad|outlier,crit_2)"), linewidth = 1.2) +
  geom_line(aes(y = bad_MAE_no_outlier_crit_2, color = "P(bad|no outlier,crit_2)"), linewidth = 1.2) +
  scale_color_manual(
    name = "Probability",
    values = c(
      "P(bad|outlier,crit_2)" = "steelblue",
      "P(bad|no outlier,crit_2)" = "tomato"
    )
  ) +
  labs(x = "Threshold value", y = "% bad") +
  theme_minimal(base_size = 14) +
  annotate(
    "text",
    x = Inf, y = Inf,                # top-right position
    label = paste0("total diff. = ", round(total_diff_crit2_MAE, 3)),
    hjust = 1.1, vjust = 1.5,
    size = 5, color = "black", fontface = "italic"
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,                # same x, slightly lower y via vjust
    label = paste0("opt threshold = ", round(best_threshold_crit2_MAE, 3)),
    hjust = 1.1, vjust = 3.0,        # increase vjust to move it below
    size = 5, color = "black", fontface = "italic"
  )

# Prob(bad|outlier,crit_3,MAE) vs Prob(bad|no outlier,crit_3,MAE)

# compute "total difference across thresholds"
# I.e. int_{threshold}|P(bad|outlier) - P(bad|not outlier|)
total_diff_crit3_MAE =  trapz(plot_dat$treshold_val,abs(plot_dat$bad_MAE_outlier_crit_3-plot_dat$bad_MAE_no_outlier_crit_3))

# compute the optima threshold where the difference between P(bad|outlier) - P(bad|no outlier) is at max
# Note I do not want abs, because I want to only choose threshold where P(bad|outlier) > P(bad|no outlier)
best_threshold_crit3_MAE = which.max(plot_dat$bad_MAE_outlier_crit_3-plot_dat$bad_MAE_no_outlier_crit_3)
if (length(best_threshold_crit3_MAE) == 0 || all(plot_dat$bad_MAE_outlier_crit_3-plot_dat$bad_MAE_no_outlier_crit_3<0)){
  best_threshold_crit3_MAE = NaN
} else {
  best_threshold_crit3_MAE = plot_dat$treshold_val[best_threshold_crit3_MAE]
}

p_bad_outlier_crit3_MAE = ggplot(plot_dat, aes(x = treshold_val)) +
  geom_line(aes(y = bad_MAE_outlier_crit_3, color = "P(bad|outlier,crit_3)"), linewidth = 1.2) +
  geom_line(aes(y = bad_MAE_no_outlier_crit_3, color = "P(bad|no outlier,crit_3)"), linewidth = 1.2) +
  scale_color_manual(
    name = "Probability",
    values = c(
      "P(bad|outlier,crit_3)" = "steelblue",
      "P(bad|no outlier,crit_3)" = "tomato"
    )
  ) +
  labs(x = "Threshold value", y = "%bad") +
  theme_minimal(base_size = 14) +
  annotate(
    "text",
    x = Inf, y = Inf,                # top-right position
    label = paste0("total diff. = ", round(total_diff_crit3_MAE, 3)),
    hjust = 1.1, vjust = 1.5,
    size = 5, color = "black", fontface = "italic"
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,                # same x, slightly lower y via vjust
    label = paste0("opt threshold = ", round(best_threshold_crit3_MAE, 3)),
    hjust = 1.1, vjust = 3.0,        # increase vjust to move it below
    size = 5, color = "black", fontface = "italic"
  )

title_MAE = textGrob("P(bad|outlier) vs P(bad|no outlier) using MAE",
                      gp = gpar(fontsize = 16, fontface = "bold")
)

#### USING RMSE ####
# Prob(bad|outlier,crit_1,RMSE) vs Prob(bad|no outlier,crit_1,RMSE)

# compute "total difference across thresholds"
# I.e. int_{threshold}|P(bad|outlier) - P(bad|not outlier|)
total_diff_crit1_RMSE =  trapz(plot_dat$treshold_val,abs(plot_dat$bad_RMSE_outlier_crit_1-plot_dat$bad_RMSE_no_outlier_crit_1))

# compute the optima threshold where the difference between P(bad|outlier) - P(bad|no outlier) is at max
# Note I do not want abs, because I want to only choose threshold where P(bad|outlier) > P(bad|no outlier)
best_threshold_crit1_RMSE = which.max(plot_dat$bad_RMSE_outlier_crit_1-plot_dat$bad_RMSE_no_outlier_crit_1)
if (length(best_threshold_crit1_RMSE) == 0 || all(plot_dat$bad_RMSE_outlier_crit_1-plot_dat$bad_RMSE_no_outlier_crit_1<0)){
  best_threshold_crit1_RMSE = NaN
} else {
  best_threshold_crit1_RMSE = plot_dat$treshold_val[best_threshold_crit1_RMSE]
}

p_bad_outlier_crit1_RMSE = ggplot(plot_dat, aes(x = treshold_val)) +
  geom_line(aes(y = bad_RMSE_outlier_crit_1, color = "P(bad|outlier,crit_1)"), linewidth = 1.2) +
  geom_line(aes(y = bad_RMSE_no_outlier_crit_1, color = "P(bad|no outlier,crit_1)"), linewidth = 1.2) +
  scale_color_manual(
    name = "Probability",
    values = c(
      "P(bad|outlier,crit_1)" = "steelblue",
      "P(bad|no outlier,crit_1)" = "tomato"
    )
  ) +
  labs(x = "Threshold value", y = "% bad") +
  theme_minimal(base_size = 14) + 
  annotate(
    "text",
    x = Inf, y = Inf,                # top-right position
    label = paste0("total diff. = ", round(total_diff_crit1_RMSE, 3)),
    hjust = 1.1, vjust = 1.5,
    size = 5, color = "black", fontface = "italic"
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,                # same x, slightly lower y via vjust
    label = paste0("opt threshold = ", round(best_threshold_crit1_RMSE, 3)),
    hjust = 1.1, vjust = 3.0,        # increase vjust to move it below
    size = 5, color = "black", fontface = "italic"
  )

# Prob(bad|outlier,crit_2,RMSE) vs Prob(bad|no outlier,crit_2,RMSE)

# compute "total difference across thresholds"
# I.e. int_{threshold}|P(bad|outlier) - P(bad|not outlier|)
total_diff_crit2_RMSE =  trapz(plot_dat$treshold_val,abs(plot_dat$bad_RMSE_outlier_crit_2-plot_dat$bad_RMSE_no_outlier_crit_2))

# compute the optima threshold where the difference between P(bad|outlier) - P(bad|no outlier) is at max
# Note I do not want abs, because I want to only choose threshold where P(bad|outlier) > P(bad|no outlier)
best_threshold_crit2_RMSE = which.max(plot_dat$bad_RMSE_outlier_crit_2-plot_dat$bad_RMSE_no_outlier_crit_2)
if (length(best_threshold_crit2_RMSE) == 0 || all(plot_dat$bad_RMSE_outlier_crit_2-plot_dat$bad_RMSE_no_outlier_crit_2<0)){
  best_threshold_crit2_RMSE = NaN
} else {
  best_threshold_crit2_RMSE = plot_dat$treshold_val[best_threshold_crit2_RMSE]
}

p_bad_outlier_crit2_RMSE = ggplot(plot_dat, aes(x = treshold_val)) +
  geom_line(aes(y = bad_RMSE_outlier_crit_2, color = "P(bad|outlier,crit_2)"), linewidth = 1.2) +
  geom_line(aes(y = bad_RMSE_no_outlier_crit_2, color = "P(bad|no outlier,crit_2)"), linewidth = 1.2) +
  scale_color_manual(
    name = "Probability",
    values = c(
      "P(bad|outlier,crit_2)" = "steelblue",
      "P(bad|no outlier,crit_2)" = "tomato"
    )
  ) +
  labs(x = "Threshold value", y = "% bad") +
  theme_minimal(base_size = 14) + 
  annotate(
    "text",
    x = Inf, y = Inf,                # top-right position
    label = paste0("total diff. = ", round(total_diff_crit2_RMSE, 3)),
    hjust = 1.1, vjust = 1.5,
    size = 5, color = "black", fontface = "italic"
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,                # same x, slightly lower y via vjust
    label = paste0("opt threshold = ", round(best_threshold_crit2_RMSE, 3)),
    hjust = 1.1, vjust = 3.0,        # increase vjust to move it below
    size = 5, color = "black", fontface = "italic"
  )

# Prob(bad|outlier,crit_3,RMSE) vs Prob(bad|no outlier,crit_3,RMSE)

# compute "total difference across thresholds"
# I.e. int_{threshold}|P(bad|outlier) - P(bad|not outlier|)
total_diff_crit3_RMSE     =  trapz(plot_dat$treshold_val,abs(plot_dat$bad_RMSE_outlier_crit_3-plot_dat$bad_RMSE_no_outlier_crit_3))

# compute the optima threshold where the difference between P(bad|outlier) - P(bad|no outlier) is at max
# Note I do not want abs, because I want to only choose threshold where P(bad|outlier) > P(bad|no outlier)
best_threshold_crit3_RMSE = which.max(plot_dat$bad_RMSE_outlier_crit_3-plot_dat$bad_RMSE_no_outlier_crit_3)
if (length(best_threshold_crit3_RMSE) == 0 || all(plot_dat$bad_RMSE_outlier_crit_3-plot_dat$bad_RMSE_no_outlier_crit_3<0)){
  best_threshold_crit3_RMSE = NaN
} else {
  best_threshold_crit3_RMSE = plot_dat$treshold_val[best_threshold_crit3_RMSE]
}

p_bad_outlier_crit3_RMSE = ggplot(plot_dat, aes(x = treshold_val)) +
  geom_line(aes(y = bad_RMSE_outlier_crit_3, color = "P(bad|outlier,crit_3)"), linewidth = 1.2) +
  geom_line(aes(y = bad_RMSE_no_outlier_crit_3, color = "P(bad|no outlier,crit_3)"), linewidth = 1.2) +
  scale_color_manual(
    name = "Probability",
    values = c(
      "P(bad|outlier,crit_3)" = "steelblue",
      "P(bad|no outlier,crit_3)" = "tomato"
    )
  ) +
  labs(x = "Threshold value", y = "%bad") +
  theme_minimal(base_size = 14) +
  annotate(
    "text",
    x = Inf, y = Inf,                # top-right position
    label = paste0("total diff. = ", round(total_diff_crit3_RMSE, 3)),
    hjust = 1.1, vjust = 1.5,
    size = 5, color = "black", fontface = "italic"
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,                # same x, slightly lower y via vjust
    label = paste0("opt threshold = ", round(best_threshold_crit3_RMSE, 3)),
    hjust = 1.1, vjust = 3.0,        # increase vjust to move it below
    size = 5, color = "black", fontface = "italic"
  )
  

title_RMSE = textGrob("P(bad|outlier) vs P(bad|no outlier) using RMSE",
                     gp = gpar(fontsize = 16, fontface = "bold")
)


pdf(paste0(out_fp,"combined_sim_bisse_bad_prob.pdf"), width = 16, height = 8)  # size in inches

# if warning shows up means one of the criterions does not have outliers, so nothing to show
grid.arrange(
  title_MAE,
  arrangeGrob(p_bad_outlier_crit1_MAE, p_bad_outlier_crit2_MAE, ncol = 2),  
  p_bad_outlier_crit3_MAE,                            
  ncol = 1,
  heights = c(0.1, 1, 1)  
)

grid.arrange(
  title_RMSE,
  arrangeGrob(p_bad_outlier_crit1_RMSE, p_bad_outlier_crit2_RMSE, ncol = 2),  
  p_bad_outlier_crit3_RMSE,                            
  ncol = 1,
  heights = c(0.1, 1, 1)  
)
dev.off()