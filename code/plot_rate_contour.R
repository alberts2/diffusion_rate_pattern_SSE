# Assuming users already save outputs from get_rules.R and get_rates.R, plot the solutions (random draws of parameter)

# See the same code in ~/Code Testing/SSA/diffusion_Rpackage/ for complete history of changes made 


plot_rate_contour <- function(variant,stat_freq,choose_rule,lower_bound,upper_bound,num_repli,batch_num,rate_dist){
  ###############################
  # LOAD R PACKAGES
  ###############################
  #
  library(ggplot2)
  library(gridExtra)
  
  ####### FILE SYSTEM ######
  ##########################
  in_fp   = "/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/"
  out_fp  = paste0(in_fp,"data/")
  plot_fp = paste0(in_fp,"plot/")
  
  setwd(in_fp)
  
  ###############################
  # FIND RULESETS AND GENERATE RATES GIVEN RULESETS
  ###############################
  #
  my_rules <- get_rules(variant,stat_freq) # Generate rules given stationary frequency vector from an SSE model variant 
  my_rates <- get_rates(my_rules,choose_rule,num_repli,lower_bound,upper_bound,rate_dist)  # Generate random rates from a given rule choice
  #
  
  ###############################
  # 2-REGION GEOSSE
  ###############################
  pdf(file = paste("plot/rate_density_batch_num_",batch_num,".pdf",sep=""),height=7, width=7)
  if (variant == "2-region GeoSSE"){
    #
    n_region    <- 3 # Number of regions
    # Convert string to variable
    param_names <- c("w_A","w_B",
                     "e_A","e_B",
                     "d_A_B","d_B_A",
                     "b_A_B")
    ###############################
    # PLOTTING
    ###############################
    #
    # Title
    title = paste("2-region GeoSSE with pi_A = ",stat_freq[1],", pi_B = ",stat_freq[2],", pi_AB = ",stat_freq[3],", under rule ",choose_rule,sep = "")
    #
    for (i in 1:(2^n_region-2)){
      for (j in (i+1):(2^n_region-1)){
        commonTheme = list(labs(color="Density",fill="Density",
                                x= colnames(my_rates)[i],
                                y= colnames(my_rates)[j],
                                theme_bw(),
                                theme(legend.position=c(0,1),
                                      legend.justification=c(0,1))))
        #
        # browser()
        x_name <- sym(param_names[i])
        y_name <- sym(param_names[j])
        #
        p <- ggplot(data=my_rates,aes_string(x = x_name,y = y_name))
        p <- p + stat_density2d(aes(fill=after_stat(level),alpha=after_stat(level)),geom='polygon',colour='black')
        p <- p + lims(x=c(min(my_rates[,i])-median(my_rates[,i]),max(my_rates[,i])+median(my_rates[,i])),y=c(min(my_rates[,j])-median(my_rates[,j]),max(my_rates[,j])+median(my_rates[,j])))
        p <- p + scale_fill_continuous(low="green",high="red")
        p <- p + guides(alpha="none")
        p <- p + geom_point() 
        p <- p + commonTheme
        p <- p + ggtitle(title)
        p <- p + theme(plot.title = element_text(hjust = 0.5))
        print(p)
      }
    }
  }
  dev.off()
  return(my_rates)
}
