# See the same code in ~/Code Testing/SSA/diffusion_Rpackage/ for complete history of changes made 

####### LOAD PACKAGES ######
############################
library(ggplot2) # Use version 3.4.4. The current version breaks ggtern.
library(ggtern)
library(gridExtra)
library(viridis)

####### FILE SYSTEM ######
##########################
in_fp   = "/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/"
out_fp  = paste0(in_fp,"data/")
plot_fp = paste0(in_fp,"plot/")

setwd(in_fp)

####### INITIAL SETUP ####
##########################
variant   = "2-region GeoSSE"
bin_size  = 0.01

###############################
# 2-REGION GEOSSE MODEL
###############################
#
if(variant == "2-region GeoSSE"){
  ###############################
  # STATIONARY FREQUENCIES SAMPLING
  ###############################
  #
  # Draw pi_A and pi_B from meshgrid with bin of size bin_size
  statio_A = statio_B <- round(seq(0,1,length.out = 1/bin_size + 1),2)
  mesh                <- expand.grid(freq_A=statio_A,freq_B=statio_B)
  mesh <- unique(mesh[,1:2]) # remove duplicated rows from rounding
  mesh$freq_AB        <- round(abs(1- mesh$freq_A - mesh$freq_B),2)
  mesh                <- round(subset(mesh,freq_AB>=0),2) # each combination must sum up to 1
  # Initialize count_rule stored inside mesh 
  mesh$count_rule     <- rep(NA,nrow(mesh))
  #
  # Construct a matrix to contain each stationary vector
  M <- matrix(NA,nrow(mesh),ncol(mesh)-1)
  for (i in 1:ncol(M)){
    for (j in 1:nrow(M)){
      M[j,i] <- mesh[j,i]
    }
  }
  #
  write.table(M,file=paste0(out_fp,"patterns_statio_batch.csv"),sep=",",row.names = FALSE,col.names=FALSE)
  #
  # browser()
  ###############################
  # GETTING RULESETS FROM MATHEMATICA FOR EACH STATIONARY VECTOR AT ONCE 
  ###############################
  #
  # Call Mathematica to return rulesets for each stationary vector from mesh. 
  setwd(paste0(in_fp,"code/"))
  mathematica_wrapper <- "/Applications/Mathematica.app/Contents/MacOS/MathKernel -script \"mathematica_geosse_batch.m\""
  system(mathematica_wrapper)
  #
  rate_name <- c("wa","wb","ea","eb","bab","dab","dba") #list of names of rate parameters in 2-region GeoSSE
  #
  ###############################
  # FILLING IN MESH$COUNT DATA
  ###############################
  setwd(out_fp) # switch directory to read mathematica output
  for (i in 1:nrow(mesh)){
    # Extract rule sets from mesh[i] stationary vector
    filename_from <- 'output_sim_'
    filename_from_upd <- paste(filename_from,i,'.txt',sep = "") #for a set of stationary distributions in matrix row 'row'
    text <- readLines(filename_from_upd, warn = FALSE)
    # Delete the text file containing rules from mesh[j] vector immediately after using 
    # Do not want to keep mathematica output in the file system after loading it to R
    file.remove(paste(filename_from,i,'.txt',sep=""))
    #
    rate_expr <- c() #empty vector to assign rate parameter that has an exact expression
    rate_cond <- c() #empty vector to assign all conditional expressions 
    #replace all occurrences of "ConditionalExpression" with " " 
    text <- gsub("ConditionalExpression","",text)
    #
    #replace all occurences of "Less" with "<","LessEqual" with "<=", "Inequality" with ""
    text <- mgsub(text,c(", Less, ",", LessEqual, ","Inequality"),c("<","<=",""))
    #
    #remove comma and whitespace before inequalities symbol
    text <- mgsub(text,c("<=, ","<, ",", <",", <="),c("<=","<","<","<="))
    #
    # Assign no solution if the given stationary vector has not solution 
    if (length(text)==0){
      # return("The given stationary frequencies do not have any solutions")
      rule_list <- list(rate_expressions = c(), rules_sets = c(), stationary_A = mesh[i,1], stationary_B = mesh[i,2], stationary_AB = mesh[i,3])
    } else { # if it has solution
      #Extract each analytical expression (not including parameters that are conditional expressions represented by inequalities) for each exact parameter
      current_string <- text #current_string, update it in each loop
      for (k in 1:length(rate_name)){
        # browser()
        count_expr <- 0 # token for the while loop
        l <- 1 # token for the while loop
        while (count_expr != 1){
          rate_expr[k] <- sub(paste("\\,",rate_name[k+l], "->.*"), "", current_string) #get the string until (excluding) the next exact parameter 
          count_expr <- str_count(rate_expr[k],"->")
          l <- l+1
        }
        # This code line above has an issue it skips over some elements from rate_name. If it does, you will have two parameters with -> in the same rate expression
        if (nchar(rate_expr[k]) == nchar(current_string)){
          if (grepl(paste(rate_name[k],"->"),current_string)==TRUE){ #check if this character "rate_name[i] ->" is in the current string
            rate_expr[k] <- current_string #check if the current parameter with the exact expression is the last one in the current string
          } else{
            # rate_expr[k] <- paste("no exact expression for parameter",rate_name[i]) #meaning that that particular parameter is a conditional expression
            rate_expr[k] <- NA #meaning that that particular parameter is a conditional expression (does not have exact expression)
          }
        }
        current_string <- gsub(paste(rate_expr[k],",",sep=""),"",current_string,fixed = TRUE) #update current_string without the string that corresponds to the previous parameter
      }
      rate_expr <- rate_expr[!is.na(rate_expr)] #remove all NAs
      for (m in 2:length(rate_expr)){
        # rate_expr[m] <- sub(".*[\\[]([^,]+)[,].*","\\1",rate_expr[m]) #remove conditional expressions from each exact parameters
        rate_expr[m] <- sub(",.*","\\1",rate_expr[m]) #remove conditional expressions from each exact parameters
        rate_expr[m] <- stri_replace_first_fixed(rate_expr[m],"[","") #remove [ at the start of the string of each exact parameter
        #Note: for the first parameter, its conditional expressions will be removed later while collecting the conditional expressions. 
      }
      #
      n_alt_cond <- str_count(rate_expr[1],pattern = '\\|{2,}') + 1 #the number of alternative conditional expressions (rules). Note: Just need to look these from one exact parameter only, since they are repeated in other exact parameters
      for (l in 1:(n_alt_cond-1)){
        cond <- gsub(".*[,]([^||]+)[||].*", "\\1", rate_expr[1]) #find a conditional expression
        rate_cond[l] <- cond 
        rate_expr[1] <- gsub(paste(cond,"||",sep=""),"",rate_expr[1],fixed = TRUE)
      }
      
      # assign the last the conditional expression
      rate_cond[n_alt_cond] <- gsub(".*,", "\\1", rate_expr[1])
      rate_cond[n_alt_cond] <- stri_replace_last_fixed(rate_cond[n_alt_cond],"]","") #remove "]" at the end
      # 
      # browser()
      # rate_expr[1] <- gsub(".*->([^\\[]+)[\\[]","\\1",rate_expr[1]) #remove the letter and -> [ at the start of the string
      rate_expr[1] <- gsub(".*\\{","\\1",rate_expr[1]) #remove { at the start of the string
      # rate_expr[1] <- gsub(".*\\[","\\1",rate_expr[1]) #remove [ at the start of the string
      rate_expr[1] <- stri_replace_first_fixed(rate_expr[1],"[","") #remove [ at the start of the string
      rate_expr[1] <- gsub(",.*","\\1",rate_expr[1]) #remove the conditional expression
      #
      # browser()
      # Return rate expressions along with sets of rules between rate parameters 
      rule_list <- list(rate_expressions = rate_expr, rules_sets = rate_cond, stationary_A = mesh[j,1], stationary_B = mesh[j,2], stationary_AB = mesh[j,3])
      #
      # browser()
      # Add "if" condition when the number of alternative rules is only 1
      if (length(rule_list$rules_sets) > 1){
        for (n in 1:(length(rule_list$rules_sets)-1)){
          rule_list$rules_sets[n] <- gsub("\\[|\\]", "",rule_list$rules_sets[n]) # remove all occurences of "[" and "]" from the condition expressions
          rule_list$rules_sets[n] <- substring(rule_list$rules_sets[n], 3, nchar(rule_list$rules_sets[n])-2) # ignore whitespace and "(" at the beginning and ")" and whitespace at the end.
        }
        rule_list$rules_sets[length(rule_list$rules_sets)] <- gsub("\\[|\\]", "",rule_list$rules_sets[length(rule_list$rules_sets)]) # remove all occurences of "[" and "]" from the condition expressions
        rule_list$rules_sets[length(rule_list$rules_sets)] <- substring(rule_list$rules_sets[length(rule_list$rules_sets)], 3, nchar(rule_list$rules_sets[length(rule_list$rules_sets)])-1) # For last conditional expression only: ignore whitespace and "(" at the beginning and ")"  at the end.
      } else{
        rule_list$rules_sets[1] <- gsub("\\[|\\]", "",rule_list$rules_sets[1]) # remove all occurences of "[" and "]" from the condition expressions
      }
    }
    mesh$count_rule[i] <- length(rule_list$rules_sets)
  }
  #
  mesh$count_rule <- factor(mesh$count_rule)
  ###############################
  # PLOTTING
  ###############################
  ##
  # fig <- ggtern(mesh,aes(x=freq_A,y=freq_B,z=freq_AB)) + geom_hex_tern(binwidth=bin_size,aes(value=count_rule),fun=median) + #the binwidth is chosen according to grid width used for sampling frequencies
  #   scale_fill_viridis(option = "A") + theme_rgbw() + theme_gridsontop() +
  #   labs(title = "Count number of rules given stationary frequencies",x=expression(hat(Pi)[A]), y=expression(hat(Pi)[B]), z=expression(hat(Pi)[AB])) +
  #   Tarrowlab("Stationary frequency B (%)") + Larrowlab("Stationary frequency A (%)") + Rarrowlab("Stationary frequency AB (%)") + theme(plot.title = element_text(hjust = 0.5)) 
  # #
 fig <-  ggtern(mesh, aes(x=freq_A, y=freq_B, z=freq_AB)) +
          geom_point(aes(color = count_rule), size = 0.5) +
          scale_color_manual(values = c("black", "#5B3794", "#FE5C66", "#FFFE9E")) + 
          labs(title = "Count number of rules given stationary frequencies",x=expression(hat(Pi)[A]), y=expression(hat(Pi)[B]), z=expression(hat(Pi)[AB])) +
          Tarrowlab("Stationary frequency B (%)") + Larrowlab("Stationary frequency A (%)") + Rarrowlab("Stationary frequency AB (%)") + theme(plot.title = element_text(hjust = 0.1)) +
          theme_rgbw() 
  #
  fig_count = paste0(plot_fp, "/plot_count.pdf")
  pdf(fig_count, height=4, width=5)
  print(fig)
  dev.off()
}
#
###############################
# BISSE MODEL
###############################
#
if (variant == "BiSSE"){
  ###############################
  # STATIONARY FREQUENCIES SAMPLING
  ###############################
  #
  # Draw pi_A and pi_B from meshgrid with bin of size bin_size
  statio_A            <- round(seq(0,1,length.out = 1/bin_size + 1),2)
  statio_B            <- round(abs(1-statio_A),2)
  mesh                <- data.frame(freq_A = c(statio_A), freq_B = c(statio_B))
  #
  # Initialize count_rule stored inside mesh 
  mesh$count_rule     <- rep(NA,nrow(mesh))
  #
  # Construct a matrix to contain each stationary vector
  M <- matrix(NA,nrow(mesh),ncol(mesh)-1)
  for (i in 1:ncol(M)){
    for (j in 1:nrow(M)){
      M[j,i] <- mesh[j,i]
    }
  }
  #
  write.table(M,file=paste0(out_fp,"patterns_statio_batch.csv"),sep=",",row.names = FALSE,col.names=FALSE)
  #
  ###############################
  # GETTING RULESETS FROM MATHEMATICA FOR EACH STATIONARY VECTOR AT ONCE 
  ###############################
  #
  # Call Mathematica to return rulesets for each stationary vector from mesh. 
  setwd(paste0(in_fp,"code/"))
  mathematica_wrapper <- "/Applications/Mathematica.app/Contents/MacOS/MathKernel -script \"mathematica_bisse_batch.m\""
  system(mathematica_wrapper)
  #
  rate_name <- c("la","lb","ma","mb","qab","qba") #list of names of rate parameters in BiSSE model
  #
  ###############################
  # FILLING IN MESH$COUNT DATA
  ###############################
  setwd(out_fp) # switch directory to read mathematica output
  for (i in 1:nrow(mesh)){
    # Extract rule sets from mesh[i] stationary vector
    filename_from <- 'output_sim_'
    filename_from_upd <- paste(filename_from,i,'.txt',sep = "") #for a set of stationary distributions in matrix row 'row'
    text <- readLines(filename_from_upd, warn = FALSE)
    # Delete the text file containing rules from mesh[j] vector immediately after using 
    # Do not want to keep mathematica output in the file system after loading it to R
    file.remove(paste(filename_from,i,'.txt',sep=""))
    #
    rate_expr <- c() #empty vector to assign rate parameter that has an exact expression
    rate_cond <- c() #empty vector to assign all conditional expressions 
    #replace all occurrences of "ConditionalExpression" with " " 
    text <- gsub("ConditionalExpression","",text)
    #
    #replace all occurences of "Less" with "<","LessEqual" with "<=", "Inequality" with ""
    text <- mgsub(text,c(", Less, ",", LessEqual, ","Inequality"),c("<","<=",""))
    #
    #remove comma and whitespace before inequalities symbol
    text <- mgsub(text,c("<=, ","<, ",", <",", <="),c("<=","<","<","<="))
    #
    # Assign no solution if the given stationary vector has not solution 
    if (length(text)==0){
      # return("The given stationary frequencies do not have any solutions")
      rule_list <- list(rate_expressions = c(), rules_sets = c(), stationary_A = mesh[i,1], stationary_B = mesh[i,2])
    } else { # if it has solution
      #Extract each analytical expression (not including parameters that are conditional expressions represented by inequalities) for each exact parameter
      current_string <- text #current_string, update it in each loop
      for (k in 1:length(rate_name)){
        # browser()
        count_expr <- 0 # token for the while loop
        l <- 1 # token for the while loop
        while (count_expr != 1){
          rate_expr[k] <- sub(paste("\\,",rate_name[k+l], "->.*"), "", current_string) #get the string until (excluding) the next exact parameter 
          count_expr <- str_count(rate_expr[k],"->")
          l <- l+1
        }
        # This code line above has an issue it skips over some elements from rate_name. If it does, you will have two parameters with -> in the same rate expression
        if (nchar(rate_expr[k]) == nchar(current_string)){
          if (grepl(paste(rate_name[k],"->"),current_string)==TRUE){ #check if this character "rate_name[i] ->" is in the current string
            rate_expr[k] <- current_string #check if the current parameter with the exact expression is the last one in the current string
          } else{
            # rate_expr[k] <- paste("no exact expression for parameter",rate_name[i]) #meaning that that particular parameter is a conditional expression
            rate_expr[k] <- NA #meaning that that particular parameter is a conditional expression (does not have exact expression)
          }
        }
        current_string <- gsub(paste(rate_expr[k],",",sep=""),"",current_string,fixed = TRUE) #update current_string without the string that corresponds to the previous parameter
      }
      rate_expr <- rate_expr[!is.na(rate_expr)] #remove all NAs
      for (m in 2:length(rate_expr)){
        # rate_expr[m] <- sub(".*[\\[]([^,]+)[,].*","\\1",rate_expr[m]) #remove conditional expressions from each exact parameters
        rate_expr[m] <- sub(",.*","\\1",rate_expr[m]) #remove conditional expressions from each exact parameters
        rate_expr[m] <- stri_replace_first_fixed(rate_expr[m],"[","") #remove [ at the start of the string of each exact parameter
        #Note: for the first parameter, its conditional expressions will be removed later while collecting the conditional expressions. 
      }
      #
      n_alt_cond <- str_count(rate_expr[1],pattern = '\\|{2,}') + 1 #the number of alternative conditional expressions (rules). Note: Just need to look these from one exact parameter only, since they are repeated in other exact parameters
      for (l in 1:(n_alt_cond-1)){
        cond <- gsub(".*[,]([^||]+)[||].*", "\\1", rate_expr[1]) #find a conditional expression
        rate_cond[l] <- cond 
        rate_expr[1] <- gsub(paste(cond,"||",sep=""),"",rate_expr[1],fixed = TRUE)
      }
      
      # assign the last the conditional expression
      rate_cond[n_alt_cond] <- gsub(".*,", "\\1", rate_expr[1])
      rate_cond[n_alt_cond] <- stri_replace_last_fixed(rate_cond[n_alt_cond],"]","") #remove "]" at the end
      # 
      # browser()
      # rate_expr[1] <- gsub(".*->([^\\[]+)[\\[]","\\1",rate_expr[1]) #remove the letter and -> [ at the start of the string
      rate_expr[1] <- gsub(".*\\{","\\1",rate_expr[1]) #remove { at the start of the string
      # rate_expr[1] <- gsub(".*\\[","\\1",rate_expr[1]) #remove [ at the start of the string
      rate_expr[1] <- stri_replace_first_fixed(rate_expr[1],"[","") #remove [ at the start of the string
      rate_expr[1] <- gsub(",.*","\\1",rate_expr[1]) #remove the conditional expression
      #
      # browser()
      # Return rate expressions along with sets of rules between rate parameters 
      rule_list <- list(rate_expressions = rate_expr, rules_sets = rate_cond, stationary_A = mesh[j,1], stationary_B = mesh[j,2])
      #
      # browser()
      # Add "if" condition when the number of alternative rules is only 1
      if (length(rule_list$rules_sets) > 1){
        for (n in 1:(length(rule_list$rules_sets)-1)){
          rule_list$rules_sets[n] <- gsub("\\[|\\]", "",rule_list$rules_sets[n]) # remove all occurences of "[" and "]" from the condition expressions
          rule_list$rules_sets[n] <- substring(rule_list$rules_sets[n], 3, nchar(rule_list$rules_sets[n])-2) # ignore whitespace and "(" at the beginning and ")" and whitespace at the end.
        }
        rule_list$rules_sets[length(rule_list$rules_sets)] <- gsub("\\[|\\]", "",rule_list$rules_sets[length(rule_list$rules_sets)]) # remove all occurences of "[" and "]" from the condition expressions
        rule_list$rules_sets[length(rule_list$rules_sets)] <- substring(rule_list$rules_sets[length(rule_list$rules_sets)], 3, nchar(rule_list$rules_sets[length(rule_list$rules_sets)])-1) # For last conditional expression only: ignore whitespace and "(" at the beginning and ")"  at the end.
      } else{
        rule_list$rules_sets[1] <- gsub("\\[|\\]", "",rule_list$rules_sets[1]) # remove all occurences of "[" and "]" from the condition expressions
      }
    }
    mesh$count_rule[i] <- length(rule_list$rules_sets)
  }
  ###############################
  # PLOTTING
  ###############################
  #
  colors <- c("green","yellow","red")
  #
  index_lower <- which(mesh$freq_A <= 0.5) # Get all the pairs that correspond to pi_A <= 0.5
  index_upper <- which(mesh$freq_A >= 0.5) # Get all the pairs that correspond to pi_A >= 0.5
  #
  mesh_lower <- mesh[index_lower,]
  mesh_upper <- mesh[index_upper,]
  #
  fig_1 <- ggplot(mesh_lower, aes(x = freq_A, y = freq_B, fill = factor(count_rule))) #factor(count_rule) so it is discrete. 
  fig_1 <- fig_1 + geom_tile()
  fig_1 <- fig_1 + xlab("Stationary frequency A") + ylab("Stationary frequency B")
  fig_1 <- fig_1 + scale_fill_manual(values=colors)
  fig_1 <- fig_1 + labs(fill = "Number of rulesets") 
  #
  fig_2 <- ggplot(mesh_upper, aes(x = freq_A, y = freq_B, fill = factor(count_rule))) 
  fig_2 <- fig_2 + geom_tile()
  fig_2 <- fig_2 + xlab("Stationary frequency A") + ylab("Stationary frequency B")
  fig_2 <- fig_2 + scale_fill_manual(values=colors) 
  fig_2 <- fig_2 + scale_y_reverse() + scale_x_reverse() # Reverse the order of x and y axes to show symmetry
  fig_2 <- fig_2 + labs(fill = "Number of rulesets") 
  #
  grid.arrange(fig_1,fig_2)
}