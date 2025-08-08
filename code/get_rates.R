# This function draws random rate values for a given set(s) of rules corresponding to a set of stationary frequencies for a 2-region GeoSSE model

# See the same code in ~/Code Testing/SSA/diffusion_Rpackage/ for complete history of changes made 

get_rates <- function(rules_sets_list,choose_rule,num_repli,lower_bound,upper_bound,rate_dist){
  #
  ############## WORKING DIRECTORY (TEMPORARY, REMOVE ONCE PORTED TO PACKAGE) ##############
  ##########################################################################################
  setwd('/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/') 
  #
  ############## CHECK WHICH SSE MODEL THAT CORRESPONDS TO THE GIVEN SETS OF RULES ##############
  ###############################################################################################
  if (names(rules_sets_list)[3] == "stationary_A" && names(rules_sets_list)[4] == "stationary_B" && names(rules_sets_list)[5] == "stationary_AB"){
    #
    # Determine the number of alternative rules for that given stationary frequency vector
    number_rules <- length(rules_sets_list$rules_sets)
    #
    # Error checking: return an error message if the choice or rule set exceeds the number of alternative rule sets from sim_no R.Data
    if (choose_rule > number_rules){
      error_rules <- "There are only "
      stop(paste(error_rules,number_rules," alternative rule sets for that set of stationary frequency",sep = ""))
    }
    # Create data frame containing rate values drawn for num_repli number of replicates
    empt_repli <- rep(NA,num_repli) # according to number of replicate draws that user sets.
    rand_draw <- data.frame(w_A = empt_repli, w_B = empt_repli, e_A = empt_repli, e_B = empt_repli, d_A_B = empt_repli, d_B_A = empt_repli, b_A_B = empt_repli, which_rule = empt_repli,  
                            stationary_A = rep(rules_sets_list$stationary_A,num_repli), stationary_B = rep(rules_sets_list$stationary_B,num_repli), stationary_AB = rep(rules_sets_list$stationary_AB,num_repli))
    
    
    # Separate each conditional expression from a given rule set and store them inside cond_string 
    no_cond_expr <- str_count(rules_sets_list$rules_sets[choose_rule],"&&") # Find the number of conditional expressions from a particular rule (from choose_rule)
    #
    cond_string <- rep(NA,no_cond_expr) # Vector that will contain each inequality from the choice of rule set
    cond_expr <- strsplit(rules_sets_list$rules_sets[choose_rule], "&&") # separate each conditional expression from the chosen rule set 
    for(i in 1:(no_cond_expr+1)){
      cond_string[i] <- cond_expr[[1]][i]
    }
    
    # Separate inequalities in a conditional expression from a chosen rule set (e.g "a < b < c" to "a < b" and "b < c")
    # NOTE: we do this, so we can use parse() to turn into math expressions.
    # browser()
    inequal_string <- c()
    for (j in 1:length(cond_string)){ # loop over all conditional expressions from the chosen rule set 
      total_char <- nchar(cond_string[j]) # total number of characters from that conditional expression 
      if (str_count(cond_string[j],"<") == 2 && str_count(cond_string[j],"<=") == 0){ # case a: only 2 "<" in cond_string[j]
        # gsub("^([^<]*<[^<]*)<.*$", "\\1", cond_string[1])
        inequa_post <- gregexpr(pattern ='<',cond_string[j]) # positions of "<' in that conditional expression
        # Getting left_substring ... < (<=) ..
        left_string <- substring(cond_string[j],1,inequa_post[[1]][2]-1)
        # Append the substring
        inequal_string <- append(inequal_string,left_string) 
        # Getting right_substring 
        right_string <- substring(cond_string[j],inequa_post[[1]][1]+1,total_char)
        inequal_string <- append(inequal_string,right_string)
        #
      } else if (str_count(cond_string[j],">") == 2 && str_count(cond_string[j],">=") == 0){ # case b: only 2 ">" in cond_string[j]
        inequa_post <- gregexpr(pattern ='>',cond_string[j]) # positions of ">' in that conditional expression
        # Getting left_substring ... > (>=) ..
        left_string <- substring(cond_string[j],1,inequa_post[[1]][2]-1)
        # Append the substring 
        inequal_string <- append(inequal_string,left_string)
        # Getting right_substring 
        right_string <- substring(cond_string[j],inequa_post[[1]][1]+1,total_char)
        inequal_string <- append(inequal_string,right_string)
        #
      } else if (str_count(cond_string[j],"<") == 2 && str_count(cond_string[j],"<=") == 1){ # case c: only 1 "<"  and 1 "<=" in cond_string[j]
        inequa_post <- gregexpr(pattern ='<',cond_string[j]) # positions of "<' in that conditional expression
        # Getting left_substring ... < (<=) ..
        left_string <- substring(cond_string[j],1,inequa_post[[1]][2]-1)
        # Append the substring
        inequal_string <- append(inequal_string,left_string) 
        # Getting right_substring 
        right_string <- substring(cond_string[j],inequa_post[[1]][1]+2,total_char)
        if (str_count(right_string,"=")==0){ # This means "<=" is at the leftmost inequality in cond_string[j]
          inequal_string <- append(inequal_string,right_string)
        } else { # This means "<=" is at the rightmost inequality in cond_string[j]
          right_string <- substring(cond_string[j],inequa_post[[1]][1]+1,total_char)
          inequal_string <- append(inequal_string,right_string)
        }
        #  
      } else if (str_count(cond_string[j],"<") == 2 && str_count(cond_string[j],"<=") == 2){ # case d:  2 "<=" in cond_string[j]
        inequa_post <- gregexpr(pattern ='<',cond_string[j]) # positions of "<' in that conditional expression
        # Getting left_substring ... < (<=) ..
        left_string <- substring(cond_string[j],1,inequa_post[[1]][2]-1)
        # Append the substring
        inequal_string <- append(inequal_string,left_string) 
        # Getting right_substring 
        right_string <- substring(cond_string[j],inequa_post[[1]][1]+2,total_char)
        inequal_string <- append(inequal_string,right_string)
        #  
      } else if (str_count(cond_string[j],">") == 2 && str_count(cond_string[j],">=") == 1){# case e: only 1 ">"  and 1 ">=" in cond_string[j]
        inequa_post <- gregexpr(pattern ='>',cond_string[j]) # positions of ">' in that conditional expression
        # Getting left_substring ... > (>=) ..
        left_string <- substring(cond_string[j],1,inequa_post[[1]][2]-1)
        # Append the substring 
        inequal_string <- append(inequal_string,left_string)
        # Getting right_substring 
        right_string <- substring(cond_string[j],inequa_post[[1]][1]+2,total_char)
        if (str_count(right_string,"=")==0){ # This means ">=" is at the leftmost inequality in cond_string[j]
          inequal_string <- append(inequal_string,right_string)
        } else { # This means ">=" is at the rightmost inequality in cond_string[j]
          right_string <- substring(cond_string[j],inequa_post[[1]][1]+1,total_char)
          inequal_string <- append(inequal_string,right_string)
        }
      } else if (str_count(cond_string[j],">") == 2 && str_count(cond_string[j],">=") == 2){# case f: 2 ">=" in cond_string[j]
        inequa_post <- gregexpr(pattern ='>',cond_string[j]) # positions of ">' in that conditional expression
        # Getting left_substring ... > (>=) ..
        left_string <- substring(cond_string[j],1,inequa_post[[1]][2]-1)
        # Append the substring 
        inequal_string <- append(inequal_string,left_string)
        # Getting right_substring 
        right_string <- substring(cond_string[j],inequa_post[[1]][1]+2,total_char)
        inequal_string <- append(inequal_string,right_string)
      } else if (str_count(cond_string[j],">") == 1 || str_count(cond_string[j],"<") == 1){ # case g: 1 "> (>=)" or 1 "<(<=)"
        inequal_string <- append(inequal_string,cond_string[j])
      }
    } 
    ############## GET THE EQUALITIES FROM RULES_SETS_LIST$RATE_EXPRESSIONS ##############
    ######################################################################################
    n_equalities <- length(rules_sets_list$rate_expressions) # Number of rate equalities 
    param_name   <- c("wa","wb","ea","eb","dab","dba","bab")
    # Label which rate that has exact expression e.g. wa = ...
    label_exact <- c()
    #
    for (r in 1:n_equalities){
      rules_sets_list$rate_expressions[r] <- gsub("->","=",rules_sets_list$rate_expressions[r]) # Replace all "->" with "="
      label_exact <- append(label_exact,strsplit(rules_sets_list$rate_expressions[r]," =")[[1]][1])
    }
    label_exact <- str_trim(label_exact) #remove the whitespace
    # Find the position of these variables from the rate_sample's column names
    index_exact <- match(label_exact,param_name)
    #
    ############## REJECTION SAMPLING TO SELECT THE RATES ##############
    ####################################################################
    for (k in 1:num_repli){
      # Model parameters in a 2-region GeoSSE system 
      # browser()
      if (rate_dist == "uniform"){
        # Check if these random parameter values satisfy all inequalities (stored in inequal_string) for the chosen rule set. 
        true_or_false <- FALSE
        # browser()
        while (all(true_or_false) == FALSE){
          # Generate another draw for the parameter values 
          if (all(true_or_false) == FALSE){
            # Sample uniformly from [lower_bound,upper_bound]
            wa  <- round(runif(1,min=lower_bound,max=upper_bound),2)
            wb  <- round(runif(1,min=lower_bound,max=upper_bound),2)
            ea  <- round(runif(1,min=lower_bound,max=upper_bound),2)
            eb  <- round(runif(1,min=lower_bound,max=upper_bound),2)
            dab <- round(runif(1,min=lower_bound,max=upper_bound),2)
            dba <- round(runif(1,min=lower_bound,max=upper_bound),2)
            bab <- round(runif(1,min=lower_bound,max=upper_bound),2)
            # # Sample uniformly from [0.01,10]
            # wa  <- round(runif(1,min=0.01,max=10),2)
            # wb  <- round(runif(1,min=0.01,max=10),2)
            # ea  <- round(runif(1,min=0.01,max=10),2)
            # eb  <- round(runif(1,min=0.01,max=10),2)
            # dab <- round(runif(1,min=0.01,max=10),2)
            # dba <- round(runif(1,min=0.01,max=10),2)
            # bab <- round(runif(1,min=0.01,max=10),2)
            #
            # # Sample uniformly from [0.01,1]
            # wa  <- round(runif(1,min=0,max=1),2)
            # wb  <- round(runif(1,min=0,max=1),2)
            # ea  <- round(runif(1,min=0,max=1),2)
            # eb  <- round(runif(1,min=0,max=1),2)
            # dab <- round(runif(1,min=0,max=1),2)
            # dba <- round(runif(1,min=0,max=1),2)
            # bab <- round(runif(1,min=0,max=1),2)
            #
            # # Sample uniformly from [0.01,10]
            # wa  <- round(runif(1,min=0.01,max=1.01),2)
            # wb  <- round(runif(1,min=0.01,max=1.01),2)
            # ea  <- round(runif(1,min=0.01,max=1.01),2)
            # eb  <- round(runif(1,min=0.01,max=1.01),2)
            # dab <- round(runif(1,min=8.01,max=9.01),2)
            # dba <- round(runif(1,min=8.01,max=9.01),2)
            # bab <- round(runif(1,min=0.01,max=1.01),2)
          }
          #
          true_or_false <- c()
          for (l in 1:length(inequal_string)){
            if(eval(parse(text=inequal_string[l]))==TRUE){
              true_or_false <- append(true_or_false,TRUE)
            } else {
              true_or_false <- append(true_or_false,FALSE)
            }
          }
        }
      }
      # browser()
      # Assign the new parameter values that satisfy conditions in inequal_string for a given draw replicate
      rand_draw$w_A[k]    <- wa
      rand_draw$w_B[k]    <- wb
      rand_draw$e_A[k]    <- ea
      rand_draw$e_B[k]    <- eb
      rand_draw$d_A_B[k]  <- dab
      rand_draw$d_B_A[k]  <- dba
      rand_draw$b_A_B[k]  <- bab
      ############## ASSIGN THE VALUES OF THE EXACT EXPRESSIONS FROM INDEX_EXACT ##############
      #########################################################################################
      for (m in 1:length(index_exact)){
        rand_draw[k,index_exact[m]] <- round(eval(parse(text=rules_sets_list$rate_expressions[m])),2)
      }
      #
      rand_draw$which_rule[k] <- choose_rule
    }
  }
  return(rand_draw)
}