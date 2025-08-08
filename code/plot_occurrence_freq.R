# Note: in Mathematica, we always assume all rates to be positive (restricted sample space)
# This mathematica is called inside get_rules.R and get_rates.R

# In general, this plot would tell: Under a system that allows for certain rate scenarios (our sample space), 
# how likely the occurrence of this stationary vector compared to the rest. 

# See the same code in ~/Code Testing/SSA/diffusion_Rpackage/ for complete history of changes made 

####### LOAD PACKAGES ######
############################
library(ggplot2) # Use version 3.4.4. The current version breaks ggtern.
library(ggtern)
library(viridis)
library(lhs)
library(mgsub)
library(stringr) #later used for counting the number of alternative conditional expressions
library(stringi)

####### FILE SYSTEM ######
##########################
in_fp   = "/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/"
out_fp  = paste0(in_fp,"data/")
plot_fp = paste0(in_fp,"plot/")

setwd(in_fp)

####### INITIAL SETUP ######
##########################
variant         = "2-region GeoSSE"
bin_size        = 0.01
num_trial       = 1000
lower_bound_vec = c(0.01,0.01,0.01,0.01,0.01,0.01,0.01)
upper_bound_vec = c(1,1,1,1,1,1,1)
tol             = 50

# For "2-region GeoSSE" model
if(variant == "2-region GeoSSE"){
  ############# ERROR CHECKING ################
  #############################################
  # Error checking if lower_bound and upper_bound for rates agree with our diffusion model's assumption (in Mathematica, all rates must be positive)
  param_names <- c("wa","wb","ea","eb","dab","dba","bab")
  #
  # Error checking for the length of the bound vectors 
  if (length(lower_bound_vec) != 7 || length(upper_bound_vec) != 7){
    stop("lower and upper bound vectors must be of length 7")
  }
  # Error checking if any of the rates take negative or 0 values
  for (i in 1:7){
    if (lower_bound_vec[i] < 0 || lower_bound_vec[i] == 0){
      stop(paste(param_names[i]," cannot take negative values or 0",sep=""))
    }
    # Error checking if upper_bound <= lower_bound, return error 
    if (lower_bound_vec[i] >= upper_bound_vec[i]){
      stop("Lower bound must be smaller than upper bound")
    }
  }
  #
  ############# RATE PARAMETER SAMPLING ################
  ######################################################
  # Sampling using latin hypercube sampling method (assuming independence of rate parameters)
  LHS_sample <- randomLHS(num_trial,7)
  for (i in 1:7){
    # for each parameter, LHS sampling from [lower_bound,upper_bound]
    LHS_sample[,i] <- lower_bound_vec[i] + (upper_bound_vec[i]-lower_bound_vec[i])*LHS_sample[,i]
  }
  LHS_sample <- round(LHS_sample,2)
  rate_sample <- as.data.frame(LHS_sample)
  colnames(rate_sample) <- c('wa','wb','ea','eb','dab','dba','bab')
  # # Generate random sample of size num_trial for rate combinations 
  # rate_sample = data.frame(wa   = runif(num_trial,lower_bound_vec[1],upper_bound_vec[1]), 
  #                          wb   = runif(num_trial,lower_bound_vec[2],upper_bound_vec[2]), 
  #                          ea   = runif(num_trial,lower_bound_vec[3],upper_bound_vec[3]), 
  #                          eb   = runif(num_trial,lower_bound_vec[4],upper_bound_vec[4]),
  #                          dab = runif(num_trial,lower_bound_vec[5],upper_bound_vec[5]), 
  #                          dba = runif(num_trial,lower_bound_vec[6],upper_bound_vec[6]), 
  #                          bab = runif(num_trial,lower_bound_vec[7],upper_bound_vec[7]))
  #
  ############# STATIONARY FREQUENCIES SAMPLING ################
  #######################################################################
  # Draw pi_A and pi_B from meshgrid with bin of size bin_size
  # 
  statio_A=statio_B   <- round(seq(0,1,length.out = 1/bin_size + 1),2)
  mesh                <- expand.grid(freq_A=statio_A,freq_B=statio_B)
  mesh <- unique(mesh[,1:2]) # remove duplicated rows from rounding
  mesh$freq_AB        <- round(abs(1- mesh$freq_A - mesh$freq_B),2)
  mesh                <- round(subset(mesh,freq_AB>=0),2) # each combination must sum up to 1
  #
  # Remove all the symmetrical pairs. Later they will be projected into mesh_project grid and copied the density from their
  # respective pairs
  # browser()
  mesh_length      <- nrow(mesh)
  for (i in 1:mesh_length){
    find_sym_A <- which(mesh[,2]==mesh[i,1]) 
    find_sym_B <- which(mesh[,1]==mesh[i,2])
    find_ind   <- intersect(find_sym_A,find_sym_B)
    # The second condition is needed so it does not remove the case where \pi_A=\pi_B
    if (length(find_ind == 1) && mesh[i,1] != mesh[i,2]){
      mesh       <-  mesh[-find_ind,]
    }
  }
  #
  ### Create the projected mesh
  mesh_project <- data.frame(freq_A = NA,freq_B = NA, freq_AB = NA,tot_prob = NA)
  # mesh_subset <- subset(mesh,mesh$freq_A<=0.5) # Collect all grid cells corresponding to region where \Pi_A <= 0.5
  
  #
  # Construct a matrix to contain each stationary vector
  M <- matrix(NA,nrow(mesh),ncol(mesh))
  for (i in 1:ncol(M)){
    for (j in 1:nrow(M)){
      M[j,i] <- mesh[j,i]
    }
  }
  #
  write.table(M,file=paste0(out_fp,"patterns_statio_batch.csv"),sep=",",row.names = FALSE,col.names=FALSE)
  #
  # browser()
  ############# GETTING RULESETS FROM MATHEMATICA FOR EACH STATIONARY VECTOR AT ONCE ################
  #############################################################################################
  # Call Mathematica to return rulesets for each stationary vector from mesh. 
  setwd(paste0(in_fp,"code/"))
  mathematica_wrapper <- "/Applications/Mathematica.app/Contents/MacOS/MathKernel -script \"mathematica_geosse_batch.m\""
  system(mathematica_wrapper)
  #
  # Initialize total_prob inside mesh 
  # (this is the total probability that describes that the given rates sample is more likely to be explained
  # by a specific stationary vector governed by rules)
  # I.e. whether that specific stationary vector is a common or rare occurence, under certain restrictions on range of rate values.
  mesh$tot_prob         <- rep(NA,nrow(mesh))
  mesh_length_trim      <- nrow(mesh)
  # browser()
  ############# FORMAT STRINGS FROM MATHEMATICA ################
  ##############################################################
  setwd(out_fp) # switch directory to read mathematica output
  # Looping over each combination of different stationary frequency vectors 
  for (j in 1:mesh_length_trim){ #looping over each stationary vector 
    # Extract rule sets from mesh[j] stationary vector
    filename_from <- 'output_sim_'
    filename_from_upd <- paste(filename_from,j,'.txt',sep = "") #for a set of stationary distributions in matrix row 'row'
    text <- readLines(filename_from_upd, warn = FALSE)
    # Delete the text file containing rules from mesh[j] vector immediately after using 
    # Do not want to keep mathematica output in the file system after loading it to R
    file.remove(paste(filename_from,j,'.txt',sep=""))
    #
    rate_name <- c("wa","wb","ea","eb","bab","dab","dba") #list of names of rate parameters in 2-region GeoSSE
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
      rule_list <- list(rate_expressions = c(), rules_sets = c(), stationary_A = mesh[j,1], stationary_B = mesh[j,2], stationary_AB = mesh[j,3])
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
    ############# SEPARATE INEQUALITIES AND EQUALITIES IN A RULE SET FROM A STATIONARY VECTOR ################
    ##########################################################################################################
    # Count the number of rule sets 
    total_prob <- 0 # Initialize total probability across rules sets for that stationary vector 
    n_rules   <- length(rule_list$rules_sets) # Get number of rules sets for that given stationary vector
    #
    ############# IN CASE THAT STATIONARY VECTOR DOES NOT HAVE A RULE SET ################
    ######################################################################################
    if (n_rules == 0){ # if there is no rule set for that specific stationary vector (under the condition that all rates must be positive)
      mesh$tot_prob[j] <- total_prob # Then, that stationary vector is impossible under rates' assumption dictated by lower_bound and upper_bound
      #### Assigning the same density to its symmetric pair
      mesh_project <- rbind(mesh_project,c(mesh$freq_B[j],mesh$freq_A[j],mesh$freq_AB[j],mesh$tot_prob[j]))
      #
      ############# IN CASE THAT STATIONARY VECTOR DOES have RULE SET(S) ###################
      ######################################################################################
    } else { # if there are rule set(s), then we loop over each rule set to compute total_prob
      #
      for (o in 1:n_rules){ #looping over each ruleset derived from a drawn set of stationary frequencies
        # Compute success probability for a rule from a stationary vector 
        # Separate each conditional expression from a given rule set and store them inside cond_string 
        no_cond_expr <- str_count(rule_list$rules_sets[o],"&&") # Find the number of conditional expressions from a particular rule (from o)
        #
        cond_string <- rep(NA,no_cond_expr) # Vector that will contain each inequality from the choice of rule set
        cond_expr <- strsplit(rule_list$rules_sets[o], "&&") # separate each conditional expression from the chosen rule set 
        for(p in 1:(no_cond_expr+1)){
          cond_string[p] <- cond_expr[[1]][p]
        }
        #
        # Separate inequalities in a conditional expression from a chosen rule set (e.g "a < b < c" to "a < b" and "b < c")
        # NOTE: we do this, so we can use parse() to turn into math expressions.
        # browser()
        inequal_string <- c()
        for (q in 1:length(cond_string)){ # loop over all conditional expressions from the chosen rule set 
          total_char <- nchar(cond_string[q]) # total number of characters from that conditional expression 
          if (str_count(cond_string[q],"<") == 2 && str_count(cond_string[q],"<=") == 0){ # case a: only 2 "<" in cond_string[q]
            # gsub("^([^<]*<[^<]*)<.*$", "\\1", cond_string[1])
            inequa_post <- gregexpr(pattern ='<',cond_string[q]) # positions of "<' in that conditional expression
            # Getting left_substring ... < (<=) ..
            left_string <- substring(cond_string[q],1,inequa_post[[1]][2]-1)
            # Append the substring
            inequal_string <- append(inequal_string,left_string) 
            # Getting right_substring 
            right_string <- substring(cond_string[q],inequa_post[[1]][1]+1,total_char)
            inequal_string <- append(inequal_string,right_string)
            #
          } else if (str_count(cond_string[q],">") == 2 && str_count(cond_string[q],">=") == 0){ # case b: only 2 ">" in cond_string[q]
            inequa_post <- gregexpr(pattern ='>',cond_string[q]) # positions of ">' in that conditional expression
            # Getting left_substring ... > (>=) ..
            left_string <- substring(cond_string[q],1,inequa_post[[1]][2]-1)
            # Append the substring 
            inequal_string <- append(inequal_string,left_string)
            # Getting right_substring 
            right_string <- substring(cond_string[q],inequa_post[[1]][1]+1,total_char)
            inequal_string <- append(inequal_string,right_string)
            #
          } else if (str_count(cond_string[q],"<") == 2 && str_count(cond_string[q],"<=") == 1){ # case c: only 1 "<"  and 1 "<=" in cond_string[q]
            inequa_post <- gregexpr(pattern ='<',cond_string[q]) # positions of "<' in that conditional expression
            # Getting left_substring ... < (<=) ..
            left_string <- substring(cond_string[q],1,inequa_post[[1]][2]-1)
            # Append the substring
            inequal_string <- append(inequal_string,left_string) 
            # Getting right_substring 
            right_string <- substring(cond_string[q],inequa_post[[1]][1]+2,total_char)
            if (str_count(right_string,"=")==0){ # This means "<=" is at the leftmost inequality in cond_string[q]
              inequal_string <- append(inequal_string,right_string)
            } else { # This means "<=" is at the rightmost inequality in cond_string[q]
              right_string <- substring(cond_string[q],inequa_post[[1]][1]+1,total_char)
              inequal_string <- append(inequal_string,right_string)
            }
            #  
          } else if (str_count(cond_string[q],"<") == 2 && str_count(cond_string[q],"<=") == 2){ # case d:  2 "<=" in cond_string[q]
            inequa_post <- gregexpr(pattern ='<',cond_string[q]) # positions of "<' in that conditional expression
            # Getting left_substring ... < (<=) ..
            left_string <- substring(cond_string[q],1,inequa_post[[1]][2]-1)
            # Append the substring
            inequal_string <- append(inequal_string,left_string) 
            # Getting right_substring 
            right_string <- substring(cond_string[q],inequa_post[[1]][1]+2,total_char)
            inequal_string <- append(inequal_string,right_string)
            #  
          } else if (str_count(cond_string[q],">") == 2 && str_count(cond_string[q],">=") == 1){# case e: only 1 ">"  and 1 ">=" in cond_string[q]
            inequa_post <- gregexpr(pattern ='>',cond_string[q]) # positions of ">' in that conditional expression
            # Getting left_substring ... > (>=) ..
            left_string <- substring(cond_string[q],1,inequa_post[[1]][2]-1)
            # Append the substring 
            inequal_string <- append(inequal_string,left_string)
            # Getting right_substring 
            right_string <- substring(cond_string[q],inequa_post[[1]][1]+2,total_char)
            if (str_count(right_string,"=")==0){ # This means ">=" is at the leftmost inequality in cond_string[q]
              inequal_string <- append(inequal_string,right_string)
            } else { # This means ">=" is at the rightmost inequality in cond_string[q]
              right_string <- substring(cond_string[q],inequa_post[[1]][1]+1,total_char)
              inequal_string <- append(inequal_string,right_string)
            }
          } else if (str_count(cond_string[q],">") == 2 && str_count(cond_string[q],">=") == 2){# case f: 2 ">=" in cond_string[q]
            inequa_post <- gregexpr(pattern ='>',cond_string[q]) # positions of ">' in that conditional expression
            # Getting left_substring ... > (>=) ..
            left_string <- substring(cond_string[q],1,inequa_post[[1]][2]-1)
            # Append the substring 
            inequal_string <- append(inequal_string,left_string)
            # Getting right_substring 
            right_string <- substring(cond_string[q],inequa_post[[1]][1]+2,total_char)
            inequal_string <- append(inequal_string,right_string)
          } else if (str_count(cond_string[q],">") == 1 || str_count(cond_string[q],"<") == 1){ # case g: 1 "> (>=)" or 1 "<(<=)"
            inequal_string <- append(inequal_string,cond_string[q])
          }
        } 
        #
        #### Get the equalities from rule_list$rate_expressions
        n_equalities <- length(rule_list$rate_expressions) # Number of rate equalities 
        # Label which rate that has exact expression e.g. wa = ...
        label_exact <- c()
        #
        for (r in 1:n_equalities){
          rule_list$rate_expressions[r] <- gsub("->","=",rule_list$rate_expressions[r]) # Replace all "->" with "="
          label_exact <- append(label_exact,strsplit(rule_list$rate_expressions[r]," =")[[1]][1])
        }
        label_exact <- str_trim(label_exact) #remove the whitespace
        # Find the position of these variables from the rate_sample's column names
        index_exact <- match(label_exact,names(rate_sample))
        #
        # browser()
        ############# COMPUTE THE PROBABILITIES FOR EACH STATIONARY VECTOR ################
        ###################################################################################
        # Initiate the counting tokens for the number of fails and number of successes
        fail_count <- 0
        success_count <- 0
        # For the first replicate
        wa  <- rate_sample[1,1]
        wb  <- rate_sample[1,2]
        ea  <- rate_sample[1,3]
        eb  <- rate_sample[1,4]
        dab <- rate_sample[1,5]
        dba <- rate_sample[1,6]
        bab <- rate_sample[1,7]
        #
        # browser()
        true_or_false <- c()
        for (r in 1:length(inequal_string)){ #Evaluate for each inequalities from a given ruleset
          if(eval(parse(text=inequal_string[r]))==TRUE){
            true_or_false <- append(true_or_false,TRUE)
          } else {
            true_or_false <- append(true_or_false,FALSE)
          }
        }
        if (all(true_or_false)== TRUE){ #if the rate parameters draw satsifies the inequalities
          # Check if the rate combination also satisfies the equalities from rate_expressions
          true_or_false_equal <- c()
          for (s in 1:length(index_exact)){
            val_mathe <- round(eval(parse(text=rule_list$rate_expressions[s])),2) #Value from the equality in s from mathematica
            # if ( abs(val_mathe - rate_sample[1,index_exact[s]]) <= tol){
            if ( abs((val_mathe - rate_sample[1,index_exact[s]])/val_mathe)*100 <= tol){
              true_or_false_equal <- append(true_or_false_equal,TRUE)
            } else{
              true_or_false_equal <- append(true_or_false_equal,FALSE)
            }
          }
          if (all(true_or_false_equal)==TRUE){
            fail_count = fail_count  # if the draw satsifies the equalities
            success_count = success_count + 1
          } else {
            fail_count = fail_count + 1  # if the draw does not satisfy the equalities 
            success_count = success_count
          }
        } else {
          fail_count = fail_count + 1  # if it does not satisfy the inequalities
          success_count = success_count 
        }
        # Model parameters in a 2-region GeoSSE system 
        # browser()
        # Check if these random parameter values satisfy all inequalities (stored in inequal_string) for the chosen rule set. 
        while(fail_count + success_count < num_trial){ #continue the checking for sampled parameters #2 to num_trial
          for (t in 2:num_trial){ # For sample replicate 2 and second from last
            # Initialize true and false string 
            true_or_false <- c()
            # Get the rate combination for replicate n
            wa  <- rate_sample[t,1] 
            wb  <- rate_sample[t,2]
            ea  <- rate_sample[t,3]
            eb  <- rate_sample[t,4]
            dab <- rate_sample[t,5]
            dba <- rate_sample[t,6]
            bab <- rate_sample[t,7]
            # Check if the given rate combination satsify the inequalities 
            for (u in 1:length(inequal_string)){
              if(eval(parse(text=inequal_string[u]))==TRUE){
                true_or_false <- append(true_or_false,TRUE)
              } else {
                true_or_false <- append(true_or_false,FALSE)
              }
            }
            #
            if (all(true_or_false)==TRUE){# If it does satisfy the inequalities 
              # Check if the rate combination also satisfies the equalities from rate_expressions
              true_or_false_equal <- c()
              for (v in 1:length(index_exact)){
                val_mathe <- round(eval(parse(text=rule_list$rate_expressions[v])),2) #Value from the equality in v from mathematica
                # if (abs(val_mathe - rate_sample[t,index_exact[v]])<= tol){
                if (abs((val_mathe - rate_sample[t,index_exact[v]])/val_mathe)*100 <= tol){
                  true_or_false_equal <- append(true_or_false_equal,TRUE)
                } else{
                  true_or_false_equal <- append(true_or_false_equal,FALSE)
                }
              }
              if (all(true_or_false_equal)==TRUE){
                fail_count = fail_count  # if the draw satsifies the equalities
                success_count = success_count + 1
              } else {
                fail_count = fail_count + 1  # if the draw does not satisfy the equalities 
                success_count = success_count
              }
            } else { # If it does not satisfy the inequalities 
              fail_count    <- fail_count + 1
              success_count <- success_count
            } 
            #
          }
        }
        success_prob <- success_count/num_trial
        # Sum the probabilities over all rules
        total_prob <- total_prob + success_prob
      }
      mesh$tot_prob[j] <- total_prob # Assigning total probability for j'th stationary vector
      mesh_project <- rbind(mesh_project,c(mesh$freq_B[j],mesh$freq_A[j],mesh$freq_AB[j],mesh$tot_prob[j]))
      #
    }
    #### Assigning the same density to its symmetric pair
    # browser()
    print(paste("Done computing for stationary vector ",j,sep=""))
  }
  # browser()
  mesh_project <- mesh_project[2:nrow(mesh_project),] # Remove the first row that contains NA values
  # Append the two mesh grids 
  mesh <- rbind(mesh,mesh_project)
  #
  mesh <- unique(mesh[,1:4]) # remove duplicated rows from projection (for the case pi_A=pi_B = 0.5)
  ############# PLOTTING ################
  #############################################################################################
  fig <- ggtern(mesh,aes(x=freq_A,y=freq_B,z=freq_AB)) + geom_hex_tern(binwidth=bin_size,aes(value=tot_prob),fun=median) + #the bidwidth is chosen according to grid width used for sampling frequencies
    scale_fill_viridis(option = "A") + theme_rgbw() + theme_gridsontop() +
    labs(title = paste("Relative occurrence of stationary frequency vector given rate samples of size ",num_trial," where \n",
                       param_names[1]," in ", "[",lower_bound_vec[1],", ",upper_bound_vec[1],"], ",
                       param_names[2]," in ", "[",lower_bound_vec[2],", ",upper_bound_vec[2],"], ",
                       param_names[3]," in ", "[",lower_bound_vec[3],", ",upper_bound_vec[3],"], ",
                       param_names[4]," in ", "[",lower_bound_vec[4],", ",upper_bound_vec[4],"], \n",
                       param_names[5]," in ", "[",lower_bound_vec[5],", ",upper_bound_vec[5],"], ",
                       param_names[6]," in ", "[",lower_bound_vec[6],", ",upper_bound_vec[6],"], ",
                       param_names[7]," in ", "[",lower_bound_vec[7],", ",upper_bound_vec[7],"], ","with tolerance = ",tol,
                       sep=""),
         x=expression(hat(Pi)[A]), y=expression(hat(Pi)[B]), z=expression(hat(Pi)[AB])) +
    Tarrowlab("Stationary frequency B (%)") + Larrowlab("Stationary frequency A (%)") + Rarrowlab("Stationary frequency AB (%)") + theme(plot.title = element_text(hjust = 0.5)) 
  #
  suppressWarnings(print(fig))
}

fig_occurence = paste0(plot_fp, "/plot_occurence.pdf")
print(fig_occurence)
pdf(fig_occurence, height=7, width=10)
print(fig)
dev.off()