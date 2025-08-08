# This function returns sets of rules describing rate relationships in a SSE model variant as a list, given a set of stationary frequencies 

# See the same code in ~/Code Testing/SSA/diffusion_Rpackage/ for complete history of changes made 

# Current issue (not solved) 23 March 2024
# If we change the condition in Mathematica so that parameters can take 0 value
# Then it will return more solutions to rate_expr. So, we might have two different expression for w_A (as an example)
# So, to be able to run this in get_rule, I need to loop over rate_expr as well, not just rate_cond


#load the required packages
get_rules <- function(variant,stat_freq){
  
  ####### LOAD PACKAGES ######
  ############################
  library(mgsub)
  library(stringr) #later used for counting the number of alternative conditional expressions
  library(stringi)
  
  ####### FILE SYSTEM ######
  ##########################
  in_fp   = "/Users/albertsoewongsono/Documents/Code Testing/rate_pattern_diffusion_SSE/"
  out_fp  = paste0(in_fp,"data/")
  plot_fp = paste0(in_fp,"plot/")
  
  setwd(in_fp)
  
  # Model choices
  if (variant == "2-region GeoSSE"){
    # list of rate parameters for 2-region GeoSSE (the same ordering as in Mathematica)
    rate_name <- c("wa","wb","ea","eb","bab","dab","dba") #list of names of rate parameters in 2-region GeoSSE
    rate_expr <- c() #empty vector to assign rate parameter that has an exact expression
    rate_cond <- c() #empty vector to assign all conditional expressions 
    # Note: the length of rate_cond implies the number of rules for a given vector of stationary distributions
    
    # Error checking 
    if (length(stat_freq) != 3){
      stop("The length of stationary frequency vector must be equal to 3 for a 2-region GeoSSE model")
    }
    
    # Construct matrix M that stores stationary frequencies, saved, and called in Mathematica
    M <- matrix(NA,1,3) 
    M[1,1] <- stat_freq[1] # Pi_A
    M[1,2] <- stat_freq[2] # Pi_B
    M[1,3] <- stat_freq[3] # Pi_AB
    # Save the matrix M
    write.table(M,file=paste0(out_fp,"patterns_statio_batch.csv"),sep=",",row.names = FALSE,col.names=FALSE)
    
    # Call Mathematica for a given M matrix
    setwd(paste0(in_fp,"code/"))
    mathematica_wrapper <- "/Applications/Mathematica.app/Contents/MacOS/MathKernel -script \"mathematica_geosse.m\""
    system(mathematica_wrapper)
    
    #load particular rules saved as a .txt output from Mathematica from a given vector of stationary distribution stored in M matrix
    #
    setwd(out_fp) # switch directory to read mathematica output
    filename_from <- 'output_sim'
    filename_from_upd <- paste(filename_from,'.txt',sep = "") #for a set of stationary distributions in matrix row 'row'
    text <- readLines(filename_from_upd, warn = FALSE)
    # browser()
    #replace all occurrences of "ConditionalExpression" with " " 
    text <- gsub("ConditionalExpression","",text)
    
    #replace all occurences of "Less" with "<","LessEqual" with "<=", "Inequality" with ""
    text <- mgsub(text,c(", Less, ",", LessEqual, ","Inequality"),c("<","<=",""))
    
    #remove comma and whitespace before inequalities symbol
    text <- mgsub(text,c("<=, ","<, ",", <",", <="),c("<=","<","<","<="))
    
    # Assign no solution if the given stationary vector has not solution 
    if (length(text)==0){
      # return("The given stationary frequencies do not have any solutions")
      rules_sets_list <- list(rate_expressions = c(), rules_sets = c(), stationary_A = stat_freq[1], stationary_B = stat_freq[2], stationary_AB = stat_freq[3])
      return(rules_sets_list)
    }
    
    #Extract each analytical expression (not including parameters that are conditional expressions represented by inequalities) for each exact parameter
    current_string <- text #current_string, update it in each loop
    for (i in 1:length(rate_name)){
      # browser()
      count_expr <- 0 # token for the while loop
      j <- 1 # token for the while loop
      while (count_expr != 1){
        rate_expr[i] <- sub(paste("\\,",rate_name[i+j], "->.*"), "", current_string) #get the string until (excluding) the next exact parameter 
        count_expr <- str_count(rate_expr[i],"->")
        j <- j+1
      }
      # This code line above has an issue it skips over some elements from rate_name. If it does, you will have two parameters with -> in the same rate expression
      if (nchar(rate_expr[i]) == nchar(current_string)){
        if (grepl(paste(rate_name[i],"->"),current_string)==TRUE){ #check if this character "rate_name[i] ->" is in the current string
          rate_expr[i] <- current_string #check if the current parameter with the exact expression is the last one in the current string
        } else{
          # rate_expr[i] <- paste("no exact expression for parameter",rate_name[i]) #meaning that that particular parameter is a conditional expression
          rate_expr[i] <- NA #meaning that that particular parameter is a conditional expression (does not have exact expression)
        }
      }
      current_string <- gsub(paste(rate_expr[i],",",sep=""),"",current_string,fixed = TRUE) #update current_string without the string that corresponds to the previous parameter
    }
    
    rate_expr <- rate_expr[!is.na(rate_expr)] #remove all NAs
    for (k in 2:length(rate_expr)){
      # rate_expr[k] <- sub(".*[\\[]([^,]+)[,].*","\\1",rate_expr[k]) #remove conditional expressions from each exact parameters
      rate_expr[k] <- sub(",.*","\\1",rate_expr[k]) #remove conditional expressions from each exact parameters
      rate_expr[k] <- stri_replace_first_fixed(rate_expr[k],"[","") #remove [ at the start of the string of each exact parameter
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
    
    # # convert these strings from rate_expr and rate_cond as mathematical expressions
    # for (l in 1:length(rate_name)){
    #   for (m in 1:length(rate_expr)){
    #     if (rate_name[l] == parse(text=gsub("->.*","\\1",rate_expr[m]))){
    #       # answ_expr <- parse(text = gsub(".*-> ConditionalExpression","\\1",rate_expr[m]))
    #       answ_expr <- parse(text = gsub(".*-> ","\\1",rate_expr[m]))
    #       rate_name[l] <- eval(answ_expr)
    #     }
    #   }
    # }
  }
  # browser()
  # Return rate expressions along with sets of rules between rate parameters 
  rules_sets_list <- list(rate_expressions = rate_expr, rules_sets = rate_cond, stationary_A = stat_freq[1], stationary_B = stat_freq[2], stationary_AB = stat_freq[3])
  #
  # browser()
  # Add "if" condition when the number of alternative rules is only 1
  if (length(rules_sets_list$rules_sets) > 1){
    for (m in 1:(length(rules_sets_list$rules_sets)-1)){
      rules_sets_list$rules_sets[m] <- gsub("\\[|\\]", "",rules_sets_list$rules_sets[m]) # remove all occurences of "[" and "]" from the condition expressions
      rules_sets_list$rules_sets[m] <- substring(rules_sets_list$rules_sets[m], 3, nchar(rules_sets_list$rules_sets[m])-2) # ignore whitespace and "(" at the beginning and ")" and whitespace at the end.
    }
    rules_sets_list$rules_sets[length(rules_sets_list$rules_sets)] <- gsub("\\[|\\]", "",rules_sets_list$rules_sets[length(rules_sets_list$rules_sets)]) # remove all occurences of "[" and "]" from the condition expressions
    rules_sets_list$rules_sets[length(rules_sets_list$rules_sets)] <- substring(rules_sets_list$rules_sets[length(rules_sets_list$rules_sets)], 3, nchar(rules_sets_list$rules_sets[length(rules_sets_list$rules_sets)])-1) # For last conditional expression only: ignore whitespace and "(" at the beginning and ")"  at the end.
  } else{
    rules_sets_list$rules_sets[1] <- gsub("\\[|\\]", "",rules_sets_list$rules_sets[1]) # remove all occurences of "[" and "]" from the condition expressions
  }
  
  # # Save the rule sets for a given stationary frequencies as list 
  # filename <- paste('rules_sets_GeoSSE_',sim_no,'.RData',sep='')
  # save(rules_sets_list,file = filename)
  
  return(rules_sets_list)
}

# get_rules("2-region GeoSSE",c(0.2,0.11,0.69))
