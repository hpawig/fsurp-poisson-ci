##################################################################
##                       Server Functions                       ##
##################################################################


library(tidyverse)

##--------------------------------------------------------------##
##                          function 1                          ##
##--------------------------------------------------------------##

# function that takes user's inputted method 
# and outputs the corresponding method's CI
# for the user's observed x

find_ci <- function(method, x, conf.level, all) {
  ci <- c()
  conf.level <- conf.level/100
  
  if (method == "1") { # Wald
    ci <- wald_CI(x, conf.level, all)
  } else if (method == "2") { # Rao's Score
    ci <- rao_score_CI(x, conf.level, all)
  } else if (method == "3") { #Wilks Likelihood Ratio
    ci <- wilksLR_CI(x, conf.level, all)
  } else if (method == "4") { # Clopper Pearson
    ci <- clopperPearson_CI(x, conf.level, all)
  } else if (method == "5") { # MST/OC
    ci <- OC(x, conf.level, all)
  } else if (method == "6") { # CG
    ci <- CG(x, conf.level, all)
  }
  

  return(ci) # ci is a data.frame of CIs from x = 0 to x = observed when all = TRUE
             # ci is a data.cframe of 1 row for CI of observed x ONLY when all = FALSE
}











