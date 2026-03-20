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

find_ci <- function(method, x, conf.level, all, digits) {
  obs.x <- x
  
  ci <- c()
  conf.level <- conf.level/100
  digits <- as.numeric(digits) + 1
  if (method == "1") { # Wald
    ci <- Wald.pois(x, conf.level, all)
    method.str <- "Wald (1812)"
    
  } else if (method == "2") { # Rao's Score
    ci <- RaoScore.pois(x, conf.level, all)
    method.str <- "Rao's Score (1948)"
    
  } else if (method == "3") { #Wilks Likelihood Ratio
    ci <- WilksLR.pois(x, conf.level, all, digits)
    method.str <- "Wilks' Likelihood Ratio (1938)"
    
  } else if (method == "4") { # Clopper Pearson/Garwood
    ci <- Garwood.pois(x, conf.level, all)
    method.str <- "Garwood (1936)"
    
  } else if (method == "5") { # MST/OC
    ci <- OC.pois(x, conf.level, all)
    method.str <- "Optimal Coverage (2017)"
    
  } else if (method == "6") { # CG
    ci <- CG.pois(x, conf.level, all)
    method.str <- "Crow & Gardner (1959)"
    
  } else if (method == "7") { # KB
    ci <- KB.pois(x, conf.level, all)
    method.str <- "Kabaila & Byrne (2001)"
    
  } else if (method == "8") { # Blaker's
    ci <- Blaker.pois(x, conf.level)
    method.str <- "Blaker (2000)"
    
  } else if (method == "9") { # CMC
    ci <- CMC.pois(x, conf.level, all)
    method.str <- "Conditional Minimal Cardinality (2023)"
    
  }
    
  
  # adding 2 new columns
  # interval: lower and upper bounds in the form of an open interval
  # input.str: contains the string that the user will see, displaying their inputs
  
  ci <- ci |> 
    mutate(interval = paste0("(",round(lower, digits),", ",round(upper,digits),")"),
           input.str = paste0("method: ", method.str, "\nx-input: ", obs.x, " observed successes",
                              "\nconfidence level: ", 100*conf.level, "%\n"))

  
  
  return(ci) # ci is a data.frame of CIs from x = 0 to x = observed when all = TRUE
             # ci is a data.frame of 1 row for CI of observed x ONLY when all = FALSE
}








