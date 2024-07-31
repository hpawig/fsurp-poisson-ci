#################################################################
##                 Large Sample Methods Script                 ##
#################################################################



##--------------------------------------------------------------##
##                        loading packages                      ##
##--------------------------------------------------------------##

library(tidyverse)



##-------------------------------------------------------------##
##                          Wald Method                        ##
##-------------------------------------------------------------##

# returns CIs for all x's from 0 to inputted x ("K") by default

wald_CI <- function(x, conf.level, all = FALSE) {

  alpha <- 1 - conf.level
  lower <- c(); upper <- c() # initialize necessary vectors
  z_star <- qnorm(1-(alpha/2)) # critical value
  # we are in n = 1 case
  
  if (all == TRUE) {
    x <- c(0:x)
    for (i in 1:length(x)) {
      
      # calculate bounds and store in respective vectors
      lb <- x[i] - z_star*sqrt(x[i])
      ub <- x[i] + z_star*sqrt(x[i])
      
      lower <- c(lower, lb)
      upper <- c(upper, ub)
    }
  }
   else {
     lower <- x - z_star*sqrt(x)
     upper <- x + z_star*sqrt(x)
   }

  
  CIs <- data.frame(x, lower, upper)

  
  return(CIs)
}


##--------------------------------------------------------------##
##                      Rao's Score Method                      ##
##--------------------------------------------------------------##

# returns CIs for x = obs
# (all == F) returns only CI for the observed (inputted) x

rao_score_CI <- function(x, conf.level, all = FALSE) {
  
  alpha <- 1 - conf.level 
  z_star <- qnorm(1-(alpha/2)) # critical value
  
  if (all == TRUE) { # want to display all CIs for x=0 to inputted x
    lower <- c(); upper <- c() # initialize necessary vectors
    x <- c(0:x)
    
    
    for (i in 1:length(x)) {
      # bound calculations
      lb <- (x[i] + (1/2)*z_star^2) - z_star*sqrt(x[i] + (1/4)*z_star^2)
      ub <- (x[i] + (1/2)*z_star^2) + z_star*sqrt(x[i] + (1/4)*z_star^2)
      
      lower <- c(lower,lb); upper <- c(upper,ub)
    }
  } else {
    lower <- (x + (1/2)*z_star^2) - z_star*sqrt(x + (1/4)*z_star^2)
    upper <- (x + (1/2)*z_star^2) + z_star*sqrt(x + (1/4)*z_star^2)
  }

  CIs <- data.frame(x,lower,upper)
  
  return(CIs)
}



##--------------------------------------------------------------##
##                Wilks' Likelihood Ratio Method                ##
##--------------------------------------------------------------##

# returns CIs for observed x by default

wilksLR_CI <- function(K, conf.level, all = FALSE) {
  alpha <- 1 - conf.level
  
  lower <- c()
  upper <- c()

  x <- c(0:K)
  
  for (i in 1:length(x)) {
    if (all == FALSE) {
      i <- length(x)
    }
    
    
    t <- K # start at MLE. t represents lambda
    
    # calculating lower limit
    lb <- c() # initialize lower bound
    
    while (length(lb) == 0) {
      
      if (x[i] == 0)  {
        lb <- 0 
        lower <- c(lower, lb)
      } else if ((2*(log(dpois(x[i],lambda=x[i]))-log(dpois(x[i],lambda=t)))) <= qchisq(1-alpha, df=1)){
        lb <- lb # keeps lower bound empty
      } else {
        lb <- t # sets lower bound to current lambda that breaks inequality
        lower <- c(lower, lb)
      }
      
      t <- t - 0.0001
    }
    
    t <- x[i] # reset lambda at MLE
    
    # calculating upper limit
    
    ub <- c() # initialize upper (same process as lower)
    
    while (length(ub) == 0) {
      if ((2*(log(dpois(x[i],lambda=x[i]))-log(dpois(x[i],lambda=t)))) <= qchisq(1-alpha, df=1)) {
        ub <- ub
      } else {
        ub <- t
        upper <- c(upper, ub)
      }
      
      t <- t + 0.0001

    }
  }
 
  CIs <- data.frame(x, lower, upper)
  if (all == FALSE) {
    CIs <- CIs |>
      filter(x == K)
  }
  
  return(CIs)
}

