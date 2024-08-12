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

# returns CIs for all x's from 0 to inputted x by default

wald_CI <- function(x, conf.level, all = FALSE) {

  alpha <- 1 - conf.level
  lower <- c(); upper <- c() # initialize necessary vectors
  z_star <- qnorm(1-(alpha/2)) # critical value
  # we are in n = 1 case
  
  if (all == TRUE) { # case to display all CIs from x = 0 to x = obs.
    x <- c(0:x)
    
    # for each x, calculate CI using Wald CI formula
    for (i in 1:length(x)) {
      
      # calculate bounds and store in respective vectors
      lb <- x[i] - z_star*sqrt(x[i])
      ub <- x[i] + z_star*sqrt(x[i])
      
      lower <- c(lower, lb)
      upper <- c(upper, ub)
    }
  }
   else { # case: only display interval for specified observed x
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
# (all == F) returns only CI for the observed (inputted) x, which is default option

rao_score_CI <- function(x, conf.level, all = FALSE) {
  
  alpha <- 1 - conf.level 
  z_star <- qnorm(1-(alpha/2)) # critical value
  
  if (all == TRUE) { # want to display all CIs for x=0 to inputted x
    lower <- c(); upper <- c() # initialize necessary vectors
    x <- c(0:x)
    
    
    for (i in 1:length(x)) {
      # bound calculations
      # using formula from Holladay (2019)
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


wilksLR_CI <- function(x, conf.level = 0.95, all = FALSE, digits = 2) {
  diff <- 1*10^(-digits) # the amount we will increment/decrement lambda by. 
                        # choices are 0.01, 0.001, 0.0001 controlled by user's 
                       # desired decimal place accuracy (input$digits)
  alpha <- (1 - conf.level)
  
  # case 1: only return ONE CI for given obs. x
  if (all == FALSE) {
    t <- x # start at MLE. t represents lambda
    # calculating lower limit
    lower <- c() # initialize lower bound vec
    while (length(lower) == 0) {
      
      if ((2*(log(dpois(x,lambda=x))-log(dpois(x,lambda=t)))) <= qchisq(1-alpha, df=1))  {
        lower <- lower # keeps lower empty
      } else if (x == 0) {
        lower <- 0
      } else {
        lower <- t # sets lower bound to current lambda that breaks inequality
      }
      
      t <- t - diff # decrement lambda until inequality is false. 
                    # Only the lambdas that satisfy the inequality are part of the CI for x
    }
    
    t <- x # reset lambda at MLE (which is observed x)
    
    # start calculating upper limit
    # (same process as lower)
    
    upper <- c() # initialize upper bound vec (same process as lower)
    while (length(upper) == 0) {
      if ((2*(log(dpois(x, lambda = x))-log(dpois(x, lambda = t)))) <= qchisq(1-alpha, df=1)) {
        upper <- upper
      } else {
        upper <- t
        
      }
      
      t <- t + diff
    }
    
    
    
  } else {
    
    # same algorithm as in the above case but User wants to display CIs from x = 0 to x = observed
    
    lower <- c()
    upper <- c()
    
    x <- c(seq(0:x))
    for (i in 1:length(x)) {
      t_l <- x[i] # start lower bound at MLE
      # calculating lower limit
      lb <- c() # initialize lower bound for current x
      while (length(lb) == 0) {
        
        if ((2*(log(dpois(x[i],lambda=x[i]))-log(dpois(x[i],lambda=t_l)))) <= qchisq(1-alpha, df=1))  {
          lb <- lb # keeps lower empty
        # } else if (x[i] == 0) {
        #   lb <- 0
        } else {
          lb <- t_l # sets lower bound to current lambda that breaks inequality
          lower <- c(lower, lb)
        }
        
        t_l <- t_l - diff
        
      }
      
      
      t_u <- x[i] # start upper bound at MLE
      # calculating upper limit
      
      ub <- c() # initialize upper (same process as lower)
      while (length(ub) == 0) {
        if ((2*(log(dpois(x[i], lambda = x[i]))-log(dpois(x[i], lambda = t_u)))) <= qchisq(1-alpha, df=1)) {
          ub <- ub
        } else {
          ub <- t_u # sets upper bound to current lambda that breaks inequality
          upper <- c(upper, ub)
        }
        
        t_u <- t_u + diff
       
      }
      
    }
  }
  
  CIs <- data.frame(x,lower,upper)
  return(CIs)
}
