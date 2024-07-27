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

# returns CI only for inputted x

wald_CI <- function(x, conf.level = 0.95) {
  alpha <- 1 - conf.level
  # we are in n = 1 case
  lambda_hat <- x
  z_star <- qnorm(1-(alpha/2)) # critical value
  
  # lower bound 
  lwr <- lambda_hat - z_star*sqrt(x)
  uppr <- lambda_hat + z_star*sqrt(x)
  
  ci <- data.frame(lwr,uppr)
  return(ci)
}


##--------------------------------------------------------------##
##                      Rao's Score Method                      ##
##--------------------------------------------------------------##

# returns CI only for inputted x

rao_score_CI <- function(x, conf.level = 0.95) {
  alpha <- 1-conf.level
  # we are in a single observation case
  lambda_hat <- mean(x) 
  # calculate critical value
  z_star <- qnorm(1-(alpha/2)) # critical value
  
  
  # bound calculations
  lwr <- (lambda_hat + (1/2)*z_star^2) - z_star*sqrt(lambda_hat + (1/4)*z_star^2)
  uppr <- (lambda_hat + (1/2)*z_star^2) + z_star*sqrt(lambda_hat + (1/4)*z_star^2)
  
  return(data.frame(lwr,uppr))
}



##--------------------------------------------------------------##
##                Wilks' Likelihood Ratio Method                ##
##--------------------------------------------------------------##

# returns CI only for inputted x

wilksLR_CI <- function(x, conf.level = 0.95) {
  x <- sum(x) # formula is for a sum 
  n <- length(x)
  t <- x # start at MLE. t represents lambda
  alpha <- 1-conf.level
  
  
  # calculating lower limit
  lwr <- c() # initialize lower bound
  while (length(lwr) == 0) {
    
    if ((2*(log(dpois(x,lambda=x))-log(dpois(x,lambda=t)))) <= qchisq(1-alpha, df=1))  {
      lwr <- lwr # keeps lwr empty
    } else {
      lwr <- t # sets lwr to current lambda that breaks inequality
    }
    
    t <- t - 0.001
  }
  
  t <- x # reset lambda at MLE
  
  # calculating upper limit
  
  uppr <- c() # initialize uppr (same process as lower)
  while (length(uppr) == 0) {
    if ((2*(log(dpois(x, lambda = x))-log(dpois(x, lambda = t)))) <= qchisq(1-alpha, df=1)) {
      uppr <- uppr
    } else {
      uppr <- t
      
    }
    
    t <- t + 0.001
  }
  return(data.frame(lwr,uppr))
}

