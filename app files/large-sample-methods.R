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

wald_CI <- function(K, conf.level, all = FALSE) {
  x <- c(0:K)
  conf.level <- conf.level/100 # conf.level must be as a %
  alpha <- 1 - conf.level
  lower <- c(); upper <- c() # initialize necessary vectors
  z_star <- qnorm(1-(alpha/2)) # critical value
  # we are in n = 1 case
  
  for (i in 1:length(x)) {
    
    # calculate bounds and store in respective vectors
    lb <- x[i] - z_star*sqrt(x[i])
    ub <- x[i] + z_star*sqrt(x[i])

    lower <- c(lower, lb)
    upper <- c(upper, ub)
  }

  
  CIs <- data.frame(x, lower, upper)
  print(CIs)
  if (all == FALSE) {
    CIs <- CIs |> 
      filter(x == K)
  }
  
  return(CIs)
}


##--------------------------------------------------------------##
##                      Rao's Score Method                      ##
##--------------------------------------------------------------##

# returns CIs for x = 0 to x = K. 
# (all == F) returns only CI for x = K (the observed x)

rao_score_CI <- function(K, conf.level, all = FALSE) {
  x <- c(0:K)
  alpha <- 1-(conf.level/100) # conf.level is a whole number
  lower <- c(); upper <- c() # initialize necessary vectors
  
  # calculate critical value
  z_star <- qnorm(1-(alpha/2)) # critical value
  
  for (i in 1:length(x)) {
    # bound calculations
    lb <- (x[i] + (1/2)*z_star^2) - z_star*sqrt(x[i] + (1/4)*z_star^2)
    ub <- (x[i] + (1/2)*z_star^2) + z_star*sqrt(x[i] + (1/4)*z_star^2)
    
    lower <- c(lower,lb); upper <- c(upper,ub)
  }
  
  CIs <- data.frame(x,lower,upper)
  if (all == FALSE) {
    CIs <- CIs |>
      filter(x == K)
  }
  
  return(CIs)
}



##--------------------------------------------------------------##
##                Wilks' Likelihood Ratio Method                ##
##--------------------------------------------------------------##

# returns CIs for only x == K by default

wilksLR_CI <- function(K, conf.level, all = FALSE) {
  lower <- c()
  upper <- c()

  alpha <- 1-(conf.level/100) # conf.level is a %
  x <- c(0:K)
  
  for (i in 1:length(x)) {
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

