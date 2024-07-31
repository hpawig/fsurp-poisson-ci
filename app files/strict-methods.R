#################################################################
##                 Strict  Methods Script                      ##
#################################################################


##--------------------------------------------------------------##
##             loading packages & R scripts                     ##
##--------------------------------------------------------------##

library(tidyverse)


##-------------------------------------------------------------##
##                 Clopper-Pearson for Poisson                 ##
##-------------------------------------------------------------##


# K = observed x...
# conf.level (%)

clopperPearson_CI <- function(K, conf.level, all = FALSE) {
  x <- c(0:K)
  alpha <- 1-(conf.level/100)
  lower <- c(); upper <- c()
  
  
  for (i in 1:length(x)) { # cycles through all x = 0, 1, ..., K
    # calculate current x's lower bound
    lb <- (1/2)*qchisq(p = alpha/2, df = 2*x[i])
    lower <- c(lower, lb)
    
    # calculate current upper bound
    ub <- (1/2)*qchisq(p = 1-alpha/2, df = 2*(x[i]+1))
    upper <- c(upper, ub) 
  }
  
  CIs <- data.frame(x,lower, upper)
  if (all == F) {
    CIs <- CIs |> 
      filter(x==K) # filters rows where x = observed (K). [only returns CI for observed x]
  }
  return(CIs)
}


##--------------------------------------------------------------##
##                  Acceptance Curve Functions                  ##
##--------------------------------------------------------------##


##----------------------------------------------------------------------------------------
##  AC Max Coordinates function finds the (x=lambda,y=coverage probability) coordinate  --
##                              where the AC's max occurs                               --
##----------------------------------------------------------------------------------------

AC_max_coords <- function(a,b) {
  # using prop 2.2a formula. 
  lambda <- prod(a:b)^(1/(b-a+1)) #prop 2.2a
  
  if (a == 0) {
    max_prob <- ppois(b, lambda) # special case of calculating cdf when a is 0
  } else {
    max_prob <- ppois(b, lambda) - ppois(a-1, lambda) # general formula
  }
  
  return(data.frame(lambda, max_prob))
}


##----------------------------------------------------------------
##                    Test Coverage function                    --
##  returns True if AC(a-b) ever rises above confidence level.  --
##                        otherwise False                       --
##----------------------------------------------------------------

test_coverage <- function(a,b, conf.level) {
  
  # this is maximum prob for current cardinality curve
  
  max_prob <- AC_max_coords(a, b)$max_prob
  
  # cases to decide how to change a-b: look at if max is above conf level
  
  if (max_prob > conf.level) {
    
    return(TRUE)
    
  } else {
    
    return(FALSE)
  }
}





##-------------------------------------------------------------##
##               Modified Stern/Optimal Coverage               ##
##-------------------------------------------------------------##


# this function returns a data frame of all CIs for x = 0 to x = K
# K = "largest number of x to create a CI for"
# all = TRUE: means create CIs for 0 to K.

OC <- function(K, conf.level, all = FALSE) {
  x <- c(0:K)

  
  # initialize necessary vectors
  lower <- c()
  upper <- c()
  
  # step 1  
  a <- 0 # starting a
  b <- 0 # starting b
  lower[1] <- 0
  
  
  
  while (a < (K+1)) {
    # step 2
    
    # this creates the current curve's function so we can find its root.    
    f <- function(lambda) { 
      if (b != 0) {
        return((ppois(b, lambda) - ppois(a-1, lambda)) - conf.level)
      } else {
        return(ppois(b, lambda) - conf.level)
      }
    }  
    
    
    if((test_coverage(a+1, b+1, conf.level) != T)) {
      
      
      # this loop finds interval (start,end) to search for current a-b's root
      start <- AC_max_coords(a, b)$lambda # start searching at current a-b's maximum lambda
      a0 <- (a+1); b0 <- (b+1)
      while (sum(dpois(a0:b0,
                       AC_max_coords(a0, b0)$lambda)) >= conf.level) {
        a0 <- a0+1 
        b0 <- b0+1
      } 
      end <- AC_max_coords(a0, b0)$lambda # end search at next a-b's maximum lambda
      
      
      # there's an "extra" +1 for indexing purposes
      lower[(b+1)+1] <- uniroot(f, c(start, end))$root # lower limit for (b+1), aka new b
      b <- b + 1
      
      
    } else {
      # step 3
      
      a <- a + 1
      b <- b + 1
      
      # setting coincidental endpoint (by prop 2.2c)
      lower[b+1] <- AC_max_coords(a, b)$lambda # current b's lower bound
      upper[a] <- lower[b+1]              # a-1's lower bound
      
    }
    
  }
  lower <- lower[1:(K+1)]
  upper <- upper[1:(K+1)]
  x <- x[1:(K+1)] # ensuring all vectors are equal in length. 
  # b/c extra lower endpoints will be cut off


  
  CIs <- data.frame(x, lower, upper)
  if (all == FALSE) {
    CIs <- CIs |> 
      filter(x == K)
  }
  
  return(CIs)
}





##----------------------------------------------------------------
##                        Crow & Gardner                         -
##----------------------------------------------------------------


# K = largest number of x to create a CI for.
# K = largest number of x to create a CI for.
CG <- function(K, conf.level, all = FALSE) {
  x <- c(0:K)
  
  # initialize necessary vectors
  lower <- c()
  upper <- c()
  
  # step 1  
  a <- 0 # starting a
  b <- 0 # starting b
  lower[1] <- 0
  
  
  
  while (a < (K+1)) {
    
    
    # this creates the current curve's function so we can find its root.    
    f <- function(lambda) { 
      if (b != 0) {
        return((ppois(b, lambda) - ppois(a-1, lambda)) - conf.level)
      } else {
        return(ppois(b, lambda) - conf.level)
      }
    }  
    
    
    
    # this loop finds interval (start,end) to search for current a-b's root
    start <- AC_max_coords(a, b)$lambda # start searching at current a-b's maximum lambda
    a0 <- (a+1); b0 <- (b+1)
    while (sum(dpois(a0:b0,
                     AC_max_coords(a0, b0)$lambda)) >= conf.level) {
      a0 <- a0+1 
      b0 <- b0+1
    } 
    end <- AC_max_coords(a0, b0)$lambda # end search at next a-b's maximum lambda
    
    
    if((test_coverage(a+1, b+2, conf.level) != T)) {
      
      
      # there's an "extra" +1 for indexing purposes
      lower[(b+2)] <- uniroot(f, c(start, end))$root # lower limit for (b+1), aka new b
      b <- b + 1
      
      
    } else if ((test_coverage(a+1, b+1, conf.level) == T)) {
      
      a <- a + 1
      b <- b + 1
      
      start <- AC_max_coords(a, b)$lambda # start searching at current a-b's maximum lambda
      a0 <- (a+1); b0 <- (b+1)
      while (sum(dpois(a0:b0,
                       AC_max_coords(a0, b0)$lambda)) >= conf.level) {
        a0 <- a0+1
        b0 <- b0+1
      }
      end <- AC_max_coords(a0, b0)$lambda # end search at next a-b's maximum lambda
      
      
      # setting coincidental endpoint
      # which will be when the previous AC hits conf.level
      lower[b+1] <- uniroot(f, c(start, end))$root  # current b's lower bound
      
      
      start <- AC_max_coords(a,b)$lambda            #  start search at (a-1), (b-1) max.
      end <- AC_max_coords(a-1, b-1)$lambda         # end search at a+1, b+1 maximum.
      upper[a] <- uniroot(f, c(start,end))$root     # a-1's upper bound. where (a+1)-(b+1) curve starts above CI
      
      
    }  else if (test_coverage(a+1, b+2, conf.level) == T) {
      start <- AC_max_coords(a, b)$lambda
      a0 <- a + 1; b0 <- b + 1
      end <- AC_max_coords(a0, b0)$lambda # end search at next a-b's maximum lambda
      
      # setting coincidental endpoint
      
      lower[b+3] <- uniroot(f, c(start, end))$root  # new b's lower bound
      lower[b+2] <- lower[b+3]  # identical lower endpoint when a + 1 and b + 2 is above conf.level
      upper[a+1] <- lower[b+3]                        # old a's upper bound
      
      
      a <- a + 1
      b <- b + 2
      
      
    } 
  }
  
  lower <- lower[1:(K+1)]
  upper <- upper[1:(K+1)]
  x <- x[1:(K+1)] # ensuring all vectors are equal in length. 
  # b/c extra lower endpoints will be cut off
  
  
  CIs <- data.frame(x,lower,upper)
  
  if (all == FALSE) {
    CIs <- CIs |> 
      filter(x == K)
  }
  return(CIs)
}








