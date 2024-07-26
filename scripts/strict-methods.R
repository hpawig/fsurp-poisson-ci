#################################################################
##                 Strict  Methods Script                      ##
#################################################################


##--------------------------------------------------------------##
##                        loading packages                      ##
##--------------------------------------------------------------##

library(tidyverse)

##-------------------------------------------------------------##
##                 Clopper-Pearson for Poisson                 ##
##-------------------------------------------------------------##


# returns only 1 confidence interval for inputted x

clopperPearson_CI <- function(x, conf.level = 0.95) {
  x <- sum(x)
  alpha <- 1-conf.level
  
  # lower bound
  lwr <- (1/2)*qchisq(p = alpha/2, df = 2*x)
  
  # upper bound
  uppr <- (1/2)*qchisq(p = 1-alpha/2, df = 2*(x+1))
  
  return(data.frame(lwr, uppr))
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






















