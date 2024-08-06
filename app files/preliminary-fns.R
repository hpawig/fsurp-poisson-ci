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


##--------------------------------------------------------------##
##                    Test Coverage function                    ##
##  returns True if AC(a-b) ever rises above confidence level.  ##
##                        otherwise False                       ##
##--------------------------------------------------------------##

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
##                Min-Tail Probability Function                ##
##-------------------------------------------------------------##



# calculate MTP for Poisson rvs

get_mtp <- function(x, lambda) {
  if (x != 0) {
    return(# general case
      min(ppois(x, lambda), # P(X<=x)
          (1-ppois(x-1, lambda)))  # P(X >= x)
    ) 
  } else {
    return(
      min(ppois(0, lambda), 1) # P(X >= 0) = 1 always
    ) 
  }
}

