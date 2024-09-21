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

##-------------------------------------------------------------##
##                     Find Roots Function                     ##
##-------------------------------------------------------------##
# function to return two roots of a given AC: root 1 and root 2

# notes:
# when AC has a = 0, we are not concerned with root 1
# root 1 is lambda when AC first rises above conf.level
# root 2 is lambda when AC first goes below conf.level

find_roots <- function(a, b, conf.level) {
  # getting root 1
  
  
  start <- 0
  end <- AC_max_coords(a,b)$lambda # end root1 search where current AC's max is
  
  # this creates the current curve's function so we can find its root.
  f <- function(lambda) {
    if (a != 0) {
      return((ppois(b, lambda) - ppois(a-1, lambda)) - conf.level)
    } else {
      return(ppois(b, lambda) - conf.level)
    }
  }
  
  if (a != 0) {
    root1 <- uniroot(f, c(start,end))$root
  } else {
    root1 <- NA # we don't need an AC's "left" root if their {a} = 0.
  }
  
  
  ## now finding root 2
  
  # this loop finds the interval of where to search for root2
  start <- AC_max_coords(a,b)$lambda # start at current AC's maximum lambda
  a0 <- a + 1; b0 <- b + 1
  while (sum(dpois(a0:b0, AC_max_coords(a0, b0)$lambda)) >= conf.level) {
    a0 <- a0+1 
    b0 <- b0+1
  } 
  end <- AC_max_coords(a0,b0)$lambda # stop search for root 2 where next AC's maximum occurs
  
  root2 <- uniroot(f, c(start,end))$root  
  
  if(is.na(root1)) {
    root1 <- root2 # this is for the case of ACs with {a} = 0
  }  
  
  
  
  # create df to return
  roots <- data.frame(root1, root2)
  return(roots)
}  




