#################################################################
##                 Strict  Methods Script                      ##
#################################################################
##                       For Poisson                           ##
#################################################################


# includes Minimal Cardinality Procedures
# utilizes functions in file "preliminary-fns.R"
source("preliminary-fns.R", encoding = "UTF-8")

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
  alpha <- 1-conf.level
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





##--------------------------------------------------------------##
##                        Crow & Gardner                        ##
##--------------------------------------------------------------##

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
    
    if((test_coverage(a+1, b+1, conf.level) == T)) { # check AC {a+1}-{b+1} first
      
      a <- a + 1
      b <- b + 1   
      
      
      # setting coincidental endpoints
      # current b's lower bound where (a+1)-(b+1) rises above conf.level
      lower[b+1] <- find_roots(a,b,conf.level)$root1
      
      # a-1's upper bound is also where (a+1)-(b+1) curve comes above CI
      upper[a] <- lower[b+1]     
      
      
      # check next AC by increasing cardinality by 1 but also ensure {a} non-decreasing       
    } else if ((test_coverage(a+1, b+2, conf.level) == F)) { 
      
      lower[(b+1)+1] <- find_roots(a,b,conf.level)$root2 # lower limit for (a+1), aka new b
      b <- b + 1   
      
    }  else if (test_coverage(a+1, b+2, conf.level) == T) {
      
      # setting coincidental endpoint
      # identical lower endpoint when AC {a+1}--{b+2} is above conf.level
      lower[(b+2)+1] <- find_roots(a,b,conf.level)$root2 # b+2's lower bound
      lower[(b+1)+1] <- lower[b+3]  # b+1's lower bound
      upper[a+1] <- find_roots(a,b,conf.level)$root2  # a's upper bound
      
      a <- a + 1
      b <- b + 2
      
    } 
  }
  
  lower <- lower[1:(K+1)]
  upper <- upper[1:(K+1)]
  x <- x[1:(K+1)] # ensuring all vectors are equal in length. 
  # b/c extra lower endpoints will be cut off
  
  
  CIs <- data.frame(x,lower,upper)
  
  # all=TRUE represents display all intervals from x=0 to x=observed/user input
  # use filter to keep rows of the data set "CI" where x = observed (K) and discards those where x != K
  # to get only 1 row for observed x's CI
  if (all == FALSE) {
    CIs <- CIs |> 
      filter(x == K)
  }
  return(CIs)
}




##--------------------------------------------------------------##
##                    Blaker's Method (2000)                    ##
##--------------------------------------------------------------##

# Process
# 1) find min tail probability (MTP) of observed x
# 2) Fix lambda; for all x, find MTPs that are as small or smaller than observed MTP at that lambda. The observed x is always among these
# 3) observed x is in acceptance set of lambda if P(x in A_lambda) > alpha. Record all lambdas
# 4) lower and upper confidence limits given obs.x are smallest and largest lambdas
#   that have acceptance sets with obs. x , respectively

blaker_CI <- function(K, conf.level, all = F){
  obs.x <- K
  if (K==0) {
    K <- 1
    obs.x <- 0
  }
  
  #Acceptance function
  #Computes probabiility of observing something with tail probability as small as x
  accept.blaker.pois <- function(x, lambda) {
    p1 <- 1
    if(x != 0) {
      p1 <- 1-ppois(x-1, lambda = lambda) # right tail probability of x  
    }
    p2 <- ppois(x, lambda = lambda) # left tail probability of x
    
    a1 <- p1 + ppois((qpois(p1, lambda)-1), lambda = lambda) 
    a2 <- p2 + (1-ppois(qpois(1-p2, lambda), lambda = lambda))
    return(min(a1,a2))
  }
  

  tol <- 0.0001 #decimal accuracy of ci
  LL <- NA; UL <- NA
  LL[1]=0 # lower for x = 0
  
  u <- 0
  # lambda included in ci for x if above prob (accept.blaker) > alpha  
  # First determine mu's in ci with grid of 10^3*tol, then and then narrow down to 
  # 10^2*tol, 10*tol and finally to tol decimal place accuracy
  while(accept.blaker.pois(x=0, lambda=u) >= (1-conf.level)){u=u+10^3*tol}; u=u-10^3*tol
  while(accept.blaker.pois(x=0, lambda=u) >= (1-conf.level)){u=u+10^2*tol}; u=u-10^2*tol 
  while(accept.blaker.pois(x=0, lambda=u) >= (1-conf.level)){u=u+10*tol}; u=u-10*tol 
  while(accept.blaker.pois(x=0, lambda=u) >= (1-conf.level)){u=u+tol}

  UL[1]=u
  
  for(x in 1:K){
    
    l=x
    u=x
    
    # lambda included in ci for x if above prob (accept.blaker) > alpha  
    # First determine mu's in ci with grid of 10^3*tol, then and then narrow down to 
    # 10^2*tol, 10*tol and finally to tol decimal place accuracy
    while((accept.blaker.pois(x=x,lambda=l) >= (1-conf.level)) && l!=0 ){l=l-10^3*tol}; l=l+10^3*tol
    while((accept.blaker.pois(x=x, lambda=l) >= (1-conf.level)) && l!=0 ){l=l-10^2*tol}; l=l+10^2*tol
    while((accept.blaker.pois(x=x, lambda=l) >= (1-conf.level)) && l!=0 ){l=l-10*tol}; l=l+10*tol
    while((accept.blaker.pois(x=x,lambda=l) >= (1-conf.level)) && l!=0 ){l=l-tol}
    
    #initially use larger grid e.g. 10^5*tol for upper endpoints comparison to lower endpoint
    #use especially large grid 10^7*tol and 10^6*tol for upper endpoints when k=1,2
    # if(k<=2){
    #   while(accept.blaker.pois(x=x, lambda=u) >= (1-conf.level)){u=u+10^7*tol}; u=u-10^7*tol
    #   while(accept.blaker.pois(x=x, lambda=u) >= (1-conf.level)){u=u+10^6*tol}; u=u-10^6*tol
    # }
    while(accept.blaker.pois(x=x, lambda=u) >= (1-conf.level)){u=u+10^5*tol}; u=u-10^5*tol
    while(accept.blaker.pois(x=x, lambda=u) >= (1-conf.level)){u=u+10^4*tol}; u=u-10^4*tol
    while(accept.blaker.pois(x=x, lambda=u) >= (1-conf.level)){u=u+10^3*tol}; u=u-10^3*tol
    while(accept.blaker.pois(x=x, lambda=u) >= (1-conf.level)){u=u+10^2*tol}; u=u-10^2*tol 
    while(accept.blaker.pois(x=x, lambda=u) >= (1-conf.level)){u=u+10*tol}; u=u-10*tol 
    while(accept.blaker.pois(x=x, lambda=u) >= (1-conf.level)){u=u+tol}
    
    LL[x+1]=l
    UL[x+1]=u-tol/100000 #subtract tol/100000 so that we have half open intervals [l,u)
    #print(x)
    
  } 
  
  
  CIs <- data.frame(x=0:K, lower=LL, upper=UL)  
  
  if (all == F) {
    CIs <- CIs |> 
      filter(x == obs.x)
  }
  
  return(CIs)
}





##-------------------------------------------------------------##
##             Conditional Minimal Cardinality (CMC)           ##
##-------------------------------------------------------------##



# Uses functions from 'preliminary-fns.R' in app files
# Based off of the CMC for Neg Binom by Dr. B.A. Holladay


# K: largest value of x of interest
# CMC method function for the Poisson distribution

CMC <-function(K, conf.level, all = F) {
  if (K <= 5) {
    obs.x <- K # temporarily store original K given by user into object obs.x
    K <- 20 # for some reason, algorithm won't output right CIs when x = 0,...,5
    # so temporarily will compute many intervals (up to x = 20) and then restore original K 
  } else {
    obs.x <- K
  }
  
  # Determining m(a) for all a <= K+1
  # For fixed a, m(a) is the smallest b such that P(a<=X<=b)>=conf.level
  # so that AC(a,m(a)) is the core of rainbow, RB(a)
  # m(a)'s are needed b/c root1(a+1,m(a+1)) determines u(a)
  a <- 0; b <- 0
  m <- c() # vector that will hold all m(a)'s to keep track of each RB core
  
  while(a <= K+1){
    
    while(test_coverage(a, b, conf.level) == F) { # while AC's coverage is BELOW conf.level
      # transition from curves in the same rainbow when AC(a,b) < conf.level
      b <- b+1
    } # once this loop ends, we have the desired b = m(a)
    
    m[a+1] <- b
    a <- a+1
    b <- b+1 #can start the search at m(a)+1=b+1 b/c m(a+1)>=m(a)+1
  }
  
  
  x <- 0:K # initialize vector of x values from 0 to K
  lower <- rep(NA, K+1) # initialize empty vectors for lower and upper bounds
  upper <- rep(NA, K+1)
  
  a <- 0; b <- 0
  lower[a+1] <- 0 # (lower bound for x=0)
  
  
  # run until lower(K) determined
  while(is.na(lower[K+1])){
    
    # staying on RB(a) until RB(a+1) core aka AC[(a+1)-m(a+1)] is above conf.level
    # so in this loop we increment b staying on RB(a) until next core rises above conf.level
    
    while((find_roots(a, b, conf.level)$root2 < find_roots(a+1, m[a+2], conf.level)$root1)) {
      
      # in this loop, we're transitioning between ACs within a RB(a)
      # root2 of AC[a-b] is the lower bound for b+1 
      lower[(b+1)+1] <- find_roots(a, b, conf.level)$root2
      b <- b + 1
      
      
    }
    
    
    
    
    
    # when transitioning between rainbows RB(a) to RB(a+1) we move from curve AC(a,b) to core of 
    # next rainbow AC(a+1, m(a+1)). The location of this transition occurs at root1(a+1,m(a+1)) 
    # and determines both the upper endpoint for a u(a) and the lower endpoints for b+1,...,m(a+1),
    # l(b+1)=...=l(m(a+1)).
    
    # now we set the coincidental endpoint that happens when transitioning RBs
    b.temp <- b
    a <- a+1
    b <- m[a+1]
    
    for(i in ((b.temp+1):min(b:K))) {
      lower[i+1] <- find_roots(a, m[a+1], conf.level)$root1
    }
    upper[(a-1)+1] <- find_roots(a, m[a+1],conf.level)$root1
  }
  
  # determine upper endpoints for remaining x:get upper(x) for a <= x <= K
  # once we have lower endpoints for x up to K we can work on upper endpoints 
  # separately. 
  # These remaining values of upper(x) are determined by upper(x)=get_roots(x+1,m(x+1))$root1 
  for(x in a:K){
    upper[x+1] <- find_roots(x,m[x+1],conf.level)$root1
    
  }
  
  K <- obs.x # restore original K...
  
  lower <- lower[1:(K+1)]
  upper <- upper[1:(K+1)]
  
  CIs <- data.frame(x = 0:K, lower, upper)
  
  if (all == F) { # indicates that user only wants to output interval for x = K
    CIs <- CIs |> 
      filter(x == K) # "filter" only keeps the row in df "CIs" where x = K
  }
  return(CIs)
}

