#################################################################
##                 Strict  Methods Script                      ##
#################################################################
##                       For Poisson                           ##
#################################################################

##--------------------------------------------------------------##
##             loading packages & R scripts                     ##
##--------------------------------------------------------------##

library(tidyverse)

# Each Strict Method is the confidence procedure for a single Poisson


# includes Minimal Cardinality Procedures
# utilizes functions in file "preliminary-fns.R"
source("preliminary-fns.R", encoding = "UTF-8")


##-------------------------------------------------------------##
##                     Garwood for Poisson                     ##
##-------------------------------------------------------------##

# AKA Clopper Pearson

# K = observed x...
# conf_level (%)
# by default, function returns one (1) confidence interval for x = given K

Garwood.pois <- function(K, conf_level, all = FALSE) {
  x <- c(0:K)
  alpha <- 1-conf_level
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
##               Modified Sterne/Optimal Coverage               ##
##-------------------------------------------------------------##

# Minimal Cardinality Procedure

# this function returns a data frame of all CIs for x = 0 to x = K
# K = "largest number of x to create a CI for"
# all = TRUE: means create CIs for 0 to K.

OC.pois <- function(K, conf_level, all = FALSE) {
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
        return((ppois(b, lambda) - ppois(a-1, lambda)) - conf_level)
      } else {
        return(ppois(b, lambda) - conf_level)
      }
    }  
    
    
    if((test_coverage(a+1, b+1, conf_level) != T)) {
      
      
      # this loop finds interval (start,end) to search for current a-b's root
      start <- AC_max_coords(a, b)$lambda # start searching at current a-b's maximum lambda
      a0 <- (a+1); b0 <- (b+1)
      while (sum(dpois(a0:b0,
                       AC_max_coords(a0, b0)$lambda)) >= conf_level) {
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

# Minimal Cardinality Procedure

# K = largest number of x to create a CI for.
CG.pois <- function(K, conf_level, all = FALSE) {
  x <- c(0:K)
  
  # initialize necessary vectors
  lower <- c()
  upper <- c()
  
  # step 1  
  a <- 0 # starting a
  b <- 0 # starting b
  lower[1] <- 0
  
  
  
  while (a < (K+1)) {
    
    if((test_coverage(a+1, b+1, conf_level) == T)) { # check AC {a+1}-{b+1} first
      
      a <- a + 1
      b <- b + 1   
      
      
      # setting coincidental endpoints
      # current b's lower bound where (a+1)-(b+1) rises above conf_level
      lower[b+1] <- find_roots(a,b,conf_level,root=1)
      
      # a-1's upper bound is also where (a+1)-(b+1) curve comes above CI
      upper[a] <- lower[b+1]     
      
      
      # check next AC by increasing cardinality by 1 but also ensure {a} non-decreasing       
    } else if ((test_coverage(a+1, b+2, conf_level) == F)) { 
      
      lower[(b+1)+1] <- find_roots(a,b,conf_level,root=2) # lower limit for (a+1), aka new b
      b <- b + 1   
      
    }  else if (test_coverage(a+1, b+2, conf_level) == T) {
      
      # setting coincidental endpoint
      # identical lower endpoint when AC {a+1}--{b+2} is above conf_level
      lower[(b+2)+1] <- find_roots(a,b,conf_level,root=2) # b+2's lower bound
      lower[(b+1)+1] <- lower[b+3]  # b+1's lower bound
      upper[a+1] <- find_roots(a,b,conf_level,root=2)  # a's upper bound
      
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

Blaker.pois <- function(K, conf_level, all = F){
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
  while(accept.blaker.pois(x=0, lambda=u) >= (1-conf_level)){u=u+10^3*tol}; u=u-10^3*tol
  while(accept.blaker.pois(x=0, lambda=u) >= (1-conf_level)){u=u+10^2*tol}; u=u-10^2*tol 
  while(accept.blaker.pois(x=0, lambda=u) >= (1-conf_level)){u=u+10*tol}; u=u-10*tol 
  while(accept.blaker.pois(x=0, lambda=u) >= (1-conf_level)){u=u+tol}
  
  UL[1]=u
  
  for(x in 1:K){
    
    l=x
    u=x
    
    # lambda included in ci for x if above prob (accept.blaker) > alpha  
    # First determine mu's in ci with grid of 10^3*tol, then and then narrow down to 
    # 10^2*tol, 10*tol and finally to tol decimal place accuracy
    while((accept.blaker.pois(x=x,lambda=l) >= (1-conf_level)) && l!=0 ){l=l-10^3*tol}; l=l+10^3*tol
    while((accept.blaker.pois(x=x, lambda=l) >= (1-conf_level)) && l!=0 ){l=l-10^2*tol}; l=l+10^2*tol
    while((accept.blaker.pois(x=x, lambda=l) >= (1-conf_level)) && l!=0 ){l=l-10*tol}; l=l+10*tol
    while((accept.blaker.pois(x=x,lambda=l) >= (1-conf_level)) && l!=0 ){l=l-tol}
    
    #initially use larger grid e.g. 10^5*tol for upper endpoints comparison to lower endpoint
    #use especially large grid 10^7*tol and 10^6*tol for upper endpoints when k=1,2
    # if(k<=2){
    #   while(accept.blaker.pois(x=x, lambda=u) >= (1-conf_level)){u=u+10^7*tol}; u=u-10^7*tol
    #   while(accept.blaker.pois(x=x, lambda=u) >= (1-conf_level)){u=u+10^6*tol}; u=u-10^6*tol
    # }
    while(accept.blaker.pois(x=x, lambda=u) >= (1-conf_level)){u=u+10^5*tol}; u=u-10^5*tol
    while(accept.blaker.pois(x=x, lambda=u) >= (1-conf_level)){u=u+10^4*tol}; u=u-10^4*tol
    while(accept.blaker.pois(x=x, lambda=u) >= (1-conf_level)){u=u+10^3*tol}; u=u-10^3*tol
    while(accept.blaker.pois(x=x, lambda=u) >= (1-conf_level)){u=u+10^2*tol}; u=u-10^2*tol 
    while(accept.blaker.pois(x=x, lambda=u) >= (1-conf_level)){u=u+10*tol}; u=u-10*tol 
    while(accept.blaker.pois(x=x, lambda=u) >= (1-conf_level)){u=u+tol}
    
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



# utilizes functions in file "preliminary-fns.R"
# K: largest value of x of interest
# Based on NB CMC code by Doi, Schilling and Holladay (2023)

CMC.pois <-function(K, conf_level, all=FALSE){
  
  #Determine m(a) for aLL.vec a<=K+1
  #For fixed a, m(a) is the smaLL.vecest b such that P(a<=X<=b)>=conf_level
  #so that AC(a,m(a)) is the core of rainbow, RB(a)
  #m(a)'s are neeeded b/c root1(a+1,m(a+1)) determines u(a)
  a <- 0; b <- 0
  m <- NA  # initialize vector of m(a)'s, the b for core(a)
  
  
  # Generate all core ACs -----------------------------------------------------------------------
  
  
  while(a <= ((K+1)+1)){ # we generate cores up to 1 more than K+1 to check for core skipping at every x.
    # Use curves from current rainbow {AC(a,b),b>=a} until core of next rainbow AC(a+1,m(a+1)) first 
    # rises above level. When transitioning between curves from same rainbow AC(a,b) to AC(a,b+1) 
    # transition when AC(a,b) faLL.vecs below level at root2(a,b) which determines lower limit for b+1
    while(test_coverage(a, b, conf_level) == F) {
      b <- b+1
    }
    
    m[a+1]=b
    a <- a+1
    b <- b+1 #can start the search at m(a)+1=b+1 b/c m(a+1)>=m(a)+1
  }
  
  
  
  # Start search from AC(0-0) -------------------------------------------------------------------
  LL.vec=rep(NA,K+1); UL.vec=rep(NA,K+1)
  a <- 0; b <- 0
  LL.vec[a+1] <- 0
  
  #Run until l(K) determined
  while(is.na(LL.vec[K+1])){
    while(find_roots(a,b,conf_level, root=2)<find_roots(a+1,m[a+2],conf_level, root=1) & is.na(LL.vec[K+1])){
      b=b+1
      LL.vec[b+1] <- find_roots(a,b-1, conf_level, root=2) 
    }
    
    #Exit loop once last needed lower limit (lower limit for x=n) is determined
    if(!is.na(LL.vec[K+1])){break}
    
    # Set Coincidental Endpoints ------------------------------------------------------------------
    
    # When transitioning between rainbows RB(a) to RB(a+1) we move from curve AC(a,b) to core of 
    # next rainbow AC(a+1, m(a+1)). The location of this transition occurs at root1(a+1,m(a+1)) 
    # and determines both the upper endpoint for a u(a) and the lower endpoints for b+1,...,m(a+1),
    # l(b+1)=...=l(m(a+1)).
    b.temp <- b
    a <- a+1
    b <- m[a+1]
    for(i in (b.temp+1):min(b,K)){
      LL.vec[i+1]=find_roots(a,m[a+1], conf_level, root=1)
    }
    UL.vec[(a-1)+1]=find_roots(a,m[a+1], conf_level, root=1)
  }
  
  
  
  # Set Remaining Upper Endpoints UL(x) ---------------------------------------------------------
  
  # Determine upper endpoints for remaining x; i.e. determine u(x) for a<=x<=n
  # B/c once we have lower endpoints for x up to n we can work on upper endpoints 
  # separately. These remaining values of u(x) are determined by u(x)=root1(x+1,m(x+1)). 
  skipped <- c() # There is also possibility of core skipping, which we will keep track of here.
  
  for(x in a:K){
    UL.vec[x+1]=find_roots(x+1,m[x+2], conf_level,root=1) # Set upper limit
    
    # SKIP CHECK: Check if next core is to be skipped.
    if (x<K & (UL.vec[x+1]>find_roots(x+2, m[x+3], conf_level,root=1))) {
      UL.vec[x+1] = find_roots(x+2, m[x+3], conf_level,root=1); skipped <- c(skipped,x)
    }
  }
  
  # Generate Results Table ----------------------------------------------------------------------
  CI <- data.frame(x=0:K,lower=LL.vec,upper=UL.vec)
  if(all==F){
    CI <- CI %>% 
      filter(x == K)
  }
  if(length(skipped) > 0) { #if there were any skips (rare)
    print(paste0("skips occurred for x=", paste(skipped, collapse = ",")))
  }
  
  return(CI)
}

###############################################################
#-------------------------------------------------------------#
# CMC EXAMPLES                                                #
#-------------------------------------------------------------#
###############################################################
# CMC.CI=CMC.pois(K=50,conf_level=.95,all=F); CMC.CI
# CMC.CI=CMC.pois(K=500,conf_level=.95,all=TRUE); CMC.CI



##-------------------------------------------------------------##
##                      Kabaila & Byrne (KB)                   ##
##-------------------------------------------------------------##


######################################################
#Poisson(theta)  
KB.pois <- function(x, conf_level=.95, all = FALSE){
  if (all == TRUE) { # indicates that user only wants to output interval for x=0 up to given x
    x <- 0:x
  }
  r <- function(x){s=1 ; while( pois((x-s):(x-1), max.pois(x-s,x-1))<=conf_level ){s=s+1}; return(s) }
  p <- function(x){q=1 ; while( pois((x+1):(x+q), max.pois(x+1,x+q))<=conf_level ){q=q+1}; return(q) }
  
  f <- function(lambda,x){return(pois((x-r(x)):(x-1),lambda)-conf_level)}
  g <- function(lambda,x){return(pois(x:(x+p(x)-1),lambda)-conf_level)}
  
  i=1
  l=x
  u=x
  
  if(x[1]==0){l[1] = 0; u[1] = uniroot(g, c(0,20),tol = 10^-10, x=x[1]) $root ; i=i+1}
  
  while(i<=length(x)){
    
    # uses 5*UL of scores method to get a good idea how far out uniroot should search for the root
    z = qnorm(1-(1-conf_level)/2,0,1); b = 4*(x[i] + (1/2)*z^2 + z*sqrt(x[i] + (1/4)*z^2))	
    
    a = max.pois((x[i]-r(x[i])),(x[i]-1)) 
    l[i] = uniroot(f, c(a,b),tol = 10^-10, x=x[i]) $root
    
    a = max.pois(x[i],(x[i]+p(x[i])-1))
    u[i] = uniroot(g, c(a,b),tol = 10^-10, x=x[i]) $root
    
    i=i+1
    
  }
  
  
  CIs <-  data.frame(x=x, lower=l, upper=u)
  
  return(CIs)
}

# KB Example
# ci.KB.pois=KB.pois(x=0:30, conf_level=.95)
# print(ci.KB.pois)

