# Length Comparisons

# Average Length (95% Confidence
# we'll calculate the average length up to K for each strict method

# Wald

# avg.length <- c()
# # for each K, calculate avg length
# for (K in 0:100) {
#   avg.length[K+1] <- (1/(K+1))*sum((wald_CI(K,0.95,all=T)$upper - wald_CI(K,0.95,all=T)$lower))
# }
# K <- 0:100
# plot(K, avg.length)

#################################################################
##       Calculating Average Length for Strict Procedures      ##
#################################################################

avgL_CP <- c(); avgL_OC <- c(); avgL_CG <- c(); avgL_CMC <- c(); #; avgL_B <- c()
# for each K, calculate avg length.
for (K in 0:150) {
  
  # CMC
  avgL_CMC[K+1] <- (1/(K+1))*sum((CMC(K,0.95,all=T)$upper - CMC(K,0.95,all=T)$lower))
  
  # Clopper-Pearson
  avgL_CP[K+1] <- (1/(K+1))*sum((clopperPearson_CI(K,0.95,all=T)$upper - clopperPearson_CI(K,0.95,all=T)$lower))
  avgL_CP[K+1] <- avgL_CP[K+1]/avgL_CMC[K+1]                        
  
  # OC/MST
  avgL_OC[K+1] <- (1/(K+1))*sum((OC(K,0.95,all=T)$upper - OC(K,0.95,all=T)$lower))
  avgL_OC[K+1] <- avgL_OC[K+1]/avgL_CMC[K+1]
  
  # Crow & Gardner
  avgL_CG[K+1] <- (1/(K+1))*sum((CG(K,0.95,all=T)$upper - CG(K,0.95,all=T)$lower))
  avgL_CG[K+1] <- avgL_CG[K+1]/avgL_CMC[K+1]
  
  # Blaker
  # avgL_B[K+1] <- (1/(K+1))*sum((blaker_CI(K,0.95,all=T)$upper - blaker_CI(K,0.95,all=T)$lower))
  # avgL_B[K+1] <- avgL_B[K+1]/avgL_CMC[K+1]
   
  
}

# create a dataframe with all running average lengths for each K = 0,...,100 [for plotting]
# each method has its own column

K <- 0:150
avg_lengths <- data.frame(K,
                          CP = avgL_CP,
                          OC = avgL_OC,
                          CG = avgL_CG
                          # ,B = avgL_B
                          )


# Now creating the Plot to visually compare strict methods' average lengths

avg_lengths_rel_CMC_plot <- avg_lengths |> 
  ggplot(mapping = aes(x = K)) +
  geom_line(mapping = aes(y = CP), color = "#00FF00", lwd = 1) + 
  geom_line(mapping = aes(y = OC), color = "#00AA55", lwd = 1) +
  geom_line(mapping = aes(y = CG), color = "#0055AA", lwd = 1) +
  scale_y_continuous(breaks = seq(0,1.5,0.1), limits = c(0.95,1.25)) +
  scale_x_continuous(limits = c(0,100)) +
  labs(
    title = "Strict Methods Comparison: Average Length for Integers 0 to K",
    subtitle = "Avg. Length Relative to CMC Method (2023)",
    x = "K",
    y = ""
  ) +
  theme_minimal(); avg_lengths_rel_CMC_plot

# save plot with fixed dimensions
# ggsave(filename = here::here("images",
#                                "avg-length-vs-K.png"),
#          width = 6,
#          height = 4,
#          units = "in",
#          dpi = 400)
  
  
  # creating a color palette
  pal <- colorRampPalette(c("green", "blue"))
  pal(4) 
  
  
  
#################################################################
##      Calculating Expected Length for Strict Procedures      ##
#################################################################
  
runtime_start <- Sys.time()
conf.level <- 0.95
# K = desired
K <- 122
lambda_vec <- seq(0,100,0.1) 
exp_width.CP.vec <- rep(0,length(lambda_vec)); exp_width.OC.vec <- rep(0,length(lambda_vec)) 
exp_width.CG.vec <- rep(0,length(lambda_vec));exp_width.CMC.vec <- rep(0, length(lambda_vec))
exp_width.B.vec <- rep(0,length(lambda_vec))



print(paste0("Code starting running at ", Sys.time()))

CI_lengths.CP <- clopperPearson_CI(K, conf.level, all=T)$upper - clopperPearson_CI(K, conf.level, all=T)$lower
CI_lengths.OC <- OC(K, conf.level, all=T)$upper - OC(K, conf.level, all=T)$lower
CI_lengths.CG <- CG(K, conf.level, all=T)$upper - CG(K, conf.level, all=T)$lower
# CI_lengths.B <- blaker_CI(K, conf.level, all=T)$upper - blaker_CI(K, conf.level, all=T)$lower
CI_lengths.CMC <- CMC(K, conf.level, all=T)$upper - CMC(K, conf.level, all=T)$lower

CI_lengths <- data.frame(CP = CI_lengths.CP,
                         CG = CI_lengths.CG,
                         OC = CI_lengths.OC,
                         CMC = CI_lengths.CMC)

  # for each lambda, calculate Expected Length
for (i in 1:length(lambda_vec)) {
  for (x in 0:K) {
    
    # Clopper-Pearson
    # current x's interval width for CP method
    length.CP <- CI_lengths$CP[x+1]
    exp_width.CP.vec[i] <- exp_width.CP.vec[i] + (dpois(x, lambda_vec[i]) * length.CP)
    
    # OC/MST
    length.OC <- CI_lengths$OC[x+1]
    exp_width.OC.vec[i] <- exp_width.OC.vec[i] + (dpois(x, lambda_vec[i]) * length.OC)
    
    # Crow & Gardner
    length.CG <- CI_lengths$CG[x+1]
    exp_width.CG.vec[i] <- exp_width.CG.vec[i] + (dpois(x, lambda_vec[i]) * length.CG)
    
    # Blaker
    # length.B <- CI_lengths$B[x+1]
    # exp_width.B.vec[i] <- exp_width.B.vec[i] + (dpois(x, lambda_vec[i]) * length.B)
    
    # CMC
    length.CMC <- CI_lengths$CMC[x+1]
    exp_width.CMC.vec[i] <-  exp_width.CMC.vec[i] + (dpois(x, lambda_vec[i]) * length.CMC)
    
  } 
}

expected_widths_table <- data.frame(lambda = lambda_vec,
                                    CP = exp_width.CP.vec,
                                    OC = exp_width.OC.vec,
                                    CG = exp_width.CG.vec,
                                    #, B = exp_width.B.vec,
                                    CMC = exp_width.CMC.vec)
runtime_end <- Sys.time()
print(paste0("Code starting running at ", runtime_start))
print(paste0("Code stopped running at ", runtime_end))


  
exp_width_rel_CMC_plot <- expected_widths_table |> 
  mutate(CP.rel = CP/CMC,
         OC.rel = OC/CMC,
         CG.rel = CG/CMC) |> 
  ggplot(mapping = aes(x = lambda)) +
  geom_line(mapping = aes(y = CP.rel), color = "#0000FF", lwd = 0.8) + 
  geom_line(mapping = aes(y = OC.rel), color = "#5555AA", lwd = 0.8) +
  geom_line(mapping = aes(y = CG.rel), color = "#AAAA55", lwd = 0.8) +
  scale_y_continuous(breaks = seq(0,1.5,0.1)) +

  labs(
    title = "Strict Methods Comparison: Expected Length",
    subtitle = "Relative Exp. Length to CMC Method (2023)",
    x = expression(lambda),
    y = ""
  ) + 
  theme_minimal(); exp_width_rel_CMC_plot
exp_width_rel_CMC_plot


  
  
  
  
  
  
  
  
  
  
  
  
  
  
  