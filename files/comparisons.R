# Length Comparisons

# Average Length (95% Confidence)
# we'll calculate the average length up to K for each strict method
library(tidyverse)

#################################################################
##       Calculating Average Length for Strict Procedures      ##
#################################################################

avgL_CP <- c(); avgL_OC <- c(); avgL_CG <- c(); avgL_CMC <- c(); avgL_B <- c()
i <- 1
# for each K, calculate avg length.
for (K in 10:100) {
  
  # CMC
  avgL_CMC[i] <- (1/(K+1))*sum((CMC.pois(K,0.95,all=T)$upper - CMC.pois(K,0.95,all=T)$lower))
  
  # Clopper-Pearson
  avgL_CP[i] <- (1/(K+1))*sum((ClopperPearson.pois(K,0.95,all=T)$upper - ClopperPearson.pois(K,0.95,all=T)$lower))
  
  # OC/MST
  avgL_OC[i] <- (1/(K+1))*sum((OC.pois(K,0.95,all=T)$upper - OC.pois(K,0.95,all=T)$lower))
  
  
  # Crow & Gardner
  avgL_CG[i] <- (1/(K+1))*sum((CG.pois(K,0.95,all=T)$upper - CG.pois(K,0.95,all=T)$lower))

  
  # Blaker
  avgL_B[i] <- (1/(K+1))*sum((Blaker.pois(K,0.95,all=T)$upper - Blaker.pois(K,0.95,all=T)$lower))
  
  i <- i + 1
  
}



# create a dataframe with all running average lengths for each K = 10,...,100 [for plotting]
# each method has its own column

avg_lengths <- data.frame(K = 10:100,
                          CP = avgL_CP,
                          OC = avgL_OC,
                          CMC = avgL_CMC,
                          CG = avgL_CG,
                          B = avgL_B
                          )


# Now creating the Plot to visually compare strict methods' average lengths

avg_lengths |> 
  ggplot(mapping = aes(x = K)) +
  geom_line(mapping = aes(y = CP, color = "#154734"), lwd = 1) +  # CP
  geom_line(mapping = aes(y = OC, color = "#3A913F"), lwd = 1) + # OC
  geom_line(mapping = aes(y = CMC, color = "#A4D65E"), lwd = 1) + # CMC
  geom_line(mapping = aes(y = CG, color = "#5CB8B2"), lwd = 1) + # CG
  geom_line(mapping = aes(y = B, color = "#ABCAE9"), lwd = 1) + # Blaker
  scale_y_continuous(breaks = seq(0,30,5)) +
  scale_x_continuous(limits = c(10,100)) +
  labs(
    title = "Strict Methods Comparison: Average Length for Integers 0 to K",
    # subtitle = "Avg. Length Relative to CMC Method (2023)",
    subtitle = "Avg. Length",
    x = "K",
    y = ""
  ) +
  scale_colour_manual(name = "method", guide = "legend",
                      values =c('#154734'='#154734','#3A913F'='#3A913F',"#A4D65E"= "#A4D65E",
                                '#5CB8B2'='#5CB8B2', '#ABCAE9'='#ABCAE9'),
                      labels =c("CP","OC","CMC", "CG", "B") ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks = element_line(size = 0.5),
    axis.ticks.length = unit(2, "mm"),
    axis.text.x = element_text(margin = margin(t = 5)),
    axis.text.y = element_text(margin = margin(r = 5))) 


ggsave(
  filename = here::here("images",
                        "avg_lengths.png"), # file path
  width = 6,
  height = 4,
  units = "in",
  dpi = 150, # change to dpi = 400 for better quality.
  bg = "white"
)
  
  

  
  
  
#################################################################
##      Calculating Expected Length for Strict Procedures      ##
#################################################################
  

conf.level <- 0.95
# K = desired
K <- 122
lambda_vec <- seq(0,100,0.1) 
exp_width.CP.vec <- rep(0,length(lambda_vec)); exp_width.OC.vec <- rep(0,length(lambda_vec)) 
exp_width.CG.vec <- rep(0,length(lambda_vec));exp_width.CMC.vec <- rep(0, length(lambda_vec))
exp_width.B.vec <- rep(0,length(lambda_vec))




CI_lengths.CP <- ClopperPearson.pois(K, conf.level, all=T)$upper - ClopperPearson.pois(K, conf.level, all=T)$lower
CI_lengths.OC <- OC.pois(K, conf.level, all=T)$upper - OC.pois(K, conf.level, all=T)$lower
CI_lengths.CG <- CG.pois(K, conf.level, all=T)$upper - CG.pois(K, conf.level, all=T)$lower
CI_lengths.B <- Blaker.pois(K, conf.level, all=T)$upper - Blaker.pois(K, conf.level, all=T)$lower
CI_lengths.CMC <- CMC.pois(K, conf.level, all=T)$upper - CMC.pois(K, conf.level, all=T)$lower

CI_lengths <- data.frame(CP = CI_lengths.CP,
                         CG = CI_lengths.CG,
                         OC = CI_lengths.OC,
                         CMC = CI_lengths.CMC,
                         B = CI_lengths.B)

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
    length.B <- CI_lengths$B[x+1]
    exp_width.B.vec[i] <- exp_width.B.vec[i] + (dpois(x, lambda_vec[i]) * length.B)
    
    # CMC
    length.CMC <- CI_lengths$CMC[x+1]
    exp_width.CMC.vec[i] <-  exp_width.CMC.vec[i] + (dpois(x, lambda_vec[i]) * length.CMC)
    
  } 
}

expected_widths_table <- data.frame(lambda = lambda_vec,
                                    CP = exp_width.CP.vec,
                                    OC = exp_width.OC.vec,
                                    CG = exp_width.CG.vec,
                                     B = exp_width.B.vec,
                                    CMC = exp_width.CMC.vec)



  
exp_width_rel_CMC_plot <- expected_widths_table |> 
  
  # divide each method by CMC length
  mutate(
    CP.rel = CP/CMC,
    OC.rel = OC/CMC,
    B.rel = B/CMC, 
    CG.rel = CG/CMC
  ) |> 
  
  ggplot(mapping = aes(x = lambda)) +
  geom_line(mapping = aes(y = CP.rel, color = "#BD8A12") , lwd = 0.8) + 
  geom_line(mapping = aes(y = OC.rel, color = "#f2c75c"), lwd = 0.8) +
  geom_line(mapping = aes(y = B.rel, color = "#3A913F"), lwd = 0.8) +
  geom_line(mapping = aes(y = CG.rel, color = "#154533"), lwd = 0.8) +
  scale_y_continuous(breaks = seq(0,1.5,0.1)) +

  labs(
    title = "Strict Methods Comparison: Expected Length",
    subtitle = "Relative Exp. Length to CMC Method (2023)",
    x = expression(lambda),
    y = ""
  ) +
  scale_color_manual(
    name = "Method",
    guide = "legend",
    values = c('#BD8A12'='#BD8A12','#f2c75c'='#f2c75c','#3A913F'='#3A913F', '#154533'='#154533'),
    labels = c('CG','B','CP','MST')
  ) +
  theme_minimal() +
  theme(
        # legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.ticks = element_line(size = 0.5),
        axis.ticks.length = unit(2, "mm"),
        axis.text.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_text(margin = margin(r = 5))
  )

# plot
exp_width_rel_CMC_plot

# save plot to images folder
ggsave(filename = here::here("images",
                             "exp_length_rel_CMC.png"),
       width = 6,
       height = 4,
       units = "in",
       dpi = 200, # change to 400 dpi for better quality
       bg = "white"
)
  
  
