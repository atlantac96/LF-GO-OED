rm(list=ls())
library(deSolve)
library(doParallel)
totalCores <- detectCores()
cluster <- makeCluster(totalCores[1]) 
registerDoParallel(cluster)
library(foreach)
library(truncnorm)
library(odin)

design<- seq(0,3,by=0.1)
num_replicates <- 5
ut_matrix <- matrix(NA, nrow = length(design), ncol = num_replicates)
S_init <- 490
I_init <- 10
R_init <- 0
dt <- 0.1
steps <- 30

for(replicate_num in 1:num_replicates){
  
  library(deSolve)
  library(doParallel)
  # simulate_sir <- function(beta, gamma) {
  #   sir_generator <- odin::odin({
  #     ## Core equations for transitions between compartments:
  #     update(S) <- S - n_SI
  #     update(I) <- I + n_SI - n_IR
  #     update(R) <- R + n_IR
  #     
  #     ## Individual probabilities of transition:
  #     p_SI <-  beta*I/N#1 - exp(-beta * I / N) # S to I
  #     p_IR <- gamma#1 - exp(-gamma) # I to R
  #     
  #     ## Draws from binomial distributions for numbers changing between
  #     ## compartments:
  #     n_SI <- rbinom(S, p_SI)
  #     n_IR <- rbinom(I, p_IR)
  #     
  #     ## Total population size
  #     N <- S + I + R
  #     
  #     ## Initial states:
  #     initial(S) <- S_ini
  #     initial(I) <- I_ini
  #     initial(R) <- 0
  #     
  #     ## User defined parameters - default in parentheses:
  #     S_ini <- user(490)
  #     I_ini <- user(10)
  #     beta <- user(0)
  #     gamma <- user(0)
  #     
  #   }, verbose = FALSE)
  #   
  #   sir <- sir_generator$new(I_ini = 10, beta = beta, gamma = gamma)
  #   res <- sir$run(1:50) 
  #   return(res)
  # }
  simulate_sir<- function(beta, gamma, S_init, I_init, R_init, dt, steps) {
    # Initialize vectors to store S, I, R
    S <- numeric(steps + 1)
    I <- numeric(steps + 1)
    R <- numeric(steps + 1)
    
    # Set initial values
    S[1] <- S_init
    I[1] <- I_init
    R[1] <- R_init
    for (i in 1:steps) {
      # Calculate increments with binomial distributions
      set.seed(123)
      delta_I <- rbinom(1, S[i], beta * I[i] / (S[i] + I[i] + R[i]))
      delta_R <- rbinom(1, I[i], gamma)
      
      # Update S, I, R using Euler's method
      S[i + 1] <- S[i] - delta_I
      I[i + 1] <- I[i] + delta_I - delta_R
      R[i + 1] <- R[i] + delta_R
    }
    
    
    # Return the result
    return(list(S = S, I = I, R = R))
  }
  
  cat("Replicate:", replicate_num, "\n")
  
  # Set a different seed for each replicate
  set.seed(replicate_num)
  
  # Priors
  beta = runif(1000, 0, 0.5) ##infection parameter
  gamma = runif(1000, 0, 0.5) ##recovery parameter
  parameters_values <- priorabc<-cbind(beta, gamma)
  sir_results <- lapply(1:nrow(priorabc), function(i) {
    as.data.frame(simulate_sir(priorabc[i, "beta"], priorabc[i, "gamma"], 
                               S_init, I_init, R_init, dt, steps))
  })
  design<- seq(0,3, by=0.1)
  
  
  trans <- vector("list", length = length(design))
  result_matrix <- matrix(NA, nrow = nrow(priorabc), ncol = 3)
  
  for (s in 1: length(design)) {
    for (i in 1:nrow(priorabc)){
      sir_res <- sir_results[[i]][s,]
      result_matrix[i, 1] <-sir_res[1,1]
      result_matrix[i, 2] <-sir_res[1,2]
      result_matrix[i, 3] <-sir_res[1,3]
    }
    trans[[s]] <- result_matrix
  }
  
  
  
  #goals: want recovered at 3rd and 4th, 5th  time instant
  z=matrix(NA, nrow=nrow(parameters_values), ncol=1)
  for(i in 1:nrow(parameters_values)){
    sirdata<-as.data.frame(simulate_sir(priorabc[i, "beta"], priorabc[i, "gamma"], 
                                        S_init, I_init, R_init, dt, steps))
    z[i] <- sum(sirdata[3:5,3])
  }
  z_recov_prior<-z
  z.abc<-list()
  
  
  
  rat <- rati <- matrix(NA, nrow = length(design), ncol = nrow(parameters_values))
  for (k in 2:length(design)) {
    cat(" Starting Replicate:", replicate_num, "with design point",k,
        "at", print(Sys.time()), "\n")
    z.abc[[k]] <- list()
    rat[k, ] <- foreach(i = 1:nrow(parameters_values), .combine = cbind) %dopar% {
      #if (trans[[k]][i, 3] == 0) {
      #rati[k,i]<-NA      
      #  }else{
      library(deSolve)
      library(truncnorm)
      library(abc)
      check<-abc(trans[[k]][i,1:3],priorabc, trans[[k]][,1:3],
                 tol=0.2,method="neuralnet", transf=c("logit", "logit"),
                 logit.bounds=rbind(c(0,0.5), c(0,0.5)))
      check$adj.values<- check$adj.values[complete.cases(check$adj.values), ]
      z.abci=matrix(NA, nrow=nrow(check$adj.values), ncol=1)
      for(l in 1:nrow(check$adj.values)){
        sirdata<-as.data.frame(simulate_sir(check$adj.values[l, 1], check$adj.values[l,2], 
                                            S_init, I_init, R_init, dt, steps))
        z.abci[l]<-sum(sirdata[3:5,3])
      }
      z.abc[[k]][[i]]<- z.abci
      

      set.seed(i)
      thetai <- matrix(NA, nrow = nrow(check$adj.values), ncol = 2)
      thetai[, 1] <- runif(nrow(check$adj.values), 0, 0.5)
      thetai[, 2] <- runif(nrow(check$adj.values), 0, 0.5)
      colnames(thetai) <- c("beta", "gamma")
      z.pi=matrix(NA, nrow=nrow(check$adj.values), ncol=1)
      for(p in 1:nrow(check$adj.values)){
        sirdata<-as.data.frame(simulate_sir(thetai[p, 1], thetai[p, 2], 
                                            S_init, I_init, R_init, dt, steps))
        z.pi[p]<-sum(sirdata[3:5,3])
      }
      zi<-z.pi
      
      library(densratio)
      set.seed(56)
      densratio_obj <- densratio::densratio(z.abc[[k]][[i]], zi, method = "RuLSIF", alpha = 0.5,
                                            sigma = 10 ^ seq(-5, 10, length.out = 10),
                                            lambda = 10 ^ seq(-5, 10, length.out = 10))
      rati[k, i] <- densratio_obj$compute_density_ratio(z_recov_prior)[i]
      # }
    }
  }
  ratna <- log(rat)
  ratna[!is.finite(ratna)] <- NA
  
  ut <- rowMeans(ratna, na.rm = TRUE)
  ut_matrix[, replicate_num] <- ut
}

save(ut_matrix, file="utmatrixstochasticsir_goal_recovered.RDa")

stopCluster(cluster)
library(ggplot2)
load("utmatrixstochasticsir.RDa")
mean_row <- rowMeans(ut_matrix, na.rm=TRUE)
sd_row <- apply(ut_matrix, 1, sd, na.rm=TRUE)
df <- data.frame(sequence = 1:nrow(ut_matrix), mean = mean_row, sd = sd_row)
pdf("sirgoalrecov.pdf", height=4, width=8)
ggplot(data=df, aes(x=design, y=mean)) + 
  labs(x = expression(paste(italic("d"))), 
       y = bquote(paste(italic(U[z](d)))), 
       color = "") +  theme_minimal() +
  geom_line() + theme(text=element_text(size=30))+
  geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd), alpha = 0.2)
dev.off()
design[which.max(df$mean)] #0.2
df$mean[which.max(df$mean)] #0.638
design[which.min(df$sd)]
df$sd[which.min(df$sd)]
df$sd[which.max(df$sd)]

data_df <- data.frame(Column1 = ut_matrix[, 1],
                      Column2 = ut_matrix[, 2],
                      Column3 = ut_matrix[, 3],
                      Column4 = ut_matrix[, 4])
library(RColorBrewer)
# Create a plot using ggplot2
library(ggplot2)
ggplot(data_df, aes(x = design)) +
  geom_line(aes(y = Column1, color = "Replicate 1"), size = 1) +
  geom_line(aes(y = Column2, color = "Replicate 2"), size = 1) +
  geom_line(aes(y = Column3, color = "Replicate 3"), size = 1) +
  geom_line(aes(y = Column4, color = "Replicate 4"), size = 1) +
  scale_color_manual(
    values = c(
      "Replicate 1" = "red",
      "Replicate 2" = "blue",
      "Replicate 3" = "green",
      "Replicate 4" = "purple"))+
  
  # "Replicate 6" = "pink"
  
  labs(
    x = expression(paste("Design, ", italic("d"))),
    y = bquote(paste("Utility, ", italic(U[z](d)))),
    color = ""
  ) +
  theme_minimal() +
  theme(legend.position = "top") +
  coord_cartesian(ylim = c(min(ut_matrix[, 1:4]), max(ut_matrix[, 1:4]) * 2.1)) +
  scale_y_continuous(
    breaks = seq(floor(min(ut_matrix[-1, 1:4])), ceiling(max(ut_matrix[-1, 1:4]) * 2.1), by = 0.002),
    labels = scales::number_format(accuracy = 0.001)
  )

design[which.max(ut)]
ut[which.max(ut)]

