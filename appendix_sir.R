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
num_replicates <- 1
ut_matrix <- matrix(NA, nrow = length(design), ncol = num_replicates)
S_init <- 498
I_init <- 2
R_init <- 0
dt <- 0.1
steps <- 30


for(replicate_num in 1:num_replicates){
  
  library(deSolve)
  library(doParallel)

  simulate_sir<- function(beta, gamma, S_init, I_init, R_init, dt, steps) {
    # Initialize vectors to store S, I, R
    S <- numeric(steps)
    I <- numeric(steps)
    R <- numeric(steps)
    
    # Set initial values
    S[1] <- S_init
    I[1] <- I_init
    R[1] <- R_init
    for (i in 2:steps) {
      # Calculate increments with binomial distributions
      set.seed(i)
      delta_I <-  beta * I[i-1]*S[i-1]/(S[i-1]+I[i-1]+R[i-1])
      delta_R <- gamma*I[i-1]

      
      # Update S, I, R using Euler's method
      S[i ] <- S[i-1] - delta_I
      I[i ] <- rpois(1,0.95*(I[i-1] + delta_I - delta_R))
      R[i] <- R[i-1] + delta_R
    }
    return(I=I)
  }
  
  library(jsonlite)
  
  # Load the JSON file into R
  data <- fromJSON("data.json")
  
  
  cat("Replicate:", replicate_num, "\n")
  
  # Set a different seed for each replicate
  set.seed(replicate_num)
  
  # Priors
  beta = rlnorm(n = 1000, meanlog = 0.50, sdlog = 0.50)
  gamma = rlnorm(n = 1000, meanlog = 0.10, sdlog = 0.50)
  parameters_values <- priorabc<-cbind(beta, gamma)
  
  z_recov_prior <- parameters_values
  sir_results <- lapply(1:nrow(priorabc), function(i) {
    as.data.frame(simulate_sir(priorabc[i,1], priorabc[i,2], 
                               S_init, I_init, R_init, dt, steps))
  })
  design<- seq(0,3, by=0.1)
  
  trans <- vector("list", length = length(design))
  result_matrix <- matrix(NA, nrow = nrow(priorabc), ncol = 1)
  
  for (s in 1: length(design)) {
    for (i in 1:nrow(priorabc)){
      sir_res <- sir_results[[i]][s,]
      result_matrix[i, 1] <-sir_res
    }
    trans[[s]] <- result_matrix
  }
  
  
  z.abc <- list()
  rat <- rati <- matrix(NA, nrow = length(design), ncol = nrow(parameters_values))
  for (k in 2:length(design)) {
    cat(" Starting Replicate:", replicate_num, "with design point",k,
        "at", print(Sys.time()), "\n")
    z.abc[[k]] <- list()
    rat[k, ] <- foreach(i = 1:nrow(parameters_values), .combine = cbind) %dopar% {
      
      library(deSolve)
      library(truncnorm)
      library(abc)
      check<-abc(trans[[k]][i,1],priorabc, trans[[k]][,1],
                 method="neuralnet", tol=0.1)
               
      z.abc[[k]][[i]] <- check$adj.values[complete.cases(check$adj.values), ]
      
      set.seed(i)
      thetai <- matrix(NA, nrow = nrow(check$unadj.values), ncol = 2)
      thetai[, 1] =  rlnorm(nrow(check$unadj.values), meanlog = 0.50, sdlog = 0.50)
      thetai[, 2] = rlnorm(nrow(check$unadj.values), meanlog = 0.10, sdlog = 0.50)
      colnames(thetai) <- c("beta", "gamma")
      zi <- thetai[complete.cases(thetai), ]
      
      library(densratio)
      set.seed(56)
      densratio_obj <- densratio::densratio(z.abc[[k]][[i]], zi, method = "RuLSIF",
                                          #  alpha=0,
                                           alpha = 0.5,
                                        sigma = 10 ^ seq(-5, 1, length.out = 10),
                                         lambda = 10 ^ seq(-5, 1, length.out = 10))
      rati[k, i] <- densratio_obj$compute_density_ratio(z_recov_prior)[i]
      # }
    }
  }
  ratna <- log(rat)
  ratna[!is.finite(ratna)] <- NA
  
  ut <- rowMeans(ratna, na.rm = TRUE)
  ut_matrix[, replicate_num] <- ut
}
library(ggplot2)
pdf("appendixsir.pdf", height=4, width=8)
plot(design, ut, type = "b",cex.lab=1.2, cex.axis=1, cex.main=1, cex.sub=1,
     xlab = expression(paste("Design, ", italic("d"))),
     ylab = bquote(paste("Utility, ", italic(U[z](d))))) 
dev.off()

plot(design,ut)
save(ut_matrix, file="utmatrixstochasticsirparaminf.RDa")

stopCluster(cluster)
library(ggplot2)
load("utmatrixstochasticsirparaminf.RDa")
mean_row <- rowMeans(ut_matrix, na.rm=TRUE)
sd_row <- apply(ut_matrix, 1, sd, na.rm=TRUE)
df <- data.frame(sequence = 1:nrow(ut_matrix), mean = mean_row, sd = sd_row)
pdf("sirparaminf.pdf", height=4, width=8)
ggplot(data=df, aes(x=design, y=mean)) + 
  labs(x = expression(paste("Design, ", italic("d"))), 
       y = bquote(paste("Utility, ", italic(U[z](d)))), 
       color = "") +  theme_minimal() +
  geom_line() + theme(text=element_text(size=20))+
  geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd), alpha = 0.2)
dev.off()
design[which.max(df$mean)] #0.4
df$mean[which.max(df$mean)] #2.066
design[which.min(df$sd)]
df$sd[which.min(df$sd)]
df$sd[which.max(df$sd)]


#lfire:
design<- seq(0,3,by=0.1)
S_init <- 49
I_init <- 1
R_init <- 0
dt <- 0.01
steps <- 30
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

set.seed(123)
beta = runif(1000, 0, 0.5) ##infection parameter
gamma = runif(1000, 0, 0.5) ##recovery parameter
parameters_values <- priorabc<-cbind(beta, gamma)

z_recov_prior <- parameters_values
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


z.abc <- list()
rat <- rati <- matrix(NA, nrow = length(design), ncol = nrow(parameters_values))
for (k in 2:length(design)) {
  cat(" Starting Replicate:", replicate_num, "with design point",k,
      "at", print(Sys.time()), "\n")
  z.abc[[k]] <- list()
  rat[k, ] <- foreach(i = 1:nrow(parameters_values), .combine = cbind) %dopar% {
    library(truncnorm)
    check<-abc(trans[[k]][i,1:3],priorabc, trans[[k]][,1:3],
               tol=0.6,method="neuralnet", transf=c("logit", "logit"),
               logit.bounds=rbind(c(0,0.5), c(0,0.5)))
    z.abc[[k]][[i]] <- check$adj.values[complete.cases(check$adj.values), ]
    
    set.seed(i)
    thetai <- matrix(NA, nrow = nrow(check$adj.values), ncol = 2)
    thetai[, 1] <- runif(nrow(check$adj.values), 0, 0.5)
    thetai[, 2] <- runif(nrow(check$adj.values), 0, 0.5)
    colnames(thetai) <- c("beta", "gamma")
    
    
  }
}
ratna <- log(rat)
ratna[!is.finite(ratna)] <- NA

ut <- rowMeans(ratna, na.rm = TRUE)
ut_matrix[, replicate_num] <- ut
plot(design, ut_matrix[,1])


