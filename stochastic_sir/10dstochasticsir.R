rm(list=ls())
library(deSolve)
library(doParallel)
library(foreach)
library(truncnorm)
library(odin)
library(abc)
library(densratio)
library(rBayesianOptimization)
library(ggplot2)

start<-Sys.time()
# Parallel setup (use one fewer than available cores to keep UI responsive)
totalCores <- parallel::detectCores()
ncores <- max(1, totalCores - 1)
cluster <- makeCluster(ncores)
registerDoParallel(cluster)

# ----------------------------
# Problem constants (same as yours)
# ----------------------------
S_init <- 490
I_init <- 10
R_init <- 0
dt <- 0.1
steps <- 30                     # your simulate_sir uses steps -> gives length steps+1
time_indices <- 1:steps   # indices we will allow for observation
# note: indices correspond to rows of simulate_sir output (1 .. steps+1)

# ----------------------------
# SIR simulator (kept same, but remove set.seed(123) inside loop to avoid identical draws each step)
# ----------------------------
simulate_sir <- function(beta, gamma, S_init, I_init, R_init, dt, steps, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  S <- numeric(steps + 1)
  I <- numeric(steps + 1)
  R <- numeric(steps + 1)
  S[1] <- S_init; I[1] <- I_init; R[1] <- R_init
  for (i in 1:steps) {
    delta_I <- rbinom(1, S[i], beta * I[i] / (S[i] + I[i] + R[i]))
    delta_R <- rbinom(1, I[i], gamma)
    S[i + 1] <- S[i] - delta_I
    I[i + 1] <- I[i] + delta_I - delta_R
    R[i + 1] <- R[i] + delta_R
  }
  return(list(S = S, I = I, R = R))
}

# ----------------------------
# Prior draws (kept similar but can be reduced for speed)
# ----------------------------
n_prior <- 1000
set.seed(1)
beta <- runif(n_prior, 0, 0.5)
gamma <- runif(n_prior, 0, 0.5)
priorabc <- cbind(beta, gamma)

# Precompute prior forward trajectories (same idea as your sir_results)
sir_results <- lapply(1:nrow(priorabc), function(i) {
  as.data.frame(simulate_sir(priorabc[i, 1], priorabc[i, 2], S_init, I_init, R_init, dt, steps, seed = 1000 + i))
})

# Precompute z_recov_prior = sum of R at indices 3:5 as you had
z_recov_prior <- sapply(1:nrow(priorabc), function(i) sum(sir_results[[i]]$R[3:5]))

# for each prior draw j, store R at all indices
# We'll construct a matrix R_prior of size n_prior x (steps+1)
R_prior_mat <- sapply(sir_results, function(df) df$R)  # columns correspond to prior draws
R_prior_mat <- t(R_prior_mat)                          # now rows = priors, cols = time indices (1..steps+1)

# ----------------------------
# Helper: make summary vector for a triple of indices (we use R only as in earlier example)
# Returns matrix n_prior x 3 (R at t1,t2,t3 for each prior)
# ----------------------------
get_sim_summary_R <- function(indices, R_prior_mat) {
  # indices: integer vector length 3 with values in 1:(steps+1)
  idx <- pmin(pmax(round(indices), 1), ncol(R_prior_mat))
  # extract columns
  mat <- R_prior_mat[, idx, drop = FALSE]  # n_prior x 3
  return(mat)
}

# ----------------------------
# Utility function for a triple (t1,t2,t3)
# Keep ABC + densratio pipeline as you wrote, but sample subset of i's to reduce cost
# Returns numeric average log density ratio (like your ut)
# ----------------------------
compute_utility_triple <- function(t1, t2, t3,t4,t5,t6,t7,t8,t9,t10,
                                   priorabc, R_prior_mat, z_recov_prior,
                                   S_init, I_init, R_init, dt, steps,
                                   tol_abc = 0.2,
                                   sample_size_i = 200,  # number of i's to loop over (lower => faster)
                                   seed_base = 1) {
  # Map to integer indices in [1, steps+1]
  idx <- pmin(pmax(round(c(t1, t2, t3,t4,t5,t6,t7,t8,t9,t10)), 1), steps + 1)
  # summary matrix: n_prior x 3 (each row: R at t1,t2,t3 for that prior)
  sumstat_mat <- get_sim_summary_R(idx, R_prior_mat)  # n_prior x 3
  
  n_prior <- nrow(priorabc)
  # choose subset of i's to loop over (like you did each i as a pseudo-observed)
  sample_i <- seq_len(n_prior)
  
  # Preallocate vector to hold density ratio estimates per sampled i
  densratio_vals <- rep(NA, length(sample_i))
  
  # Loop (parallelize with foreach across sampled i's)
  densratio_vals <- foreach(ii = seq_along(sample_i),
                            .combine = c,
                            .packages = c("abc", "densratio"),
                            .export = c("simulate_sir", "get_sim_summary_R")) %dopar% {
                              i <- sample_i[ii]
                              tryCatch({
                                # observed summary for this i (vector length 3)
                                obs_summary <- sumstat_mat[i, ]
                                
                                # Run abc: target = obs_summary, param = priorabc, sumstat = sumstat_mat
                                # Keep method and transforms from your original code
                                abc_res <- abc::abc(target = obs_summary,
                                                    param = priorabc,
                                                    sumstat = sumstat_mat,
                                                    tol = tol_abc,
                                                    method = "neuralnet",
                                                    transf = c("logit", "logit"),
                                                    logit.bounds = rbind(c(0,0.5), c(0,0.5)))
                                if (is.null(abc_res) || is.null(abc_res$adj.values) || nrow(abc_res$adj.values) == 0) {
                                  return(NA_real_)
                                }
                                adj_vals <- abc_res$adj.values
                                if (nrow(adj_vals) < 2) return(NA_real_)
                                
                                # compute z.abci: for each adjusted param draw, simulate and compute sum of R[3:5]
                                z_abci <- numeric(nrow(adj_vals))
                                for (l in 1:nrow(adj_vals)) {
                                  sdat <- simulate_sir(adj_vals[l,1], adj_vals[l,2], S_init, I_init, R_init, dt, steps, seed = seed_base + 1000 + l)
                                  z_abci[l] <- sum(sdat$R[3:5])
                                }
                                # generate synthetic thetai from prior-like draws, same shape, and compute z.pi
                                set.seed(seed_base + i + 5000)
                                thetai <- cbind(runif(nrow(adj_vals), 0, 0.5), runif(nrow(adj_vals), 0, 0.5))
                                z_pi <- numeric(nrow(thetai))
                                for (p in 1:nrow(thetai)) {
                                  sdat2 <- simulate_sir(thetai[p,1], thetai[p,2], S_init, I_init, R_init, dt, steps, seed = seed_base + 2000 + p)
                                  z_pi[p] <- sum(sdat2$R[3:5])
                                }
                                
                                # densratio: try-catch robust
                                dr_obj <- tryCatch({
                                  densratio::densratio(z_abci, z_pi, method = "RuLSIF",
                                                       alpha = 0.5,
                                                       sigma = 10 ^ seq(-5, 10, length.out = 5),
                                                       lambda = 10 ^ seq(-5, 10, length.out = 5))
                                }, error = function(e) NULL)
                                if (is.null(dr_obj)) return(NA_real_)
                                
                                # compute density ratio evaluated at z_recov_prior (same indexing as original)
                                densvals <- tryCatch({
                                  v <- dr_obj$compute_density_ratio(z_recov_prior)
                                  if (length(v) < 1) return(NA_real_)
                                  # pick the i-th element like you did: rati[k,i] <- densratio_obj$compute_density_ratio(z_recov_prior)[i]
                                  as.numeric(v[i])
                                }, error = function(e) NA_real_)
                                return(densvals)
                              }, error = function(e) {
                                return(NA_real_)
                              })
                            } # end foreach
  
  # densratio_vals is vector of density-ratio estimates for the sampled i's
  # convert to log space like you did in original code (with NA handling)
  
  lr <- log(densratio_vals)
  lr[!is.finite(lr)] <- NA
  if (all(is.na(lr))) return(NA_real_)
  return(mean(lr, na.rm = TRUE))   # same as ut = rowMeans(log(rat), na.rm=TRUE) conceptually
}



# cat("Test utility eval (single triple):\n")
# test_val <- compute_utility_triple(2, 8, 14,1,4,3,7,9,10,6,
#                                    priorabc, R_prior_mat, z_recov_prior,
#                                    S_init, I_init, R_init, dt, steps,
#                                    tol_abc = 0.2, sample_size_i = n_prior, seed_base = 42)
# print(test_val)

# ----------------------------
# Bayesian optimization wrapper (using rBayesianOptimization)
# We'll scale t inputs into [1, steps+1] domain. BO will propose continuous values; we round inside compute_utility_triple.
# ----------------------------
# Wrapper expecting named args t1, t2, t3
bo_fun <- function(t1, t2, t3,t4,t5,t6,t7,t8,t9,t10) {
  # Do a small average (reps) to reduce noise; optionally increase reps
  reps <- 2
  vals <- numeric(reps)
  for (r in 1:reps) {
    vals[r] <- compute_utility_triple(t1, t2, t3,t4,t5,t6,t7,t8,t9,t10,
                                      priorabc, R_prior_mat, z_recov_prior,
                                      S_init, I_init, R_init, dt, steps,
                                      tol_abc = 0.2,
                                      sample_size_i = 1000,    # tradeoff: higher->better, slower
                                      seed_base = 100 * r)
  }
  # rBayesianOptimization expects Score; we maximize Score
  # If utility is NA or -Inf, return a very small number to avoid crashing the GP
  if (all(is.na(vals))) {
    return(list(Score = -1e10))
  } else {
    return(list(Score = mean(vals, na.rm = TRUE)))
  }
}

# Bounds: allow continuous proposals between 1 and steps+1
bounds <- list(t1 = c(1, steps), t2 = c(1, steps), t3 = c(1, steps),
               t4 = c(1, steps), t5 = c(1, steps), t6 = c(1, steps),
               t7 = c(1, steps), t8 = c(1, steps), t9 = c(1, steps),
               t10=c(1, steps))

# Run BO (tune init_points and n_iter to taste; keep small here)
set.seed(12345)
bo_res <- BayesianOptimization(FUN = bo_fun,
                               bounds = bounds,
                               init_points = 10,   # initial random evals
                               n_iter = 1,       # sequential BO iterations
                               acq = "ucb",
                               kappa = 2.576,
                               verbose = TRUE)
end<-Sys.time()
print(end-start)

print(bo_res)
# Best Parameters Found: 
#t1 = 23.06849	t2 = 6.443661	t3 = 29.59353	
#t4 = 3.285587	t5 = 3.601644	t6 = 24.07446	
#t7 = 6.691329	t8 = 11.95243	t9 = 27.29791	t10 = 5.523388	Value = 0.2163195 


# Stop cluster
stopCluster(cluster)
