rm(list = ls())

# --- libs -------------------------------------------------------------------
# required_pkgs <- c("deSolve","truncnorm","abc","densratio","rBayesianOptimization",
#                    "doParallel","foreach")
# for (p in required_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(deSolve)
library(truncnorm)
library(abc)
library(densratio)
library(rBayesianOptimization)
library(doParallel)
library(foreach)

start<-Sys.time()
# --- parallel setup (used optionally by you) --------------------------------
totalCores <- parallel::detectCores()
cl <- makeCluster(max(1, totalCores - 1))   # leave 1 core for system
registerDoParallel(cl)

# --- piecewise stimulus -----------------------------------------------------
xi_time <- function(t, zetas, T) {
  # zetas: numeric vector length 3
  seg_len <- T / 3
  if (t < seg_len) {
    return(zetas[1])
  } else if (t < 2 * seg_len) {
    return(zetas[2])
  } else {
    return(zetas[3])
  }
}

# --- FHN Euler-Maruyama SDE simulator --------------------------------------
fitzhugh_nagumo <- function(time, y0, parms, zetas) {
  n <- length(time) - 1
  dt <- diff(time)
  y <- matrix(NA, nrow = 2, ncol = n + 1)
  y[, 1] <- y0
  
  with(as.list(parms), {
    for (i in 1:n) {
      sigma <- 0.1
      dB <- rnorm(1, mean = 0, sd = sqrt(dt[i]) * sigma)
      current <- xi_time(time[i], zetas, max(time))
      
      du_dt <- gamma * (y[1, i] - y[1, i]^3 / 3 + y[2, i] + current)
      dv_dt <- -(1 / gamma) * (y[1, i] - theta0 + theta1 * y[2, i]) + dB
      
      y[, i + 1] <- y[, i] + c(du_dt, dv_dt) * dt[i]
    }
    return(t(y))   # columns: time steps (u,v) rows -> time points
  })
}

# --- summary stats: spike rate and spike duration ---------------------------
compute_spikes <- function(parms, uspike, time, y_initial, zetas) {
  solution <- fitzhugh_nagumo(time, y_initial, parms, zetas)
  # solution rows: time points, columns: (u,v)
  u_t <- solution[, 1]
  
  t_spike <- integer(0)
  start_spikes <- integer(0)
  end_spikes <- integer(0)
  in_spike <- FALSE
  
  for (t in 2:length(u_t)) {
    if (!in_spike) {
      if (uspike <= u_t[t] && uspike > u_t[t - 1]) {
        t_spike <- c(t_spike, t)
        start_spikes <- c(start_spikes, t)
        in_spike <- TRUE
      }
    } else {
      if (u_t[t] < uspike) {
        end_spikes <- c(end_spikes, t - 1)
        in_spike <- FALSE
      }
    }
  }
  if (in_spike) end_spikes <- c(end_spikes, length(u_t))
  
  num_spikes <- length(t_spike)
  spike_dur <- if (length(start_spikes) > 0 && length(end_spikes) > 0) mean(end_spikes - start_spikes) else NA
  # spike_rate scaled like in your original code
  return(list(spike_rate = num_spikes / 0.2, spike_duration = spike_dur))
}

# --- problem setup ----------------------------------------------------------
uspike <- 0.5
time <- seq(0.1, 200, length = 1000)   # same as your code
y_initial <- c(u = 0, v = 0)
gamma_val <- 3.0

# Prior sampling helpers
sample_prior_thetas <- function(N) {
  theta0 <- rtruncnorm(N, a = 0, b = 1, mean = 0.4, sd = 0.3)
  theta1 <- rtruncnorm(N, a = 0, b = 1, mean = 0.4, sd = 0.4)
  cbind(theta0, theta1)
}


utility_estimate <- function(zetas,
                             N_sim = 1000,     # prior samples used to construct ABC proposals
                             tol = 0.05,      # abc tolerance for rejection
                             cur_for_zratio = 0.2 # used when mapping adjusted thetas to spike-rate (you used cur=0.2 earlier)
) {
  if (any(is.na(zetas)) || length(zetas) != 3) return(NA_real_)
  # store replicate utilities
  T_max <- max(time)
  
    # 1) draw pseudo-true parameter and simulate observed dataset (y_obs)
    theta_true <- sample_prior_thetas(N_sim)
    parms_true <- list(gamma = gamma_val, theta0 = theta_true[1], theta1 = theta_true[2])
    y_obs_stats <- compute_spikes(parms_true, uspike, time, y_initial, zetas)
    obs_sumstat <- c(y_obs_stats$spike_rate, y_obs_stats$spike_duration)
    if (any(is.na(obs_sumstat))) {
      replicate_utils<- NA
      
    }
    
    # 2) draw N_sim prior proposals and simulate summaries for given design zetas
    prior_thetas <- sample_prior_thetas(N_sim)
    sim_stats <- matrix(NA, nrow = N_sim, ncol = 2)
    for (i in 1:N_sim) {
      p <- list(gamma = gamma_val, theta0 = prior_thetas[i,1], theta1 = prior_thetas[i,2])
      cs <- compute_spikes(p, uspike, time, y_initial, zetas)
      sim_stats[i, ] <- c(cs$spike_rate, cs$spike_duration)
    }
    # remove rows with NA summaries
    valid_idx <- complete.cases(sim_stats)
    if (sum(valid_idx) < max(10, floor(0.1 * N_sim))) {
      replicate_utils<- NA
      
    }
    sim_stats <- sim_stats[valid_idx, , drop = FALSE]
    prior_thetas_valid <- prior_thetas[valid_idx, , drop = FALSE]
    
    # 3) run ABC rejection to get adjusted posterior samples for observed summary
    # abc::abc requires obs summary as vector, param matrix and sumstat matrix
    abc_fit <- tryCatch({
      abc::abc(target = obs_sumstat,
               param = prior_thetas_valid,
               sumstat = sim_stats,
               tol = tol,
               method = "rejection")
    }, error = function(e) NULL)
    
    if (is.null(abc_fit)) {
      replicate_utils <- NA
    }
    
    # Adjusted values (unadj/adj depending on method)
    post_vals <- abc_fit$unadj.values
    if (is.null(post_vals) || nrow(post_vals) == 0) {
      # try unadjusted if adj empty
      post_vals <- abc_fit$unadj.values
    }
    if (is.null(post_vals) || nrow(post_vals) == 0) {
      replicate_utils <- NA
      
    }
    
    # 4) Map posterior parameter draws to spike-rate values at the reference design 'cur_for_zratio'
    #    (This mimics your existing densratio-based utility approach.)
    #    Note: we evaluate both posterior-generated spike-rates and a fresh sample of prior spike-rates
    # Posterior predictive spike rates:
    post_N <- nrow(post_vals)
    z_post <- numeric(post_N)
    for (j in 1:post_N) {
      p <- list(gamma = gamma_val, theta0 = post_vals[j,1], theta1 = post_vals[j,2])
      cs <- compute_spikes(p, uspike, time, y_initial, c(cur_for_zratio, cur_for_zratio, cur_for_zratio))
      z_post[j] <- cs$spike_rate
    }
    # Prior predictive spike rates (same number)
    prior_sample_for_z <- sample_prior_thetas(post_N)
    z_prior <- numeric(post_N)
    for (j in 1:post_N) {
      p <- list(gamma = gamma_val, theta0 = prior_sample_for_z[j,1], theta1 = prior_sample_for_z[j,2])
      cs <- compute_spikes(p, uspike, time, y_initial, c(cur_for_zratio, cur_for_zratio, cur_for_zratio))
      z_prior[j] <- cs$spike_rate
    }
    # remove NA
    valid_z_idx <- complete.cases(z_post, z_prior)
    if (sum(valid_z_idx) < max(10, floor(0.2 * post_N))) {
      replicate_utils <- NA
      
    }
    z_post <- z_post[valid_z_idx]
    z_prior <- z_prior[valid_z_idx]
    
    # 5) densratio estimation and utility calculation
    dens_obj <- tryCatch({
      densratio::densratio(z_post, z_prior, method = "RuLSIF",
                           alpha = 0.5,
                           sigma = 10 ^ seq(-5, 10, length.out = 20),
                           lambda = 10 ^ seq(-5, 10, length.out = 20))
    }, error = function(e) NULL)
    
    if (is.null(dens_obj)) {
      replicate_utils <- NA
      
    }
    # compute density ratios at prior points, avoid zeros
    ratios <- dens_obj$compute_density_ratio(z_prior)
    ratios[is.na(ratios) | ratios <= 0] <- NA
    if (all(is.na(ratios))) {
      replicate_utils <- NA
      
    }
    mean_log_ratio <- mean(log(ratios), na.rm = TRUE)
    replicate_utils <- mean_log_ratio
 
  
  # return mean over replicates (or NA if all NA)
  if (all(is.na(replicate_utils))) return(NA_real_)
  return(mean(replicate_utils, na.rm = TRUE))
}

# --- Bayesian Optimization wrapper for rBayesianOptimization ----------------
# rBayesianOptimization expects an R function with named numeric arguments and returns list(Score= , Pred=)
BO_objective <- function(zeta1, zeta2, zeta3) {
  zetas <- c(zeta1, zeta2, zeta3)
  # tune these counts for speed vs accuracy
  N_sim_per_eval <- 1000    # number prior particles per evaluation (reduce to speed up)
  tol <- 0.06

  set.seed(123)            # for reproducibility of stochastic sims within function
  util <- utility_estimate(zetas, N_sim = N_sim_per_eval, tol = tol,
                           cur_for_zratio = 0.2)
  # rBayesianOptimization maximizes Score. If util is NA, return a low value.
  if (is.na(util) || !is.finite(util)) util <- -1e6
  return(list(Score = util, Pred = util))
}

# --- Run Bayesian Optimization ----------------------------------------------
# Bounds for each zeta: [0, 0.8] as you specified earlier
bounds <- list(zeta1 = c(0, 0.8),
               zeta2 = c(0, 0.8),
               zeta3 = c(0, 0.8))

# Initial grid (random initial points)
init_points <- 50
n_iter <- 25   # total additional BO iterations (increase for more thorough search)
set.seed(42)
BO_res <- BayesianOptimization(
  FUN = BO_objective,
  bounds = bounds,
  init_points = init_points,
  n_iter = n_iter,
  acq = "ucb",            # acquisition: "ucb", "ei", "poi"
  kappa = 2.576,          # controls exploration for UCB
  eps = 0.0,
  verbose = TRUE
)
end<-Sys.time()
print(end-start)
cat("BayesOpt finished. Best found:\n")
print(BO_res$Best_Par)
print(BO_res$Best_Value)
#zeta1 = 0.1661272	zeta2 = 0.07198441	zeta3 = 0.09097489	Value = 0.4802889 
# Save results
save(BO_res, file = "fhn_bo_optionA_results.RData")

# --- stop cluster -----------------------------------------------------------
stopCluster(cl)
registerDoSEQ()
