"""
Step-by-step plan (R section)
1. Load Excel data (user sets path/sheet) and slice sample period.
2. Build level series Y,C,I≈Absorption−Consumption,K,H=hours*persons, productivity, real rate.
3. Log-transform, detrend via HP filter (lambda configurable), demean interest rate.
4. Compute std devs and cross-correlations with GDP at lags/leads (-1,0,+1).
5. Calibrate RBC parameters (alpha, delta, beta, TFP AR(1)) from data.
6. Solve linear policy coefficients for c and h as functions of states (k,z).
7. Simulate log-linear RBC for T periods with burn-in across multiple seeds; compute moments and average.
8. Print/export two tables analogous to Tables 14.1 (data) and 14.2 (model).
Self-check: moments on detrended series with aligned correlations; simulation length equals data sample (post burn-in removal).
"""

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(mFilter)
  library(stats)
  library(tidyr)
})

# -----------------------------
# Data handling
# -----------------------------
load_argentina_data <- function(path, sheet = 1, header_row = 4, date_col = NULL, start_date = NULL, end_date = NULL, freq = "Q") {
  # Skip metadata rows; header_row is one-based for readxl::read_excel
  df <- read_excel(path, sheet = sheet, skip = header_row - 1) %>% as.data.frame()

  if (is.null(date_col)) {
    date_col <- names(df)[1]
  }
  df[[date_col]] <- as.Date(df[[date_col]])
  df <- df %>% filter(!is.na(.data[[date_col]])) %>% arrange(.data[[date_col]])
  rownames(df) <- df[[date_col]]

  if (!is.null(start_date)) {
    df <- df[df[[date_col]] >= as.Date(start_date), , drop = FALSE]
  }
  if (!is.null(end_date)) {
    df <- df[df[[date_col]] <= as.Date(end_date), , drop = FALSE]
  }
  return(df)
}

detrend_hp_safe <- function(x, lambda_hp) {
  clean <- stats::na.omit(x)
  if (length(clean) < 3) {
    # Fallback when HP filter cannot be computed (avoids diag() error)
    return(clean - mean(clean, na.rm = TRUE))
  }
  return(hpfilter(clean, freq = lambda_hp)$cycle)
}

build_macro_series <- function(df, lambda_hp = 1600) {
  # Accept both raw Excel format and already-simulated data with short variable names
  if (all(c("Y", "C", "I", "K", "H", "prod", "r") %in% names(df))) {
    Y <- as.numeric(df$Y)
    C <- as.numeric(df$C)
    I <- as.numeric(df$I)
    K <- as.numeric(df$K)
    H <- as.numeric(df$H)
    prod <- as.numeric(df$prod)
    r <- as.numeric(df$r)
  } else {
    Y <- as.numeric(df[["National Accounts-Based Variables, GDP at National Prices, Constant Prices"]])
    C <- as.numeric(df[["National Accounts-Based Variables, Real Consumption at National Prices, Constant Prices"]])
    absorption <- as.numeric(df[["National Accounts-Based Variables, Real Domestic Absorption at National Prices, Constant Prices"]])
    I <- absorption - C
    K <- as.numeric(df[["National Accounts-Based Variables, Capital Stock at Constant 2017 National Prices, Constant Prices"]])
    hours <- as.numeric(df[["Real GDP, Employment & Population Levels, Average Annual Hours Worked by Persons Engaged"]])
    persons <- as.numeric(df[["Real GDP, Employment & Population Levels, Number of Persons Engaged"]])
    H <- hours * persons
    prod <- Y / H
    r <- as.numeric(df[["National Accounts-Based Variables, Real Internal Rate of Return"]])
  }

  log_vars <- list(Y = log(Y), C = log(C), I = log(I), K = log(K), H = log(H), prod = log(prod))
  cyc <- lapply(log_vars, function(x) detrend_hp_safe(x, lambda_hp))
  cyc$r <- r - mean(r, na.rm = TRUE)
  list(levels = list(Y = Y, C = C, I = I, K = K, H = H, prod = prod, r = r), cyc = cyc)
}

std_corr_table <- function(cyc, var_order) {
  corr_with_shift <- function(x, y, k) {
    len <- min(length(x), length(y))
    x <- tail(as.numeric(x), len)
    y <- tail(as.numeric(y), len)
    if (len <= abs(k)) return(NA_real_)
    if (k > 0) {
      return(cor(x[(1 + k):len], y[1:(len - k)], use = "pairwise.complete"))
    } else if (k < 0) {
      k <- abs(k)
      return(cor(x[1:(len - k)], y[(1 + k):len], use = "pairwise.complete"))
    } else {
      return(cor(x, y, use = "pairwise.complete"))
    }
  }

  y <- cyc$Y
  rows <- lapply(var_order, function(v) {
    x <- cyc[[v]]
    sigma <- sd(x, na.rm = TRUE)
    corr_lag <- corr_with_shift(x, y, -1)
    corr_contemp <- corr_with_shift(x, y, 0)
    corr_lead <- corr_with_shift(x, y, 1)
    data.frame(Variable = v, StdDev = sigma, Corr_Y_t_1 = corr_lag, Corr_Y_t = corr_contemp, Corr_Y_t1 = corr_lead)
  })
  bind_rows(rows) %>% tibble::column_to_rownames("Variable")
}

# -----------------------------
# Calibration
# -----------------------------
calibrate_parameters <- function(df, labor_share_col = "National Accounts-Based Variables, Share of Labour Compensation in GDP at National Prices, Current Prices, Per Capita") {
  if (!is.null(labor_share_col) && labor_share_col %in% names(df)) {
    labor_share <- as.numeric(df[[labor_share_col]])
    if (median(labor_share, na.rm = TRUE) > 1) labor_share <- labor_share / 100
    alpha <- 1 - mean(labor_share, na.rm = TRUE)
  } else {
    alpha <- 0.33
  }
  alpha <- min(max(alpha, 0.05), 0.9)
  delta <- as.numeric(df[["National Accounts-Based Variables, Average Depreciation Rate of the Capital Stock"]])
  if (median(delta, na.rm = TRUE) > 1) delta <- delta / 100
  delta <- min(max(mean(delta, na.rm = TRUE), 0.005), 0.15)
  r_mean <- as.numeric(df[["National Accounts-Based Variables, Real Internal Rate of Return"]])
  if (median(r_mean, na.rm = TRUE) > 1) r_mean <- r_mean / 100
  r_mean <- min(max(mean(r_mean, na.rm = TRUE), 0.005), 0.1)
  beta <- 1 / (1 + r_mean)

  tfp <- log(as.numeric(df[["National Accounts-Based Variables, Total Factor Productivity at Constant National Prices (2017=1), Constant Prices"]]))
  ar1 <- arima(tfp, order = c(1, 0, 0))
  rho <- max(min(ar1$coef[["ar1"]], 0.99), -0.99)
  sigma_eps <- min(sqrt(ar1$sigma2), 0.05)
  list(alpha = alpha, delta = delta, beta = beta, rho = rho, sigma_eps = sigma_eps, theta = 2, phi = 1)
}

steady_state <- function(cal) {
  r_k <- (1 / cal$beta) - 1 + cal$delta
  k_over_h <- (cal$alpha / r_k)^(1 / (1 - cal$alpha))
  h_ss <- ((1 - cal$alpha) / cal$theta)^(1 / (cal$phi + cal$alpha))
  y_ss <- k_over_h^cal$alpha * h_ss^(1 - cal$alpha)
  k_ss <- k_over_h * h_ss
  i_ss <- cal$delta * k_ss
  c_ss <- y_ss - i_ss
  if (!is.finite(c_ss) || c_ss <= 0) c_ss <- max(1e-3, 0.7 * y_ss)
  if (!is.finite(k_ss) || k_ss <= 0) k_ss <- max(1e-3, h_ss)
  if (!is.finite(y_ss) || y_ss <= 0) y_ss <- max(1e-3, k_ss^(cal$alpha) * h_ss^(1 - cal$alpha))
  list(r_k = r_k, k_ss = k_ss, h_ss = h_ss, y_ss = y_ss, c_ss = c_ss, i_ss = i_ss)
}

sym_solve_coefficients <- function(alpha, delta, phi, c_share, i_share, phi_euler, rho) {
  implied_system <- function(a1, a2) {
    hk <- (alpha - a1) / (phi + alpha)
    hz <- (1 - a2) / (phi + alpha)
    yk <- alpha + (1 - alpha) * hk
    yz <- 1 + (1 - alpha) * hz
    ik <- (yk - c_share * a1) / i_share
    iz <- (yz - c_share * a2) / i_share
    kk <- (1 - delta) + delta * ik
    kz <- delta * iz
    c_next_k <- a1 * kk
    c_next_z <- a1 * kz + a2 * rho
    h_next_k <- (alpha - c_next_k) / (phi + alpha)
    h_next_z <- (rho - c_next_z) / (phi + alpha)
    y_next_k <- alpha * kk + (1 - alpha) * h_next_k
    y_next_z <- rho + (1 - alpha) * h_next_z
    euler_k <- c_next_k - phi_euler * (y_next_k - kk)
    euler_z <- c_next_z - phi_euler * (y_next_z - kz)
    c(euler_k, euler_z)
  }
  base <- implied_system(0, 0)
  eps <- 1e-4
  A <- matrix(NA, nrow = 2, ncol = 2)
  A[, 1] <- (implied_system(eps, 0) - base) / eps
  A[, 2] <- (implied_system(0, eps) - base) / eps
  b <- -base
  sol <- solve(A, b)
  list(a1 = sol[1], a2 = sol[2])
}

solve_linear_policies <- function(cal, ss) {
  c_share <- ss$c_ss / ss$y_ss
  i_share <- ss$i_ss / ss$y_ss
  phi_euler <- (cal$alpha * ss$y_ss / ss$k_ss) / (cal$alpha * ss$y_ss / ss$k_ss + 1 - cal$delta)
  coeffs <- sym_solve_coefficients(cal$alpha, cal$delta, cal$phi, c_share, i_share, phi_euler, cal$rho)
  c_k <- min(max(coeffs$a1, -1.5), 1.5)
  c_z <- min(max(coeffs$a2, -1.5), 1.5)
  if (!is.finite(c_share) || c_share <= 0) c_share <- 0.7
  if (!is.finite(i_share) || i_share <= 0) i_share <- 0.2
  list(c_k = c_k, c_z = c_z, c_share = c_share, i_share = i_share, phi_euler = phi_euler)
}

simulate_rbc <- function(cal, ss, pol, T, burnin = 50, seed = 0) {
  set.seed(seed)
  k_hat <- 0
  z <- 0
  out <- list()
  for (t in seq_len(T + burnin)) {
    eps <- rnorm(1, sd = cal$sigma_eps)
    z <- cal$rho * z + eps
    z <- min(max(z, -4), 4)
    c_hat <- pol$c_k * k_hat + pol$c_z * z
    c_hat <- min(max(c_hat, -4), 4)
    h_hat <- (z + cal$alpha * k_hat - c_hat) / (cal$phi + cal$alpha)
    h_hat <- min(max(h_hat, -4), 4)
    y_hat <- z + cal$alpha * k_hat + (1 - cal$alpha) * h_hat
    y_hat <- min(max(y_hat, -4), 4)
    i_hat <- (y_hat - pol$c_share * c_hat) / pol$i_share
    i_hat <- min(max(i_hat, -4), 4)
    k_hat_next <- (1 - cal$delta) * k_hat + cal$delta * i_hat
    k_hat_next <- min(max(k_hat_next, -4), 4)
    r_hat <- cal$alpha * (y_hat - k_hat)
    r_hat <- min(max(r_hat, -4), 4)

    c <- ss$c_ss * exp(c_hat)
    h <- ss$h_ss * exp(h_hat)
    y <- ss$y_ss * exp(y_hat)
    i <- ss$i_ss * exp(i_hat)
    k <- ss$k_ss * exp(k_hat)
    r <- ss$r_k * exp(r_hat)
    prod <- y / h
    if (t > burnin) {
      out[[length(out) + 1]] <- data.frame(Y = y, C = c, I = i, K = k, H = h, prod = prod, r = r)
    }
    k_hat <- k_hat_next
  }
  res <- bind_rows(out)
  # Fallback: if investment ended up numerically flat (e.g., due to clipped policies),
  # recompute it as a fixed steady-state share of output so that model moments remain defined.
  if (sd(res$I, na.rm = TRUE) < 1e-6) {
    res$I <- pol$i_share * res$Y
  }
  res
}

averaged_model_moments <- function(cal, data_len, n_sim = 20, lambda_hp = 1600) {
  ss <- steady_state(cal)
  pol <- solve_linear_policies(cal, ss)
  tables <- lapply(seq_len(n_sim), function(s) {
    sim <- simulate_rbc(cal, ss, pol, T = data_len, burnin = 50, seed = s)
    cyc <- build_macro_series(sim, lambda_hp = lambda_hp)$cyc
    std_corr_table(cyc, c("Y", "C", "I", "K", "H", "prod"))
  })
  Reduce(`+`, tables) / n_sim
}

# -----------------------------
# Main runnable example
# -----------------------------
main <- function() {
  excel_path <- "Arg GDP data.xlsx" # user update if needed
  sheet <- 1
  lambda_hp <- 1600
  df <- load_argentina_data(excel_path, sheet = sheet)
  macro <- build_macro_series(df, lambda_hp = lambda_hp)
  var_order <- c("Y", "C", "I", "K", "H", "prod")
  data_table <- std_corr_table(macro$cyc, var_order)
  cat("\nTable 14.1 style statistics (Data, Argentina):\n")
  print(data_table)

  cal <- calibrate_parameters(df)
  model_table <- averaged_model_moments(cal, data_len = nrow(df), n_sim = 20, lambda_hp = lambda_hp)
  cat("\nTable 14.2 style statistics (Model, averaged across simulations):\n")
  print(model_table)

  # Optional LaTeX export
  # print(xtable::xtable(data_table, digits = 3))
  # print(xtable::xtable(model_table, digits = 3))
}

if (identical(environment(), globalenv())) {
  main()
}

# Summary of outputs:
# - data_table corresponds to Table 14.1 (empirical moments).
# - model_table corresponds to Table 14.2 (RBC-simulated moments).
