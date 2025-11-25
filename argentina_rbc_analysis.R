############################################################
# RBC para Argentina: Tablas 14.1 (Datos) y 14.2 (Modelo)
############################################################

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(mFilter)
  library(stats)
  library(tidyr)
})

############################################################
# 1. Carga de datos y series macro
############################################################

load_argentina_data <- function(path, sheet = 1, header_row = 4,
                                date_col = NULL,
                                start_date = NULL,
                                end_date   = NULL,
                                freq = "Q") {
  df <- read_excel(path, sheet = sheet, skip = header_row - 1) %>%
    as.data.frame()
  
  if (is.null(date_col)) {
    date_col <- names(df)[1]
  }
  df[[date_col]] <- as.Date(df[[date_col]])
  df <- df %>%
    filter(!is.na(.data[[date_col]])) %>%
    arrange(.data[[date_col]])
  rownames(df) <- df[[date_col]]
  
  if (!is.null(start_date)) {
    df <- df[df[[date_col]] >= as.Date(start_date), , drop = FALSE]
  }
  if (!is.null(end_date)) {
    df <- df[df[[date_col]] <= as.Date(end_date), , drop = FALSE]
  }
  df
}

detrend_hp_safe <- function(x, lambda_hp) {
  clean <- stats::na.omit(x)
  if (length(clean) < 3L) {
    return(clean - mean(clean, na.rm = TRUE))
  }
  hpfilter(clean, freq = lambda_hp)$cycle
}

build_macro_series <- function(df, lambda_hp = 1600) {
  # Caso 1: datos simulados
  if (all(c("Y","C","I","K","H","prod","r") %in% names(df))) {
    Y    <- as.numeric(df$Y)
    C    <- as.numeric(df$C)
    I    <- as.numeric(df$I)
    K    <- as.numeric(df$K)
    H    <- as.numeric(df$H)
    prod <- as.numeric(df$prod)
    r    <- as.numeric(df$r)
  } else {
    # Caso 2: Excel PWT/NA
    Y <- as.numeric(df[["National Accounts-Based Variables, GDP at National Prices, Constant Prices"]])
    C <- as.numeric(df[["National Accounts-Based Variables, Real Consumption at National Prices, Constant Prices"]])
    absorption <- as.numeric(df[["National Accounts-Based Variables, Real Domestic Absorption at National Prices, Constant Prices"]])
    I <- absorption - C
    K <- as.numeric(df[["National Accounts-Based Variables, Capital Stock at Constant 2017 National Prices, Constant Prices"]])
    hours   <- as.numeric(df[["Real GDP, Employment & Population Levels, Average Annual Hours Worked by Persons Engaged"]])
    persons <- as.numeric(df[["Real GDP, Employment & Population Levels, Number of Persons Engaged"]])
    H    <- hours * persons
    prod <- Y / H
    r    <- as.numeric(df[["National Accounts-Based Variables, Real Internal Rate of Return"]])
  }
  
  log_vars <- list(
    Y    = log(Y),
    C    = log(C),
    I    = log(I),
    K    = log(K),
    H    = log(H),
    prod = log(prod)
  )
  cyc <- lapply(log_vars, function(x) detrend_hp_safe(x, lambda_hp))
  cyc$r <- r - mean(r, na.rm = TRUE)
  
  list(
    levels = list(Y = Y, C = C, I = I, K = K, H = H, prod = prod, r = r),
    cyc    = cyc
  )
}

############################################################
# 2. Momentos: std y correlaciones cruzadas con Y
############################################################

safe_cor <- function(x, y) {
  ok <- complete.cases(x, y)
  x2 <- x[ok]; y2 <- y[ok]
  if (length(x2) < 2L) return(NA_real_)
  if (sd(x2) < 1e-10 || sd(y2) < 1e-10) return(0)  # evita NA si una serie es casi plana
  cor(x2, y2)
}

std_corr_table <- function(cyc, var_order) {
  corr_with_shift <- function(x, y, k) {
    len <- min(length(x), length(y))
    x <- tail(as.numeric(x), len)
    y <- tail(as.numeric(y), len)
    if (len <= abs(k)) return(NA_real_)
    if (k > 0) {
      # corr(x_t, y_{t-1})
      safe_cor(x[(1 + k):len], y[1:(len - k)])
    } else if (k < 0) {
      k <- abs(k)
      # corr(x_t, y_{t+1})
      safe_cor(x[1:(len - k)], y[(1 + k):len])
    } else {
      # contemporánea
      safe_cor(x, y)
    }
  }
  
  y <- cyc$Y
  rows <- lapply(var_order, function(v) {
    x <- cyc[[v]]
    sigma <- sd(x, na.rm = TRUE)
    corr_lag      <- corr_with_shift(x, y, -1)
    corr_contemp  <- corr_with_shift(x, y,  0)
    corr_lead     <- corr_with_shift(x, y,  1)
    data.frame(
      Variable   = v,
      StdDev     = sigma,
      Corr_Y_t_1 = corr_lag,
      Corr_Y_t   = corr_contemp,
      Corr_Y_t1  = corr_lead
    )
  })
  bind_rows(rows) %>% tibble::column_to_rownames("Variable")
}

############################################################
# 3. Calibración RBC (desde datos, pero con sigma_eps > 0)
############################################################

calibrate_parameters <- function(
    df,
    labor_share_col = "National Accounts-Based Variables, Share of Labour Compensation in GDP at National Prices, Current Prices, Per Capita"
) {
  # alpha
  if (!is.null(labor_share_col) && labor_share_col %in% names(df)) {
    labor_share <- as.numeric(df[[labor_share_col]])
    if (median(labor_share, na.rm = TRUE) > 1) labor_share <- labor_share / 100
    alpha <- 1 - mean(labor_share, na.rm = TRUE)
  } else {
    alpha <- 0.33
  }
  
  # delta
  delta <- as.numeric(df[["National Accounts-Based Variables, Average Depreciation Rate of the Capital Stock"]])
  if (median(delta, na.rm = TRUE) > 1) delta <- delta / 100
  delta <- mean(delta, na.rm = TRUE)
  
  # r medio y beta
  r_series <- as.numeric(df[["National Accounts-Based Variables, Real Internal Rate of Return"]])
  if (median(r_series, na.rm = TRUE) > 1) r_series <- r_series / 100
  r_mean <- mean(r_series, na.rm = TRUE)
  beta <- 1 / (1 + r_mean)
  
  # TFP AR(1) -> rho y sigma_eps
  tfp_raw <- as.numeric(
    df[[
      "National Accounts-Based Variables, Total Factor Productivity at Constant National Prices (2017=1), Constant Prices"
    ]]
  )
  tfp <- log(tfp_raw)
  tfp <- tfp[is.finite(tfp)]
  
  rho <- 0.95
  sigma_eps <- 0.01
  
  if (length(tfp) > 5L && sd(tfp) > 1e-6) {
    fit <- try(arima(tfp, order = c(1,0,0)), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      if (!is.null(fit$coef[["ar1"]])) {
        rho <- max(min(fit$coef[["ar1"]], 0.99), -0.99)
      }
      sig <- sqrt(fit$sigma2)
      if (is.finite(sig)) {
        sigma_eps <- min(max(sig, 0.001), 0.02)  # entre 0.1% y 2% aprox
      }
    }
  }
  
  list(
    alpha = alpha,
    delta = delta,
    beta  = beta,
    rho   = rho,
    sigma_eps = sigma_eps,
    theta = 2,   # escala desutilidad trabajo
    phi   = 1    # parámetro tipo Frisch
  )
}

steady_state <- function(cal) {
  alpha <- cal$alpha
  delta <- cal$delta
  beta  <- cal$beta
  theta <- cal$theta
  phi   <- cal$phi
  
  r_k     <- (1 / beta) - 1 + delta
  k_over_h <- (alpha / r_k)^(1 / (1 - alpha))
  h_ss    <- ((1 - alpha) / theta)^(1 / (phi + alpha))
  y_ss    <- k_over_h^alpha * h_ss^(1 - alpha)
  k_ss    <- k_over_h * h_ss
  i_ss    <- delta * k_ss
  c_ss    <- y_ss - i_ss
  
  list(
    r_k  = r_k,
    k_ss = k_ss,
    h_ss = h_ss,
    y_ss = y_ss,
    c_ss = c_ss,
    i_ss = i_ss
  )
}

############################################################
# 4. Políticas lineales: c_hat = a1 k_hat + a2 z
############################################################

sym_solve_coefficients <- function(alpha, delta, phi,
                                   c_share, i_share,
                                   phi_euler, rho) {
  implied_system <- function(a1, a2) {
    # h_t(k,z)
    hk <- (alpha - a1) / (phi + alpha)
    hz <- (1     - a2) / (phi + alpha)
    
    # y_t = z + alpha k + (1-alpha) h
    yk <- alpha + (1 - alpha) * hk
    yz <- 1     + (1 - alpha) * hz
    
    # i_t de la RC: y = c_share c + i_share i
    ik <- (yk - c_share * a1) / i_share
    iz <- (yz - c_share * a2) / i_share
    
    # k_{t+1}
    kk <- (1 - delta) + delta * ik
    kz <- delta * iz
    
    # consumo y trabajo en t+1
    c_next_k <- a1 * kk
    c_next_z <- a1 * kz + a2 * rho
    
    h_next_k <- (alpha - c_next_k) / (phi + alpha)
    h_next_z <- (rho  - c_next_z) / (phi + alpha)
    
    y_next_k <- alpha * kk + (1 - alpha) * h_next_k
    y_next_z <- rho + (1 - alpha) * h_next_z
    
    # Euler linealizada (misma forma que en tu Python)
    euler_k <- c_next_k - phi_euler * (y_next_k - kk)
    euler_z <- c_next_z - phi_euler * (y_next_z - kz)
    c(euler_k, euler_z)
  }
  
  base <- implied_system(0, 0)
  eps  <- 1e-4
  A <- matrix(NA_real_, 2, 2)
  A[, 1] <- (implied_system(eps, 0) - base) / eps
  A[, 2] <- (implied_system(0, eps) - base) / eps
  b <- -base
  sol <- solve(A, b)
  list(a1 = as.numeric(sol[1]), a2 = as.numeric(sol[2]))
}

solve_linear_policies <- function(cal, ss) {
  c_share <- ss$c_ss / ss$y_ss
  i_share <- ss$i_ss / ss$y_ss
  phi_euler <- (cal$alpha * ss$y_ss / ss$k_ss) /
    (cal$alpha * ss$y_ss / ss$k_ss + 1 - cal$delta)
  
  coeffs <- sym_solve_coefficients(
    cal$alpha, cal$delta, cal$phi,
    c_share, i_share,
    phi_euler, cal$rho
  )
  
  list(
    c_k       = coeffs$a1,
    c_z       = coeffs$a2,
    c_share   = c_share,
    i_share   = i_share,
    phi_euler = phi_euler
  )
}

############################################################
# 5. Simulación RBC LOG-LINEAL (igual al Python)
############################################################

simulate_rbc <- function(cal, ss, pol,
                         T_sim, burnin = 5, seed = 0) {
  set.seed(seed)
  
  alpha <- cal$alpha
  delta <- cal$delta
  rho   <- cal$rho
  phi   <- cal$phi
  
  c_share <- pol$c_share
  i_share <- pol$i_share
  a1 <- pol$c_k
  a2 <- pol$c_z
  
  k_hat <- 0
  z     <- 0
  
  out <- vector("list", T_sim)
  t_out <- 0L
  
  for (t in seq_len(T_sim + burnin)) {
    eps <- rnorm(1, sd = cal$sigma_eps)
    z   <- rho * z + eps
    
    c_hat <- a1 * k_hat + a2 * z
    c_hat <- max(min(c_hat,  0.5), -0.5)
    
    h_hat <- (z + alpha * k_hat - c_hat) / (phi + alpha)
    h_hat <- max(min(h_hat,  0.5), -0.5)
    
    y_hat <- z + alpha * k_hat + (1 - alpha) * h_hat
    y_hat <- max(min(y_hat,  0.5), -0.5)
    
    i_hat <- (y_hat - c_share * c_hat) / i_share
    i_hat <- max(min(i_hat,  0.5), -0.5)
    
    k_hat_next <- (1 - delta) * k_hat + delta * i_hat
    r_hat <- alpha * (y_hat - k_hat)
    
    # niveles (no imponen Y=C+I exactamente, como en la log-lineal)
    C <- ss$c_ss * exp(c_hat)
    H <- ss$h_ss * exp(h_hat)
    Y <- ss$y_ss * exp(y_hat)
    I <- ss$i_ss * exp(i_hat)
    K <- ss$k_ss * exp(k_hat)
    R <- ss$r_k  * exp(r_hat)
    prod <- Y / H
    
    if (t > burnin) {
      t_out <- t_out + 1L
      out[[t_out]] <- data.frame(
        Y    = Y,
        C    = C,
        I    = I,
        K    = K,
        H    = H,
        prod = prod,
        r    = R
      )
    }
    
    k_hat <- k_hat_next
  }
  
  dplyr::bind_rows(out[seq_len(t_out)])
}

############################################################
# 6. Momentos del modelo promediados
############################################################

averaged_model_moments <- function(cal, data_len,
                                   n_sim = 20,
                                   lambda_hp = 1600) {
  ss  <- steady_state(cal)
  pol <- solve_linear_policies(cal, ss)
  
  tables <- lapply(seq_len(n_sim), function(s) {
    sim <- simulate_rbc(cal, ss, pol,
                        T_sim = data_len,
                        burnin = 5,
                        seed = s - 1L)
    cyc <- build_macro_series(sim, lambda_hp = lambda_hp)$cyc
    std_corr_table(cyc, c("Y", "C", "I", "K", "H", "prod"))
  })
  
  # Promedio elemento a elemento
  mats <- lapply(tables, as.matrix)
  arr  <- simplify2array(mats)
  avg_mat <- apply(arr, c(1, 2), mean, na.rm = TRUE)
  
  out <- tables[[1]]
  out[,] <- avg_mat
  out
}

############################################################
# 7. MAIN
############################################################

main <- function() {
  excel_path <- "Arg GDP data.xlsx"  # ajustá path
  sheet      <- 1
  lambda_hp  <- 1600
  
  df     <- load_argentina_data(excel_path, sheet = sheet)
  macro  <- build_macro_series(df, lambda_hp = lambda_hp)
  var_order  <- c("Y", "C", "I", "K", "H", "prod")
  
  data_table <- std_corr_table(macro$cyc, var_order)
  cat("\nTable 14.1 style statistics (Data, Argentina):\n")
  print(data_table)
  
  cal <- calibrate_parameters(df)
  cat("\nCalibration used (RBC):\n")
  print(cal)
  
  cal$alpha     <- 1/3
  cal$delta     <- 0.10
  cal$beta      <- 1/(1 + 0.01)
  cal$rho       <- 0.95   # VARPHI del libro
  cal$sigma_eps <- 0.01
  cal$theta     <- 2
  cal$phi       <- 1
  
  
  model_table <- averaged_model_moments(
    cal,
    data_len = nrow(df),
    n_sim    = 20,
    lambda_hp = lambda_hp
  )
  cat("\nTable 14.2 style statistics (Model, averaged across simulations):\n")
  print(model_table)
}

if (identical(environment(), globalenv())) {
  main()
}
