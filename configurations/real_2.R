# real_2: based on the original, archived 175_3 ACTG175 implementation.
# Provenance: corresponding working Code/real_2.R; latest real RData has no embedded source.
# QR2 uses full-training nuisance predictions, not inner OOF predictions.
# Outer cross-fitting, tuning CV, and q normalization/truncation are retained.
# QR1/QR2/TR pi and O-TR p probabilities are clipped to [0.01, 0.99].
# TR q restored to 175_2 KNN density and normalization (2026-09-19).
# MQR1 and QR1 are enabled with their original implementations.
# Internal MTR2/O-TR2 code retains its original switches and RNG order.
# Report MQR1, MQR2, O-TR', QR1, QR2, TR, and O-TR with named output columns.
# Run from the folder containing this file and real_2_data_generate.R.
rm(list=ls())
require(glmnet);require(itertools);require(foreach);require(doParallel);require(RCAL);require(truncnorm);require(mvtnorm);require(corpcor);require(MASS);require(matrixStats);require(ranger);require(caret);require(cobalt);require(densratio);require(butcher);require(randomForestSRC);require(logger);require(filelock);require(sn);require(doRNG)
library(mice)
library(survival)
library(dplyr);library(tidyr);library(purrr)
library(speff2trial); library(ks)
#--------------------functions-------------------------
source("real_2_data_generate.R")

# Nonparametric implementations ported from jasa_sparse6.R.
lasso_var_selection <- function(X, y, task = 0) {
  # std_result <- standardize_data(X)
  # X_scaled <- std_result$train
  # if (task == 1) {  # 分类任务
  #   fit <- cv.glmnet(x = X_scaled, y = y, family = "binomial", alpha = 1)
  # } else {          # 回归任务
  #   fit <- cv.glmnet(x = X_scaled, y = y, family = "gaussian", alpha = 1)
  # }
  # beta <- coef(fit, s = "lambda.min")
  # selected_vars <- rownames(beta)[which(abs(beta)>1e-6)][-1]
  # if(length(selected_vars)==0){
  # rf <- ranger(y ~ ., data = data.frame(X=X, y=y), importance = "impurity", num.trees=50)
  # importance<-importance(rf)
  # selected_vars <- names(importance)[importance > mean(importance)- sd(importance)]#
  # # }
  # if(length(selected_vars)==0){selected_vars=colnames(X)}
  selected_vars=colnames(X)  
  selected_vars <- gsub("^X\\.", "", selected_vars)  # 去除前缀
  selected_vars <- gsub("^M\\.", "", selected_vars)
  return(list(selected_vars = selected_vars))#, model = fit
}
select_different_M_components_ks <- function(
    X, M, A,
    max_selected = min(5L, ncol(M)),
    min_selected = 1L,
    rel_score_min = 0.10,
    ntree = 200,
    nodesize = 5,
    mtry = NULL,
    seed = 123) {
  
  X <- as.matrix(X)
  M <- as.matrix(M)
  A <- as.numeric(A)
  
  n <- nrow(X)
  p <- ncol(X)
  d <- ncol(M)
  
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(p))
  if (is.null(colnames(M))) colnames(M) <- paste0("M", seq_len(d))
  
  idx1 <- which(A == 1)
  idx0 <- which(A == 0)
  
  if (length(idx1) < 20 || length(idx0) < 20) {
    stop("Too few observations in one treatment arm for mediator screening.")
  }
  
  max_selected <- min(max_selected, d)
  min_selected <- min(min_selected, max_selected)
  
  if (is.null(mtry)) {
    mtry <- max(1L, round(sqrt(p)))
  }
  mtry <- min(mtry, p)
  
  formula_multi <- as.formula(
    paste0(
      "cbind(", paste(colnames(M), collapse = ","), ") ~ ",
      paste(colnames(X), collapse = "+")
    )
  )
  
  fit_arm <- function(idx, arm_seed) {
    
    dat <- data.frame(
      X[idx, , drop = FALSE],
      M[idx, , drop = FALSE],
      check.names = FALSE
    )
    
    fit <- rfsrc(
      formula_multi,
      data = dat,
      ntree = ntree,
      block.size = d,
      importance = FALSE,
      forest = TRUE,
      splitrule = "mv",
      nodesize = min(nodesize, length(idx) - 1L),
      mtry = mtry,
      bootstrap = "by.root",
      seed = arm_seed
    )
    
    oob_pred <- get.mv.predicted(fit, oob = TRUE)
    if (d == 1) oob_pred <- matrix(oob_pred, ncol = 1)
    
    pred <- predict(
      fit,
      newdata = data.frame(X, check.names = FALSE)
    )
    
    pred_all <- sapply(seq_len(d), function(j) {
      pred$regrOutput[[j]][[1]]
    })
    if (d == 1) pred_all <- matrix(pred_all, ncol = 1)
    
    list(
      pred_all = pred_all,
      oob_pred = oob_pred
    )
  }
  
  fit1 <- fit_arm(idx1, seed + 1)
  fit0 <- fit_arm(idx0, seed + 2)
  
  # Conditional-mean discrepancy evaluated on the same empirical X distribution
  mean_diff2 <- colMeans(
    (fit1$pred_all - fit0$pred_all)^2,
    na.rm = TRUE
  )
  
  # OOB residual variances, avoiding in-sample residual shrinkage
  resid1 <- M[idx1, , drop = FALSE] - fit1$oob_pred
  resid0 <- M[idx0, , drop = FALSE] - fit0$oob_pred
  
  var1 <- apply(resid1, 2, var, na.rm = TRUE)
  var0 <- apply(resid0, 2, var, na.rm = TRUE)
  
  positive_var <- c(var0, var1)
  positive_var <- positive_var[is.finite(positive_var) & positive_var > 0]
  
  if (length(positive_var) == 0) {
    stop("Unable to estimate residual variances for mediator screening.")
  }
  
  var_floor <- max(1e-8, 1e-6 * median(positive_var))
  
  var1[!is.finite(var1) | var1 < var_floor] <- var_floor
  var0[!is.finite(var0) | var0 < var_floor] <- var_floor
  
  # Symmetric Gaussian KL / Jeffreys discrepancy:
  # combines conditional-mean and marginal conditional-variance changes
  score <- 0.5 * (
    var1 / var0 +
      var0 / var1 - 2 +
      mean_diff2 * (1 / var1 + 1 / var0)
  )
  
  score[!is.finite(score)] <- 0
  
  ranked <- order(score, decreasing = TRUE)
  ranked <- order(score, decreasing = TRUE)
  
  finite_scores <- score[is.finite(score)]
  if (d < 3L) {
    # ACTG175 has only two mediators: a bulk noise threshold is not
    # estimable. Retain all components allowed by max_selected.
    return(ranked[seq_len(min(max_selected, d))])
  }
  if (length(finite_scores) < 3L) {
    stop("Too few finite mediator screening scores.")
  }
  
  # Robust empirical noise level from the bulk of mediator scores
  score_med <- median(finite_scores)
  score_mad <- mad(finite_scores, center = score_med, constant = 1)
  
  if (!is.finite(score_mad) || score_mad <= 0) {
    score_mad <- IQR(finite_scores) / 1.349
  }
  if (!is.finite(score_mad) || score_mad <= 0) {
    score_mad <- max(1e-8, 1e-6 * max(finite_scores))
  }
  
  noise_threshold <- score_med + 2 * score_mad
  
  eligible <- ranked[
    is.finite(score[ranked]) &
      score[ranked] > noise_threshold
  ]
  
  if (length(eligible) == 0L) {
    # Keep only the strongest component if no score clearly separates from noise
    selected <- ranked[1L]
  } else {
    selected <- eligible[
      seq_len(min(max_selected, length(eligible)))
    ]
  }
  
  if (length(selected) < min_selected) {
    selected <- ranked[seq_len(min(min_selected, length(ranked)))]
  }
  
  k_select <- length(selected)
  selected
}

multi_y_var_selection <- function(X, M) {
  # 对多维Y的每一列进行变量选择，返回并集
  M=as.matrix(M)
  # 确保输入格式
  if (nrow(X) != nrow(M)) {
    stop("X和M必须有相同的行数")
  }
  
  # 获取列名
  y_names <- colnames(M)
  if (is.null(y_names)) {
    y_names <- paste0("M", 1:ncol(M))
  }
  
  # 初始化存储所有选择变量的列表
  all_selected_vars <- list()
  
  # 遍历M的每一列
  for (i in seq_along(y_names)) {
    y_current <- M[, i]
    
    # 调用已有的lasso_var_selection函数
    rf <- ranger(y ~ ., data = data.frame(X=X, y=y_current), importance = "impurity",
                 num.trees=100, num.threads=1)
    importance<-importance(rf)
    selected_vars <- names(importance)[importance > mean(importance)+sd(importance)]#
    # if(length(selected_vars)==0){selected_vars=colnames(X)} else{
    selected_vars <- gsub("^X\\.", "", selected_vars)  # 去除前缀
    selected_vars <- gsub("^M\\.", "", selected_vars)
    # }
    
    # 存储结果
    all_selected_vars[[y_names[i]]] <- selected_vars
  }
  
  # 取所有选择结果的并集
  union_vars <- unique(unlist(all_selected_vars))
  
  # 如果并集为空，返回所有X的变量名
  if (length(union_vars) == 0) {
    union_vars <- colnames(X)
  }
  
  # 清理变量名（与原始函数一致）
  union_vars <- gsub("^X\\.", "", union_vars)
  union_vars <- gsub("^M\\.", "", union_vars)
  
  return(list(
    union_vars = union_vars,
    selected_by_y = all_selected_vars
  ))
}
generate_num_tree <- function(n) {num_trees=50;return(num_trees)}
generate_tune_grid <- function(n_train, d, task = 0) {
  if (n_train < 1 || d < 1) stop("n_train and d must be positive integers.")
  if (!task %in% c(0, 1)) stop("task must be 0 (regression) or 1 (classification).")
  
  if (task == 0) {
    mtry0 <- max(1, floor(d / 3))
    splitrule_candidates <- c("variance", "extratrees")
    node0 <- max(5, floor(0.03 * n_train))
    node_cap <- 30
  } else {
    mtry0 <- max(1, floor(sqrt(d)))
    splitrule_candidates <- c("gini", "extratrees")
    node0 <- max(1, floor(0.01 * n_train))
    node_cap <- 10
  }
  
  mtry_candidates <- unique(pmax(1, pmin(d, round(c(0.5, 1, 1.5) * mtry0))))
  min_node_size_candidates <- unique(pmax(1, pmin(node_cap, round(c(0.5, 1) * node0))))
  
  tuneGrid <- expand.grid(
    mtry = mtry_candidates,
    min.node.size = min_node_size_candidates,
    splitrule = splitrule_candidates,
    stringsAsFactors = FALSE
  )
  
  tuneGrid <- unique(tuneGrid)
  rownames(tuneGrid) <- NULL
  tuneGrid
}

generate_tune_grid_otr_p <- function(n_train, d) {

  ## Start from the original tuning grid
  base_grid <- generate_tune_grid(
    n_train = n_train,
    d = d,
    task = 1
  )

  ## Extend only the upper range of mtry
    extra_mtry <- unique(
      pmax(
        1,
        pmin(
          d,
          round(c(
            2 * sqrt(d),
            3 * sqrt(d),
            d / 3,
            d / 2,
            d
          ))
        )
      )
    )

  ## Keep exactly the same node-size and split-rule candidates
  ## as in the original grid
  extra_grid <- expand.grid(
    mtry = extra_mtry,
    min.node.size = unique(base_grid$min.node.size),
    splitrule = unique(base_grid$splitrule),
    stringsAsFactors = FALSE
  )
    

  ## Original grid + additional larger-mtry candidates
  tune_grid <- unique(rbind(base_grid, extra_grid))

  tune_grid <- tune_grid[
    order(
      tune_grid$mtry,
      tune_grid$min.node.size,
      tune_grid$splitrule
    ),
  ]

  rownames(tune_grid) <- NULL

  return(tune_grid)
}                    

fit_mv_rf_cf <- function(X_train, M_train, K = 3,
                         mtry_grid = NULL, nodesize_grid = c(3, 5, 10, 50),
                         ntree_tune = 50, ntree_final = 50,
                         seed = 123) {
  X_train <- as.matrix(X_train)
  M_train <- as.matrix(M_train)
  n <- nrow(X_train)
  p <- ncol(X_train)
  d <- ncol(M_train)
  
  if (is.null(colnames(X_train))) colnames(X_train) <- paste0("X", seq_len(p))
  if (is.null(colnames(M_train))) colnames(M_train) <- paste0("M", seq_len(d))
  
  if (is.null(mtry_grid)) {
    mtry_grid <- unique(pmax(1L, pmin(p, round(c(sqrt(p), p/10, p/5, p/3)))))
  }
  nodesize_grid <- unique(pmax(1L, pmin(n - 1L, nodesize_grid)))
  tune_grid <- expand.grid(mtry = mtry_grid, nodesize = nodesize_grid)
  
  K <- min(K, n)
  if (K < 2) stop("Not enough observations for RF inner CV.")
  
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (has_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (has_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  
  set.seed(seed)
  fold_id <- sample(rep(seq_len(K), length.out = n))
  
  formula_multi <- as.formula(
    paste0(
      "cbind(", paste(colnames(M_train), collapse = ","), ") ~ ",
      paste(colnames(X_train), collapse = "+")
    )
  )
  
  cv_loss <- matrix(NA_real_, nrow = nrow(tune_grid), ncol = K)
  
  
  
  for (v in seq_len(K)) {
    idx_val <- which(fold_id == v)
    idx_fit <- which(fold_id != v)
    
    X_fit <- X_train[idx_fit, , drop = FALSE]
    X_val <- X_train[idx_val, , drop = FALSE]
    M_fit <- M_train[idx_fit, , drop = FALSE]
    M_val <- M_train[idx_val, , drop = FALSE]
    
    # Standardize validation MSE across mediator dimensions so that
    # high-variance mediators do not dominate hyperparameter selection.
    M_scale <- apply(M_fit, 2, sd)
    M_scale[!is.finite(M_scale) | M_scale <= 0] <- 1
    
    dat_fit <- data.frame(X_fit, M_fit, check.names = FALSE)
    
    for (g in seq_len(nrow(tune_grid))) {
      fit_g <- rfsrc(
        formula_multi,
        data = dat_fit,
        ntree = ntree_tune,
        block.size = d,
        importance = FALSE,
        forest = TRUE,
        splitrule = "mv",
        nodesize = tune_grid$nodesize[g],
        mtry = tune_grid$mtry[g],
        bootstrap = "by.root",
        seed = seed + 1000 * v
      )
      
      pred_g <- predict(
        fit_g,
        newdata = data.frame(X_val, check.names = FALSE)
      )
      
      pred_val <- sapply(seq_len(d), function(j) {
        pred_g$regrOutput[[j]][[1]]
      })
      if (d == 1) pred_val <- matrix(pred_val, ncol = 1)
      
      err_std <- sweep(pred_val - M_val, 2, M_scale, "/")
      cv_loss[g, v] <- mean(err_std^2)
    }
  }
  
  mean_cv_loss <- rowMeans(cv_loss)
  best_id <- which.min(mean_cv_loss)
  best_mtry <- tune_grid$mtry[best_id]
  best_nodesize <- tune_grid$nodesize[best_id]
  
  dat_full <- data.frame(X_train, M_train, check.names = FALSE)
  
  fit_final <- rfsrc(
    formula_multi,
    data = dat_full,
    ntree = ntree_final,
    block.size = d,
    importance = TRUE,
    forest = TRUE,
    splitrule = "mv",
    nodesize = best_nodesize,
    mtry = best_mtry,
    bootstrap = "by.root",
    seed = seed + 99999
  )
  
  oob_pred <- get.mv.predicted(fit_final, oob = TRUE)
  colnames(oob_pred) <- colnames(M_train)
  
  list(
    fit = fit_final,
    oob_pred = oob_pred,
    best_mtry = best_mtry,
    best_nodesize = best_nodesize,
    cv_loss = mean_cv_loss,
    tune_grid = tune_grid
  )
}

dmvnorm_eigen_log <- function(x, sigma) {
  x <- as.matrix(x)
  sigma <- as.matrix(sigma)

  if (ncol(x) != nrow(sigma) || nrow(sigma) != ncol(sigma)) {
    stop("Dimension mismatch in dmvnorm_eigen_log.")
  }
  if (any(!is.finite(x)) || any(!is.finite(sigma))) {
    stop("Non-finite input in dmvnorm_eigen_log.")
  }

  # Work with a plain numeric, exactly symmetric matrix.
  sigma_num <- matrix(
    as.numeric(sigma),
    nrow = nrow(sigma),
    ncol = ncol(sigma)
  )
  sigma_num <- 0.5 * (sigma_num + t(sigma_num))

  eig <- eigen(
    sigma_num,
    symmetric = TRUE
  )

  values <- as.numeric(eig$values)
  vectors <- eig$vectors

  if (any(!is.finite(values)) || any(!is.finite(vectors))) {
    stop("Non-finite eigendecomposition in dmvnorm_eigen_log.")
  }

  # Numerical safeguard only: floor eigenvalues at machine-precision scale.
  eig_floor <- max(
    1e-10,
    max(abs(values)) * ncol(sigma_num) * .Machine$double.eps
  )
  values <- pmax(values, eig_floor)

  # If sigma = V diag(lambda) V', then
  # x' Sigma^{-1} x = sum_j (x' v_j)^2 / lambda_j.
  z <- x %*% vectors
  quad <- rowSums(
    sweep(z^2, 2, values, "/")
  )

  logdet <- sum(log(values))
  d <- ncol(x)

  out <- -0.5 * (
    d * log(2 * base::pi) +
      logdet +
      quad
  )

  if (any(!is.finite(out))) {
    stop("Non-finite log density in dmvnorm_eigen_log.")
  }

  as.numeric(out)
}


stabilize_cov_eigen <- function(sigma) {
  sigma <- as.matrix(sigma)

  sigma_num <- matrix(
    as.numeric(sigma),
    nrow = nrow(sigma),
    ncol = ncol(sigma)
  )
  sigma_num <- 0.5 * (sigma_num + t(sigma_num))

  eig <- eigen(
    sigma_num,
    symmetric = TRUE
  )

  values <- as.numeric(eig$values)
  vectors <- eig$vectors

  if (any(!is.finite(values)) || any(!is.finite(vectors))) {
    stop("Non-finite eigendecomposition in stabilize_cov_eigen.")
  }

  eig_floor <- max(
    1e-10,
    max(abs(values)) * ncol(sigma_num) * .Machine$double.eps
  )
  values <- pmax(values, eig_floor)

  sigma_pd <- vectors %*%
    (values * t(vectors))

  sigma_pd <- 0.5 * (sigma_pd + t(sigma_pd))
  dimnames(sigma_pd) <- dimnames(sigma)

  sigma_pd
}

rmvnorm_eigen_custom <- function(n, mean, sigma) {
  mean <- as.numeric(mean)
  sigma <- stabilize_cov_eigen(sigma)

  eig <- eigen(
    sigma,
    symmetric = TRUE
  )

  values <- pmax(
    as.numeric(eig$values),
    max(1e-10, max(abs(eig$values)) * ncol(sigma) * .Machine$double.eps)
  )
  vectors <- eig$vectors

  z <- matrix(
    rnorm(n * length(mean)),
    nrow = n,
    ncol = length(mean)
  )

  z <- sweep(
    z,
    2,
    sqrt(values),
    "*"
  )

  draws <- z %*% t(vectors)
  sweep(draws, 2, mean, "+")
}

compute_shrinkage_cov <- function(residuals, lambda = "auto",
                                  ensure_positive_definite = TRUE) {
  if (identical(lambda, "auto")) {
    shrink_cov <- cov.shrink(residuals, verbose = FALSE)
  } else {
    shrink_cov <- cov.shrink(residuals, lambda = lambda, verbose = FALSE)
  }
  
  lambda_val <- attr(shrink_cov, "lambda")
  cov_matrix <- as.matrix(shrink_cov)
  class(cov_matrix) <- "matrix"
  cov_matrix <- 0.5 * (cov_matrix + t(cov_matrix))
  
  if (ensure_positive_definite) {
    eigen_vals <- eigen(cov_matrix, symmetric = TRUE, only.values = TRUE)$values
    if (any(eigen_vals <= 0)) {
      diag(cov_matrix) <- diag(cov_matrix) + 1e-6
    }
  }
  
  return(list(
    matrix = cov_matrix,
    lambda = lambda_val,
    lambda.estimated = attr(shrink_cov, "lambda.estimated"),
    diag_values = diag(cov_matrix)
  ))
}

# Legacy TR-only density helpers; QR retains the newer normalization.
tr_legacy_multi_y_var_selection <- function(X, M) {
  M=as.matrix(M)
  if (nrow(X) != nrow(M)) {
    stop("X和M必须有相同的行数")
  }
  y_names <- colnames(M)
  if (is.null(y_names)) {
    y_names <- paste0("M", seq_len(ncol(M)))
  }
  all_selected_vars <- list()
  for (i in seq_along(y_names)) {
    y_current <- M[, i]
    rf <- ranger(y ~ ., data = data.frame(X=X, y=y_current),
                 importance = "impurity", num.trees=50)
    importance<-importance(rf)
    selected_vars <- names(importance)[importance > mean(importance)+sd(importance)]
    if(length(selected_vars)==0){
      selected_vars=colnames(X)
    } else {
      selected_vars <- gsub("^X\\.", "", selected_vars)
      selected_vars <- gsub("^M\\.", "", selected_vars)
    }
    all_selected_vars[[y_names[i]]] <- selected_vars
  }
  union_vars <- unique(unlist(all_selected_vars))
  if (length(union_vars) == 0) {
    union_vars <- colnames(X)
  }
  union_vars <- gsub("^X\\.", "", union_vars)
  union_vars <- gsub("^M\\.", "", union_vars)
  return(list(
    union_vars = union_vars,
    selected_by_y = all_selected_vars
  ))
}

tr_legacy_normalized_truncate_q <- function(log_dens_A0, log_dens_A1, treated = NULL,
                                  lower_prob = 0.01, upper_prob = 0.99) {
  log_q <- log_dens_A0 - log_dens_A1
  q_raw <- exp(pmax(pmin(log_q, 700), -700))
  if (is.null(treated)) {
    treated <- rep(TRUE, length(q_raw))
  }
  treated <- as.logical(treated)
  ref <- q_raw[treated & is.finite(q_raw)]
  # ref <- q_raw[is.finite(q_raw)]
  if (length(ref) < 2) {
    return(q_raw)
  }
  q_bounds <- quantile(ref, probs = c(lower_prob, upper_prob), na.rm = TRUE)
  q_trunc <- pmax(pmin(q_raw, q_bounds[2]), q_bounds[1])
  norm_const <- mean(q_trunc[treated & is.finite(q_trunc)], na.rm = TRUE)
  if (!is.finite(norm_const) || norm_const <= 0) {
    norm_const <- 1
  }
  q_trunc / norm_const
}

select_mu_M_components <- function(X, M, y, top_m = min(5, ncol(M))) {
  top_m <- max(1L, min(top_m, ncol(M)))
  dat <- data.frame(X, M, y = y)
  rf <- ranger(y ~ ., data = dat, importance = "impurity", num.trees = 100)
  imp <- importance(rf)
  m_names <- colnames(M)
  m_imp <- imp[names(imp) %in% m_names]
  if (length(m_imp) == 0 || all(!is.finite(m_imp))) {
    return(seq_len(top_m))
  }
  selected_names <- names(sort(m_imp, decreasing = TRUE))[
    seq_len(min(top_m, length(m_imp)))
  ]
  match(selected_names, m_names)
}

knn_cond_log_density <- function(train_X, train_M, query_X, query_M,
                                 k = NULL, k_candidates = c(40, 60, 100, 150),
                                 cv_frac = 0.25, min_bw = 1e-3) {
  train_X <- as.matrix(train_X)
  query_X <- as.matrix(query_X)
  train_M <- as.matrix(train_M)
  query_M <- as.matrix(query_M)
  if (ncol(train_M) != ncol(query_M)) {
    stop("train_M and query_M must have the same number of columns")
  }
  n_train <- nrow(train_X)
  x_center <- colMeans(train_X)
  x_scale <- apply(train_X, 2, sd)
  x_scale[!is.finite(x_scale) | x_scale <= 0] <- 1
  train_X_std <- sweep(sweep(train_X, 2, x_center, "-"), 2, x_scale, "/")
  query_X_std <- sweep(sweep(query_X, 2, x_center, "-"), 2, x_scale, "/")

  estimate_log_density <- function(ref_X_std, ref_M, eval_X_std, eval_M, k_use) {
    leave_one_out <- nrow(ref_X_std) == nrow(eval_X_std) &&
      ncol(ref_X_std) == ncol(eval_X_std) &&
      nrow(ref_M) == nrow(eval_M) &&
      ncol(ref_M) == ncol(eval_M) &&
      identical(unname(ref_X_std), unname(eval_X_std)) &&
      identical(unname(ref_M), unname(eval_M))
    k_target <- min(k_use, nrow(ref_X_std) - as.integer(leave_one_out))
    if (k_target < 1) {
      stop("Not enough reference observations for leave-one-out KNN density estimation")
    }
    k_search <- min(nrow(ref_X_std), k_target + as.integer(leave_one_out))
    nn_idx <- FNN::get.knnx(ref_X_std, eval_X_std, k = k_search)$nn.index
    m_dim <- ncol(ref_M)
    log_density <- rep(NA_real_, nrow(eval_M))
    for (i in seq_len(nrow(eval_M))) {
      idx_i <- nn_idx[i, ]
      if (leave_one_out) {
        idx_i <- idx_i[idx_i != i]
      }
      if (length(idx_i) < k_target) {
        stop("Insufficient non-self neighbors for leave-one-out KNN density estimation")
      }
      idx_i <- idx_i[seq_len(k_target)]
      neigh_M <- ref_M[idx_i, , drop = FALSE]
      local_bw <- apply(neigh_M, 2, sd) * k_target^(-1 / (m_dim + 4))
      local_bw[!is.finite(local_bw) | local_bw < min_bw] <- min_bw
      log_norm <- sum(log(local_bw)) + m_dim * 0.5 * log(2 * pi)
      z <- sweep(neigh_M, 2, eval_M[i, ], "-")
      z <- sweep(z, 2, local_bw, "/")
      log_kernel <- -0.5 * rowSums(z^2) - log_norm
      max_log <- max(log_kernel)
      log_density[i] <- max_log + log(mean(exp(log_kernel - max_log)))
    }
    log_density
  }

  if (is.null(k)) {
    val_step <- max(2, round(1 / cv_frac))
    val_idx <- seq(1, n_train, by = val_step)
    if (length(val_idx) < 10 || length(val_idx) >= n_train) {
      k <- min(80, n_train)
    } else {
      train_idx <- setdiff(seq_len(n_train), val_idx)
      k_grid <- k_candidates[k_candidates < length(train_idx)]
      if (length(k_grid) == 0) {
        k <- min(80, n_train)
      } else {
        cv_scores <- sapply(k_grid, function(k_try) {
          mean(estimate_log_density(
            train_X_std[train_idx, , drop = FALSE],
            train_M[train_idx, , drop = FALSE],
            train_X_std[val_idx, , drop = FALSE],
            train_M[val_idx, , drop = FALSE],
            k_try
          ), na.rm = TRUE)
        })
        k <- k_grid[which.max(cv_scores)]
      }
    }
  }

  k <- min(k, n_train)
  estimate_log_density(train_X_std, train_M, query_X_std, query_M, k)
}

# Restore the 175_2 density-ratio module only. Keep its random screening
# separate from the RNG stream used by the other nuisance models and MC tau1.
estimate_tr_q_legacy <- function(X, M, A, Y, train_index, query_index,
                                 selected_vars_mu) {
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (has_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (has_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  train_A1 <- train_index[A[train_index] == 1]
  train_A0 <- train_index[A[train_index] == 0]
  selected_M <- match(intersect(selected_vars_mu, colnames(M)), colnames(M))
  selected_M <- selected_M[!is.na(selected_M)]
  if (length(selected_M) == 0) {
    selected_M <- select_mu_M_components(
      X[train_A1, , drop = FALSE], M[train_A1, , drop = FALSE], Y[train_A1],
      top_m = min(5, ncol(M))
    )
  }
  x_vars <- tr_legacy_multi_y_var_selection(
    X[train_index, , drop = FALSE], M[train_index, selected_M, drop = FALSE]
  )$union_vars
  log_dens_A1 <- knn_cond_log_density(
    X[train_A1, x_vars, drop = FALSE], M[train_A1, selected_M, drop = FALSE],
    X[query_index, x_vars, drop = FALSE], M[query_index, selected_M, drop = FALSE]
  )
  log_dens_A0 <- knn_cond_log_density(
    X[train_A0, x_vars, drop = FALSE], M[train_A0, selected_M, drop = FALSE],
    X[query_index, x_vars, drop = FALSE], M[query_index, selected_M, drop = FALSE]
  )
  tr_legacy_normalized_truncate_q(
    log_dens_A0, log_dens_A1, treated = A[query_index] == 1
  )
}

normalized_truncate_q <- function(
    log_q,
    ref_log_q,
    lower_prob = 0.01,
    upper_prob = 0.99) {
  
  log_q <- as.numeric(log_q)
  ref_log_q <- as.numeric(ref_log_q)
  
  if (length(ref_log_q) < 2L) {
    stop("Too few reference observations for q normalization/truncation.")
  }
  
  if (any(!is.finite(log_q)) || any(!is.finite(ref_log_q))) {
    stop("Non-finite log density-ratio values.")
  }
  
  # No truncation: pure self-normalization
  if (lower_prob == 0 && upper_prob == 1) {
    
    m <- max(ref_log_q)
    log_norm_const <- m + log(mean(exp(ref_log_q - m)))
    
    q_out <- exp(log_q - log_norm_const)
    
    if (any(!is.finite(q_out))) {
      stop("Non-finite normalized density-ratio values.")
    }
    
    attr(q_out, "lower_bound") <- 0
    attr(q_out, "upper_bound") <- Inf
    attr(q_out, "norm_const") <- exp(log_norm_const)
    
    return(q_out)
  }
  
  # Quantile truncation
  q_raw <- exp(log_q)
  q_ref <- exp(ref_log_q)
  
  q_bounds <- quantile(
    q_ref,
    probs = c(lower_prob, upper_prob),
    names = FALSE,
    na.rm = TRUE
  )
  
  lower <- q_bounds[1]
  upper <- q_bounds[2]
  
  q_trunc <- pmax(pmin(q_raw, upper), lower)
  q_ref_trunc <- pmax(pmin(q_ref, upper), lower)
  
  norm_const <- mean(q_ref_trunc)
  
  if (!is.finite(norm_const) || norm_const <= 0) {
    stop("Invalid q normalization constant.")
  }
  
  q_out <- q_trunc / norm_const
  
  attr(q_out, "lower_bound") <- lower
  attr(q_out, "upper_bound") <- upper
  attr(q_out, "norm_const") <- norm_const
  
  q_out
}


normalized_truncate_q_logstable <- function(
    log_q,
    ref_log_q,
    lower_prob = 0.01,
    upper_prob = 0.99) {

  log_q <- as.numeric(log_q)
  ref_log_q <- as.numeric(ref_log_q)

  if (length(ref_log_q) < 2L) {
    stop("Too few reference observations for q normalization/truncation.")
  }

  if (any(!is.finite(log_q)) || any(!is.finite(ref_log_q))) {
    stop("Non-finite log density-ratio values.")
  }

  # No truncation: stable self-normalization on the log scale.
  if (lower_prob == 0 && upper_prob == 1) {
    m <- max(ref_log_q)
    log_norm_const <- m + log(mean(exp(ref_log_q - m)))

    log_q_norm <- log_q - log_norm_const

    if (any(log_q_norm > log(.Machine$double.xmax))) {
      stop("Normalized q overflow in log-stable normalization.")
    }

    q_out <- exp(log_q_norm)

    if (any(!is.finite(q_out))) {
      stop("Non-finite normalized density-ratio values.")
    }

    attr(q_out, "lower_bound") <- 0
    attr(q_out, "upper_bound") <- Inf
    attr(q_out, "norm_const") <- exp(min(log_norm_const, log(.Machine$double.xmax)))

    return(q_out)
  }

  # Because exp(.) is monotone, truncating q at its empirical quantiles
  # is equivalent to truncating log(q) at the corresponding log-quantiles.
  log_bounds <- quantile(
    ref_log_q,
    probs = c(lower_prob, upper_prob),
    names = FALSE,
    na.rm = TRUE
  )

  lower_log <- log_bounds[1]
  upper_log <- log_bounds[2]

  log_q_trunc <- pmax(pmin(log_q, upper_log), lower_log)
  ref_log_q_trunc <- pmax(pmin(ref_log_q, upper_log), lower_log)

  # log(mean(exp(ref_log_q_trunc))) computed stably.
  m <- max(ref_log_q_trunc)
  log_norm_const <- m + log(mean(exp(ref_log_q_trunc - m)))

  log_q_norm <- log_q_trunc - log_norm_const

  # After common quantile truncation this should be moderate, but keep
  # an explicit numerical guard.
  log_q_norm <- pmin(
    pmax(log_q_norm, log(.Machine$double.xmin)),
    log(.Machine$double.xmax)
  )

  q_out <- exp(log_q_norm)

  if (any(!is.finite(q_out)) || any(q_out <= 0)) {
    stop("Invalid q after log-stable normalization.")
  }

  attr(q_out, "lower_bound") <- exp(max(lower_log, log(.Machine$double.xmin)))
  attr(q_out, "upper_bound") <- exp(min(upper_log, log(.Machine$double.xmax)))
  attr(q_out, "norm_const") <- exp(
    min(max(log_norm_const, log(.Machine$double.xmin)),
        log(.Machine$double.xmax))
  )

  q_out
}

make_stratified_fold_id <- function(A, K = 3, seed = 123) {
  A <- as.numeric(A)
  K <- min(K, min(table(A)))
  if (K < 2) stop("Not enough observations for inner CV.")
  
  set.seed(seed)
  fold_id <- integer(length(A))
  
  for (a in sort(unique(A))) {
    idx <- sample(which(A == a))
    fold_id[idx] <- rep(seq_len(K), length.out = length(idx))
  }
  
  fold_id
}     



# The QR basis is supplied by S_poly; never expand these columns again.
el_drm_raw_q <- function(M) {
  M <- as.matrix(M)
  if (ncol(M) < 1L || any(!is.finite(M))) {
    stop("QR requires at least one finite, precomputed basis column.")
  }
  if (is.null(colnames(M))) colnames(M) <- paste0("B", seq_len(ncol(M)))
  if (anyDuplicated(colnames(M))) stop("QR basis column names must be unique.")
  M
}

el_drm_q_transform <- function(M, basis_names = NULL) {
  raw_q <- el_drm_raw_q(M)
  
  if (!is.null(basis_names)) {
    if (!all(basis_names %in% colnames(raw_q))) {
      stop("Some selected basis terms are missing from raw q basis.")
    }
    
    raw_q <- raw_q[, basis_names, drop = FALSE]
  }
  
  q_center <- colMeans(raw_q)
  q_scale <- apply(raw_q, 2, sd)
  q_scale[!is.finite(q_scale) | q_scale <= 0] <- 1
  
  list(
    center = q_center,
    scale = q_scale,
    basis_names = colnames(raw_q)
  )
}

el_drm_apply_q <- function(M, q_transform) {
  raw_q <- el_drm_raw_q(M)
  
  raw_q <- raw_q[
    ,
    q_transform$basis_names,
    drop = FALSE
  ]
  
  sweep(
    sweep(raw_q, 2, q_transform$center, "-"),
    2,
    q_transform$scale,
    "/"
  )
}

el_drm_design <- function(X) {
  X <- as.matrix(X)
  cbind(`(Intercept)` = 1, X)
}

el_drm_unpack_theta <- function(theta, x_dim, q_dim) {
  theta0_end <- x_dim * q_dim
  list(
    theta0 = matrix(theta[seq_len(theta0_end)], nrow = x_dim, ncol = q_dim),
    theta1 = matrix(theta[theta0_end + seq_len(theta0_end)], nrow = x_dim, ncol = q_dim)
  )
}

el_drm_solve_lambda <- function(base, q_mat, lambda_start,
                                maxit = 50, tol = 1e-7, min_den = 1e-10) {
  lambda <- lambda_start
  for (iter in seq_len(maxit)) {
    den <- as.vector(base + q_mat %*% lambda)
    if (any(!is.finite(den)) || any(den <= min_den)) {
      lambda <- lambda * 0.5
      next
    }
    inv_den <- 1 / den
    score <- as.vector(crossprod(q_mat, inv_den))
    if (max(abs(score)) < tol) {
      return(lambda)
    }
    jac <- -crossprod(q_mat * as.vector(inv_den^2), q_mat)
    step <- tryCatch(qr.solve(jac, score), error = function(e) MASS::ginv(jac) %*% score)
    step <- as.vector(step)
    step_scale <- 1
    score_norm <- sum(score^2)
    accepted <- FALSE
    for (ls_iter in seq_len(20)) {
      lambda_try <- lambda - step_scale * step
      den_try <- as.vector(base + q_mat %*% lambda_try)
      if (all(is.finite(den_try)) && all(den_try > min_den)) {
        score_try <- as.vector(crossprod(q_mat, 1 / den_try))
        if (sum(score_try^2) <= score_norm || step_scale < 1e-4) {
          lambda <- lambda_try
          accepted <- TRUE
          break
        }
      }
      step_scale <- step_scale * 0.5
    }
    if (!accepted) {
      lambda <- lambda - step_scale * step
    }
  }
  lambda
}

el_drm_solve_nuisance <- function(theta, x_design, q_mat, A,
                                  inner_maxit = 60, inner_tol = 1e-5,
                                  lambda_maxit = 50, lambda_tol = 1e-7) {
  n_obs <- nrow(x_design)
  q_dim <- ncol(q_mat)
  x_dim <- ncol(x_design)
  theta_list <- el_drm_unpack_theta(theta, x_dim, q_dim)
  beta0 <- x_design %*% theta_list$theta0
  beta1 <- x_design %*% theta_list$theta1
  beta_obs <- beta0
  beta_obs[A == 1, ] <- beta1[A == 1, ]
  eta <- beta_obs %*% t(q_mat)
  alpha <- rep(0, n_obs)
  lambda <- rep(0, q_dim)
  p_vec <- rep(1 / n_obs, n_obs)
  for (iter in seq_len(inner_maxit)) {
    base_log <- matrixStats::colLogSumExps(sweep(eta, 1, alpha, "+"))
    base <- exp(pmin(base_log, 700))
    den <- as.vector(base + q_mat %*% lambda)
    if (any(!is.finite(den)) || any(den <= 0)) {
      stop("EL-DRM denominator became non-positive")
    }
    p_new <- 1 / den
    alpha_new <- -matrixStats::rowLogSumExps(sweep(eta, 2, log(p_new), "+"))
    base_log_new <- matrixStats::colLogSumExps(sweep(eta, 1, alpha_new, "+"))
    base_new <- exp(pmin(base_log_new, 700))
    lambda_new <- el_drm_solve_lambda(base_new, q_mat, lambda,
                                      maxit = lambda_maxit, tol = lambda_tol)
    den_new <- as.vector(base_new + q_mat %*% lambda_new)
    if (any(!is.finite(den_new)) || any(den_new <= 0)) {
      stop("EL-DRM lambda update produced non-positive denominator")
    }
    p_final <- 1 / den_new
    diff <- max(abs(alpha_new - alpha), abs(lambda_new - lambda), abs(p_final - p_vec))
    alpha <- alpha_new
    lambda <- lambda_new
    p_vec <- p_final
    if (diff < inner_tol) {
      break
    }
  }
  list(p = p_vec, alpha = alpha, lambda = lambda, eta = eta,
       theta0 = theta_list$theta0, theta1 = theta_list$theta1)
}

el_drm_profile <- function(theta, x_design, q_mat, A,
                           inner_maxit = 60, inner_tol = 1e-5,
                           lambda_maxit = 50, lambda_tol = 1e-7) {
  nuisance <- el_drm_solve_nuisance(theta, x_design, q_mat, A,
                                    inner_maxit = inner_maxit,
                                    inner_tol = inner_tol,
                                    lambda_maxit = lambda_maxit,
                                    lambda_tol = lambda_tol)
  n_obs <- nrow(x_design)
  diag_eta <- nuisance$eta[cbind(seq_len(n_obs), seq_len(n_obs))]
  loglik <- sum(log(nuisance$p)) + sum(nuisance$alpha) + sum(diag_eta)
  log_weights <- sweep(sweep(nuisance$eta, 1, nuisance$alpha, "+"),
                       2, log(nuisance$p), "+")
  weights <- exp(pmax(pmin(log_weights, 700), -700))
  expected_q <- weights %*% q_mat
  resid_q <- q_mat - expected_q
  grad0 <- matrix(0, nrow = ncol(x_design), ncol = ncol(q_mat))
  grad1 <- grad0
  idx0 <- which(A == 0)
  idx1 <- which(A == 1)
  if (length(idx0) > 0) {
    grad0 <- crossprod(x_design[idx0, , drop = FALSE], resid_q[idx0, , drop = FALSE])
  }
  if (length(idx1) > 0) {
    grad1 <- crossprod(x_design[idx1, , drop = FALSE], resid_q[idx1, , drop = FALSE])
  }
  list(loglik = loglik, gradient = c(as.vector(grad0), as.vector(grad1)),
       nuisance = nuisance)
}

el_drm_fit_marginal_p <- function(q_mat, train_A,
                                  marginal_maxit = 80,
                                  inner_maxit = 60, inner_tol = 1e-5,
                                  lambda_maxit = 50, lambda_tol = 1e-7) {
  x_marginal <- matrix(1, nrow = nrow(q_mat), ncol = 1)
  train_A <- as.numeric(train_A)
  q_dim <- ncol(q_mat)
  theta_init <- rep(0, 2 * q_dim)
  last_profile <- NULL
  objective <- function(theta) {
    prof <- tryCatch(el_drm_profile(theta, x_marginal, q_mat, train_A,
                                    inner_maxit = inner_maxit,
                                    inner_tol = inner_tol,
                                    lambda_maxit = lambda_maxit,
                                    lambda_tol = lambda_tol),
                     error = function(e) NULL)
    if (is.null(prof) || !is.finite(prof$loglik)) {
      return(1e100)
    }
    last_profile <<- prof
    -prof$loglik
  }
  gradient <- function(theta) {
    prof <- tryCatch(el_drm_profile(theta, x_marginal, q_mat, train_A,
                                    inner_maxit = inner_maxit,
                                    inner_tol = inner_tol,
                                    lambda_maxit = lambda_maxit,
                                    lambda_tol = lambda_tol),
                     error = function(e) NULL)
    if (is.null(prof) || !is.finite(prof$loglik) ||
        any(!is.finite(prof$gradient))) {
      return(rep(0, length(theta)))
    }
    last_profile <<- prof
    -prof$gradient
  }
  opt <- optim(
    theta_init, objective, gradient,
    method = "L-BFGS-B",
    control = list(maxit = marginal_maxit)
  )
  
  final_profile <- el_drm_profile(
    opt$par, x_marginal, q_mat, train_A,
    inner_maxit = inner_maxit,
    inner_tol = inner_tol,
    lambda_maxit = lambda_maxit,
    lambda_tol = lambda_tol
  )
  
  p_mar <- final_profile$nuisance$p
  p_mar / sum(p_mar)
}

el_drm_approx_loss_grad <- function(theta, x_design, q_mat, A, p_mar) {
  theta_list <- el_drm_unpack_theta(theta, ncol(x_design), ncol(q_mat))
  log_p <- log(pmax(p_mar, .Machine$double.xmin))
  loss <- 0
  grad0 <- matrix(0, nrow = ncol(x_design), ncol = ncol(q_mat))
  grad1 <- grad0
  add_group <- function(idx, theta_k) {
    if (length(idx) == 0) {
      return(list(loss = 0, grad = matrix(0, nrow = ncol(x_design), ncol = ncol(q_mat))))
    }
    x_g <- x_design[idx, , drop = FALSE]
    q_obs <- q_mat[idx, , drop = FALSE]
    beta_g <- x_g %*% theta_k
    logits <- sweep(beta_g %*% t(q_mat), 2, log_p, "+")
    log_norm <- matrixStats::rowLogSumExps(logits)
    weights <- exp(sweep(logits, 1, log_norm, "-"))
    expected_q <- weights %*% q_mat
    list(
      loss = sum(log_norm - rowSums(beta_g * q_obs)),
      grad = crossprod(x_g, expected_q - q_obs)
    )
  }
  group0 <- add_group(which(A == 0), theta_list$theta0)
  group1 <- add_group(which(A == 1), theta_list$theta1)
  loss <- group0$loss + group1$loss
  grad0 <- group0$grad
  grad1 <- group1$grad
  list(loss = loss, gradient = c(as.vector(grad0), as.vector(grad1)))
}

el_drm_approx_eval_loss <- function(theta, x_eval, q_eval, A_eval, q_ref, p_ref) {
  theta_list <- el_drm_unpack_theta(theta, ncol(x_eval), ncol(q_ref))
  log_p <- log(pmax(p_ref, .Machine$double.xmin))
  add_group <- function(idx, theta_k) {
    if (length(idx) == 0) {
      return(0)
    }
    x_g <- x_eval[idx, , drop = FALSE]
    q_g <- q_eval[idx, , drop = FALSE]
    beta_g <- x_g %*% theta_k
    logits <- sweep(beta_g %*% t(q_ref), 2, log_p, "+")
    sum(matrixStats::rowLogSumExps(logits) - rowSums(beta_g * q_g))
  }
  add_group(which(A_eval == 0), theta_list$theta0) +
    add_group(which(A_eval == 1), theta_list$theta1)
}

# Select supplied basis columns involving only the retained mediators.
# This only subsets S_poly; it never constructs polynomial terms.
sdr_existing_mediator_basis <- function(S_poly, mediator_names, selected_names) {
  column_names <- colnames(S_poly)
  contains <- vapply(mediator_names, function(nm) {
    vapply(strsplit(column_names, "_", fixed = TRUE), function(parts) nm %in% parts, logical(1))
  }, logical(length(column_names)))
  retained <- rowSums(contains[, mediator_names %in% selected_names, drop = FALSE]) > 0L
  excluded <- rowSums(contains[, !mediator_names %in% selected_names, drop = FALSE]) > 0L
  S_poly[, retained & !excluded, drop = FALSE]
}

select_el_drm_basis <- function(
    X, M, A, linear_names,
    stability_reps = 10L,
    subsample_frac = 0.8,
    cv_folds = 5L,
    alpha = 0.5,
    stability_cutoff = 0.8,
    linear_penalty_factor = 0.25,
    seed = 123) {
  
  X <- as.matrix(X)
  M <- as.matrix(M)
  A <- as.numeric(A)
  
  full_basis <- el_drm_raw_q(M)
  
  # M contains already-expanded S_poly columns; original mediator terms
  # remain mandatory, and all supplied higher-order terms are candidates.
  if (!all(linear_names %in% colnames(full_basis))) stop("Missing linear mediator columns.")
  quadratic_names <- setdiff(colnames(full_basis), linear_names)
  full_basis <- full_basis[, c(linear_names, quadratic_names), drop = FALSE]

  if (length(linear_names) == 0L) {
    stop("No linear mediator basis found.")
  }
  
  if (length(quadratic_names) == 0L) {
    return(list(
      selected_basis = linear_names,
      linear_basis = linear_names,
      selected_quadratic = character(0),
      quadratic_frequency = numeric(0),
      successful_reps = 0L
    ))
  }
  
  selection_count <- setNames(
    rep(0L, length(quadratic_names)),
    quadratic_names
  )
  
  idx1 <- which(A == 1)
  idx0 <- which(A == 0)
  
  if (length(idx1) < 20L || length(idx0) < 20L) {
    stop("Too few observations for EL-DRM basis screening.")
  }
  
  successful_reps <- 0L
  
  for (r in seq_len(stability_reps)) {
    
    set.seed(seed + r)
    
    sub1 <- sample(
      idx1,
      size = floor(subsample_frac * length(idx1)),
      replace = FALSE
    )
    
    sub0 <- sample(
      idx0,
      size = floor(subsample_frac * length(idx0)),
      replace = FALSE
    )
    
    idx_sub <- c(sub1, sub0)
    
    X_sub <- X[idx_sub, , drop = FALSE]
    B_sub <- full_basis[idx_sub, , drop = FALSE]
    A_sub <- A[idx_sub]
    
    Z_sub <- cbind(
      X_sub,
      B_sub
    )
    
    # Linear mediator terms receive a weaker penalty,
    # but are always retained in the final basis anyway.
    penalty_factor <- c(
      rep(1, ncol(X_sub)),
      rep(linear_penalty_factor, length(linear_names)),
      rep(1, length(quadratic_names))
    )
    
    fold_id <- make_stratified_fold_id(
      A_sub,
      K = cv_folds,
      seed = seed + 1000L + r
    )
    
    fit_cv <- tryCatch(
      cv.glmnet(
        x = Z_sub,
        y = A_sub,
        family = "binomial",
        alpha = alpha,
        foldid = fold_id,
        penalty.factor = penalty_factor,
        standardize = TRUE,
        type.measure = "deviance",
        parallel = FALSE
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit_cv)) {
      next
    }
    
    beta_hat <- tryCatch(
      coef(fit_cv, s = "lambda.1se"),
      error = function(e) NULL
    )
    
    if (is.null(beta_hat)) {
      next
    }
    
    beta_vec <- as.vector(beta_hat)
    beta_names <- rownames(beta_hat)
    
    selected_now <- beta_names[
      is.finite(beta_vec) & abs(beta_vec) > 1e-8
    ]
    
    selected_quad_now <- intersect(
      selected_now,
      quadratic_names
    )
    
    successful_reps <- successful_reps + 1L
    
    if (length(selected_quad_now) > 0L) {
      selection_count[selected_quad_now] <-
        selection_count[selected_quad_now] + 1L
    }
  }
  
  if (successful_reps < max(3L, ceiling(stability_reps / 2))) {
    stop(
      "Too few successful EL-DRM basis-screening repetitions: ",
      successful_reps, "/", stability_reps
    )
  }
  
  selection_frequency <-
    selection_count / successful_reps
  
  selected_quadratic <- names(selection_frequency)[
    selection_frequency >= stability_cutoff
  ]
  
  selected_basis <- c(
    linear_names,
    selected_quadratic
  )
  
  list(
    selected_basis = selected_basis,
    linear_basis = linear_names,
    selected_quadratic = selected_quadratic,
    quadratic_frequency = sort(
      selection_frequency,
      decreasing = TRUE
    ),
    successful_reps = successful_reps
  )
}

# Screen the tensor terms already present in the final EL-DRM model.
# Products here are screening regressors only: q's supplied basis is unchanged.
select_el_drm_interactions <- function(
    X, M, A, linear_names, initial_x_vars,
    stability_reps = 10L, subsample_frac = 0.8, cv_folds = 5L,
    alpha = 0.5, stability_cutoff = 0.8,
    linear_penalty_factor = 0.25, seed = 123) {
  X <- as.matrix(X); M <- el_drm_raw_q(M); A <- as.numeric(A)
  stopifnot(nrow(X) == nrow(M), length(A) == nrow(X),
            all(linear_names %in% colnames(M)),
            all(initial_x_vars %in% colnames(X)),
            !anyDuplicated(colnames(X)), all(is.finite(X)))
  px <- ncol(X); pb <- ncol(M)
  # Include ALL X_poly columns here, including ones missed by the RF screen.
  pair_x <- rep(seq_len(px), times = pb)
  pair_b <- rep(seq_len(pb), each = px)
  products <- X[, pair_x, drop = FALSE] * M[, pair_b, drop = FALSE]
  Z <- cbind(X, M, products)
  rm(products)
  # Internal unique names avoid ambiguity from underscore-rich polynomial names.
  colnames(Z) <- paste0("screen_", seq_len(ncol(Z)))
  penalty <- c(rep(1, px), ifelse(colnames(M) %in% linear_names,
                                linear_penalty_factor, 1), rep(1, length(pair_x)))
  count <- integer(ncol(Z))
  idx1 <- which(A == 1); idx0 <- which(A == 0)
  if (min(length(idx1), length(idx0)) < 20L) stop("Too few observations for interaction screening.")
  successful_reps <- 0L
  for (r in seq_len(stability_reps)) {
    set.seed(seed + r)
    idx <- c(sample(idx1, floor(subsample_frac * length(idx1))),
             sample(idx0, floor(subsample_frac * length(idx0))))
    fold_id <- make_stratified_fold_id(A[idx], K = cv_folds, seed = seed + 1000L + r)
    fit <- tryCatch(glmnet::cv.glmnet(
      x = Z[idx, , drop = FALSE], y = A[idx], family = "binomial",
      alpha = alpha, foldid = fold_id, penalty.factor = penalty,
      standardize = TRUE, type.measure = "deviance", parallel = FALSE
    ), error = function(e) NULL)
    if (is.null(fit)) next
    beta <- as.matrix(coef(fit, s = "lambda.1se"))[-1L, 1L]
    if (length(beta) != ncol(Z) || any(!is.finite(beta))) next
    count <- count + as.integer(abs(beta) > 1e-8)
    successful_reps <- successful_reps + 1L
  }
  if (successful_reps < max(3L, ceiling(stability_reps / 2))) {
    stop("Too few successful interaction-screening repetitions: ", successful_reps)
  }
  frequency <- count / successful_reps
  main_frequency <- setNames(frequency[px + seq_len(pb)], colnames(M))
  pair_frequency <- frequency[px + pb + seq_along(pair_x)]
  stable_pairs <- which(pair_frequency >= stability_cutoff)
  retained_basis <- unique(c(linear_names,
    names(main_frequency)[main_frequency >= stability_cutoff],
    colnames(M)[pair_b[stable_pairs]]))
  # Preserve RF-selected X, and rescue X contributing stable main/interaction terms.
  retained_x <- unique(c(initial_x_vars,
    colnames(X)[frequency[seq_len(px)] >= stability_cutoff],
    colnames(X)[pair_x[stable_pairs]]))
  observed_pairs <- which(pair_frequency > 0)
  interaction_frequency <- data.frame(
    X = colnames(X)[pair_x[observed_pairs]],
    basis = colnames(M)[pair_b[observed_pairs]],
    frequency = pair_frequency[observed_pairs],
    retained = pair_frequency[observed_pairs] >= stability_cutoff,
    stringsAsFactors = FALSE)
  list(selected_basis = retained_basis, selected_x = retained_x,
       linear_basis = linear_names,
       selected_quadratic = setdiff(retained_basis, linear_names),
       quadratic_frequency = sort(main_frequency[!names(main_frequency) %in% linear_names], decreasing = TRUE),
       interaction_frequency = interaction_frequency,
       selected_interactions = interaction_frequency[interaction_frequency$retained, , drop = FALSE],
       candidate_interactions = length(pair_x), successful_reps = successful_reps)
}

el_drm_fit_density_ratio <- function(
    train_X, train_M, train_A,
    basis_names = NULL,
    ridge_grid = c(
      0, 0.01, 0.03, 0.05, 0.1,
      0.2, 0.3, 0.5, 1
    ),
    ridge_cv_folds = 3,
    ridge_cv_seed = 123,
    eta = NULL,
    outer_maxit = 80,
    marginal_maxit = 80,
    inner_maxit = 60,
    inner_tol = 1e-5,
    lambda_maxit = 50,
    lambda_tol = 1e-7) {
  train_X <- as.matrix(train_X)
  train_M <- as.matrix(train_M)
  train_A <- as.numeric(train_A)
  
  # Fit theta with ridge penalty:
  # loss(theta)/n + ridge_lambda/2 * ||theta||_2^2
  fit_theta <- function(x_fit, q_fit, A_fit, p_fit,
                        ridge_lambda = 0, eta_fit = NULL) {
    n_fit <- nrow(x_fit)
    x_dim_fit <- ncol(x_fit)
    q_dim_fit <- ncol(q_fit)
    
    # # Do not penalize intercept-related coefficients
    # penalty_mask <- rep(1, 2 * x_dim_fit * q_dim_fit)
    
    # intercept_idx0 <- 1 + (0:(q_dim_fit - 1)) * x_dim_fit
    # intercept_idx1 <- x_dim_fit * q_dim_fit + intercept_idx0
    
    # penalty_mask[c(intercept_idx0, intercept_idx1)] <- 0
    
    if (is.null(eta_fit)) {
      eta_fit <- 1 / max(q_dim_fit + ridge_lambda, 1)
    }
    
    theta_cur <- rep(0, 2 * x_dim_fit * q_dim_fit)
    y_cur <- theta_cur
    t_cur <- 1
    eta_cur <- eta_fit    
    
    objective <- function(theta, prof = NULL) {
      if (is.null(prof)) {
        prof <- el_drm_approx_loss_grad(theta, x_fit, q_fit, A_fit, p_fit)
      }
      if (is.null(prof) || !is.finite(prof$loss)) return(Inf)
      prof$loss / n_fit + 0.5 * ridge_lambda * sum(theta^2)
      # prof$loss / n_fit + 0.5 * ridge_lambda * sum((penalty_mask * theta)^2)
    }
    
    prof_cur <- el_drm_approx_loss_grad(theta_cur, x_fit, q_fit, A_fit, p_fit)
    objective_prev <- objective(theta_cur, prof_cur)
    
    for (iter in seq_len(outer_maxit)) {
      
      # Evaluate objective/gradient at the FISTA extrapolation point
      prof_y <- el_drm_approx_loss_grad(y_cur, x_fit, q_fit, A_fit, p_fit)
      objective_y <- objective(y_cur, prof_y)
      
      # Restart FISTA if extrapolation itself is unstable or worse
      if (is.null(prof_y) || !is.finite(objective_y) ||
          any(!is.finite(prof_y$gradient)) ||
          objective_y > objective_prev + 1e-10) {
        
        y_cur <- theta_cur
        t_cur <- 1
        
        prof_y <- el_drm_approx_loss_grad(y_cur, x_fit, q_fit, A_fit, p_fit)
        objective_y <- objective(y_cur, prof_y)
        
        if (is.null(prof_y) || !is.finite(objective_y) ||
            any(!is.finite(prof_y$gradient))) {
          break
        }
      }
      
      grad <- prof_y$gradient / n_fit + ridge_lambda * y_cur
      
      # Backtracking line search
      eta_trial <- eta_cur
      accepted <- FALSE
      
      for (bt in seq_len(20)) {
        theta_trial <- y_cur - eta_trial * grad
        
        prof_trial <- tryCatch(
          el_drm_approx_loss_grad(theta_trial, x_fit, q_fit, A_fit, p_fit),
          error = function(e) NULL
        )
        
        objective_trial <- if (is.null(prof_trial)) {
          Inf
        } else {
          objective(theta_trial, prof_trial)
        }
        
        # Armijo sufficient-decrease condition
        if (is.finite(objective_trial) &&
            objective_trial <= objective_y -
            1e-4 * eta_trial * sum(grad^2)) {
          accepted <- TRUE
          break
        }
        
        eta_trial <- eta_trial / 2
      }
      
      if (!accepted) {
        break
      }
      
      theta_next <- theta_trial
      objective_next <- objective_trial
      
      # Keep the accepted smaller step size for subsequent iterations
      eta_cur <- eta_trial
      
      if (max(abs(theta_next - theta_cur)) < 1e-5) {
        theta_cur <- theta_next
        objective_prev <- objective_next
        break
      }
      
      t_next <- (1 + sqrt(1 + 4 * t_cur^2)) / 2
      y_next <- theta_next +
        ((t_cur - 1) / t_next) * (theta_next - theta_cur)
      
      theta_cur <- theta_next
      y_cur <- y_next
      t_cur <- t_next
      objective_prev <- objective_next
    }
    
    prof_final <- el_drm_approx_loss_grad(theta_cur, x_fit, q_fit, A_fit, p_fit)
    
    list(
      theta = theta_cur,
      objective = objective(theta_cur, prof_final),
      loss = prof_final$loss / n_fit,
      ridge_lambda = ridge_lambda,
      # eta = eta_fit
      eta = eta_cur  
    )
  }
  
  # ---------------- inner CV for ridge lambda ----------------
  fold_id <- make_stratified_fold_id(
    train_A,
    K = ridge_cv_folds,
    seed = ridge_cv_seed
  )
  K_inner <- max(fold_id)
  
  cv_loss <- matrix(
    Inf,
    nrow = length(ridge_grid),
    ncol = K_inner
  )
  
  
  for (v in seq_len(K_inner)) {
    idx_val <- which(fold_id == v)
    idx_fit <- which(fold_id != v)
    
    X_fit <- train_X[idx_fit, , drop = FALSE]
    X_val <- train_X[idx_val, , drop = FALSE]
    M_fit <- train_M[idx_fit, , drop = FALSE]
    M_val <- train_M[idx_val, , drop = FALSE]
    A_fit <- train_A[idx_fit]
    A_val <- train_A[idx_val]
    
    x_fit <- el_drm_design(X_fit)
    x_val <- el_drm_design(X_val)
    
    # All transformations must be estimated using inner-training data only
    q_transform_v <- el_drm_q_transform(
      M_fit,
      basis_names = basis_names
    )
    q_fit <- el_drm_apply_q(M_fit, q_transform_v)
    q_val <- el_drm_apply_q(M_val, q_transform_v)
    
    p_fit <- el_drm_fit_marginal_p(
      q_fit, A_fit,
      marginal_maxit = marginal_maxit,
      inner_maxit = inner_maxit,
      inner_tol = inner_tol,
      lambda_maxit = lambda_maxit,
      lambda_tol = lambda_tol
    )
    
    for (l in seq_along(ridge_grid)) {
      ridge_l <- ridge_grid[l]
      
      fit_l <- tryCatch(
        fit_theta(
          x_fit = x_fit,
          q_fit = q_fit,
          A_fit = A_fit,
          p_fit = p_fit,
          ridge_lambda = ridge_l,
          eta_fit = NULL
        ),
        error = function(e) NULL
      )
      
      if (!is.null(fit_l) && all(is.finite(fit_l$theta))) {
        
        cv_loss[l, v] <- tryCatch(
          el_drm_approx_eval_loss(
            theta = fit_l$theta,
            x_eval = x_val,
            q_eval = q_val,
            A_eval = A_val,
            q_ref = q_fit,
            p_ref = p_fit
          ) / length(idx_val),
          error = function(e) Inf
        )
        
      }
    }
  }
  
  mean_cv_loss <- rowMeans(cv_loss)
  
  if (all(!is.finite(mean_cv_loss))) {
    stop("All ridge parameters failed in EL-DRM inner CV.")
  }
  
  mean_cv_loss[!is.finite(mean_cv_loss)] <- Inf
  best_idx <- which.min(mean_cv_loss)
  best_ridge <- ridge_grid[best_idx]
  
  
  # ---------------- refit on all outer-training data ----------------
  x_design <- el_drm_design(train_X)
  q_transform <- el_drm_q_transform(
    train_M,
    basis_names = basis_names
  )
  q_mat <- el_drm_apply_q(train_M, q_transform)
  
  p_mar <- el_drm_fit_marginal_p(
    q_mat, train_A,
    marginal_maxit = marginal_maxit,
    inner_maxit = inner_maxit,
    inner_tol = inner_tol,
    lambda_maxit = lambda_maxit,
    lambda_tol = lambda_tol
  )
  
  fit <- fit_theta(
    x_fit = x_design,
    q_fit = q_mat,
    A_fit = train_A,
    p_fit = p_mar,
    ridge_lambda = best_ridge,
    eta_fit = eta
  )
  
  full_train_log_ratio <- el_drm_predict_log_ratio(
    list(
      theta = fit$theta,
      p_mar = p_mar,
      q_transform = q_transform,
      q_train = q_mat,
      x_dim = ncol(x_design),
      q_dim = ncol(q_mat)
    ),
    train_X,
    train_M
  )
  
  list(
    theta = fit$theta,
    p_mar = p_mar,
    q_transform = q_transform,
    q_train = q_mat,
    x_dim = ncol(x_design),
    q_dim = ncol(q_mat),
    basis_names = q_transform$basis_names,
    eta = fit$eta,
    objective = fit$objective,
    ridge_lambda = best_ridge,
    ridge_cv_loss = mean_cv_loss,
    ridge_grid = ridge_grid,
    full_train_log_ratio = full_train_log_ratio
  )
}

el_drm_predict_log_ratio <- function(fit, query_X, query_M) {
  x_query <- el_drm_design(query_X)
  q_query <- el_drm_apply_q(query_M, fit$q_transform)
  theta_list <- el_drm_unpack_theta(fit$theta, fit$x_dim, fit$q_dim)
  beta0_query <- x_query %*% theta_list$theta0
  beta1_query <- x_query %*% theta_list$theta1
  eta0_train <- beta0_query %*% t(fit$q_train)
  eta1_train <- beta1_query %*% t(fit$q_train)
  log_p <- log(pmax(fit$p_mar, .Machine$double.xmin))
  alpha0 <- -matrixStats::rowLogSumExps(sweep(eta0_train, 2, log_p, "+"))
  alpha1 <- -matrixStats::rowLogSumExps(sweep(eta1_train, 2, log_p, "+"))
  alpha0 - alpha1 + rowSums((beta0_query - beta1_query) * q_query)
}

predict_caret_regression_final <- function(fit, newdata) {
  pred <- predict(fit$finalModel, data = newdata)
  if (!is.null(pred$predictions)) {
    as.numeric(pred$predictions)
  } else {
    as.numeric(predict(fit, newdata = newdata))
  }
}


# Compatibility with the ACTG175 endpoint regressions.
var_selection <- lasso_var_selection

#---------------Choose data---------------
data("ACTG175")
df0 <- ACTG175
df <- df0[df0$offtrt == 0 & df0$r == 1, ]

# RCT
# oprior: only 23 1;  zprior: always 1;  
useless <- c("pidnum", "offtrt", "r", "days", "cens", "oprior", "zprior")
X_maybe <- c("z30", "strat")
X_conti <- c( c("age", "wtkg", "karnof", "cd40", "cd80"))#, c("preanti")
X_binar <- c( c("hemo", "homo", "drugs", "race", "gender", "str2", "symptom") )
X_vars <- c(X_conti, X_binar)

A_vars <- c("treat")
A_vars_arm <- c("arms")

M_vars <- c("cd420", "cd820")
Y_vars <- c("cd496")

df=df[df$arms!=3,]
A_all=df$treat #0: 单药; 1: 联合治疗

X_all <- df[,c(X_conti,X_binar)]  
M_all <- log(df[,M_vars])
# colnames(M_all) <- c("M1","M2")
Y_all <- log(df[,Y_vars])

S_all <- cbind(X_all,M_all)

# X_expanded <- expand_covariates(df, X_conti, X_binar)
# X_expanded <- expand_covariates_ck(df, X_conti, X_binar)
# S_expanded <- cbind(X_expanded, M_all)

result <- expand_covariates_with_mediation(
  data = df,
  X_conti = X_conti,
  X_binar = X_binar,
  M_conti = M_vars,
  add_conti_binar_interactions = TRUE,
  remove_zero_var = TRUE
)
X_expanded <- result$X_expanded
M_expanded <- result$M_expanded
S_expanded <- cbind(X_expanded, M_expanded)

# QR uses log mediators, as in 175_2 and the original 175_3.
# Transform the mediators BEFORE creating powers and interactions.
df_sdr_log <- df
df_sdr_log[, M_vars] <- M_all
stopifnot(all(is.finite(as.matrix(df_sdr_log[, M_vars, drop = FALSE]))))
result_sdr_log <- expand_covariates_with_mediation(
  data = df_sdr_log,
  X_conti = X_conti,
  X_binar = X_binar,
  M_conti = M_vars,
  add_conti_binar_interactions = TRUE,
  remove_zero_var = TRUE
)
stopifnot(identical(result_sdr_log$X_expanded, X_expanded))
S_sdr_log_expanded <- cbind(X_expanded, result_sdr_log$M_expanded)

print(ncol(X_expanded))
print(ncol(M_expanded))

N_all=length(Y_all);sw1=1;sw2=1
n_mc=50L
basis_stability_reps=5L
Kmb=5;
# N=floor(N_all/Kmb)*Kmb
N=floor(N_all/Kmb)*Kmb
current_time <- format(Sys.time(), "%Y%m%d_%H_%M")
seed=2026;limit=31;tmax=100
set.seed(seed)

#range of tuning parameters
rho.range0_default=c(seq(from=0.3,to=2,length.out=96),50,100,1000,10000)
rho.range02_default=c(seq(from=0.3,to=2,length.out=96),50,100,1000,10000)
rho.range0b_default=c(seq(from=0.5,to=4,length.out=96),50,100,1000,10000)
rho.range1_default=c(seq(from=0.5,to=3,length.out=96),50,100,1000,10000)
rho.range12_default=c(seq(from=0.4,to=2,length.out=96),50,100,1000,10000)

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1",
  R_RANGER_NUM_THREADS = "1",
  RF_CORES = "1"
)
options(ranger.num.threads = 1L, rf.cores = 1L)
ctrl_pi=trainControl(method="cv", number=5, summaryFunction=twoClassSummary, classProbs=TRUE,
                     allowParallel=FALSE)
ctrl_mu=trainControl(method="cv", number=5, summaryFunction=defaultSummary,classProbs=FALSE,
                     savePredictions=FALSE, allowParallel=FALSE)


# ctrl_pi <- trainControl(method = "oob", classProbs=TRUE)
# ctrl_mu <- trainControl(method = "oob")  # 回归任务只需这个

#--------------For record---------------------
log_dir <- "logs_cd"
if (!dir.exists(log_dir)) {
  dir.create(log_dir)
}
log_file<-file.path(log_dir, paste0("real_2_cd_", current_time, ".txt"))
log_appender(appender_file(log_file))
log_formatter(formatter_paste)
log_threshold(INFO)  # 设置日志级别为 INFO
log_info("Logger initialized and ready to record logs.")
cat("Simulation begun\n")

#------preparation--------------------------
num_mu_est=10;num_pi_est=10;num_q_est=10;num_p_est=6;num_totest=num_pi_est
grid=data.frame(est_mu=1:num_mu_est,est_pi=1:num_pi_est)
t=1;k=1
if (!is.null(getOption("cl"))) { stopCluster(getOption("cl")) }
cl <- makeCluster(limit)  # 使用核心数减 1，保留一个核心给操作系统
registerDoParallel(cl)
replication_ids <- 1:100
results_par=foreach (t=replication_ids, .combine=rbind, .packages=c("ranger", "MASS", "matrixStats","glmnet","RCAL", "mvtnorm", "caret", "densratio","randomForestSRC","butcher","logger","filelock","survival", "dplyr","tidyr","purrr","ks","corpcor")) %dorng% {
  # for (t in 1:tmax){
  set.seed(seed+t)
  
  index_out=(t-1)*N+1:N
  if (t>1) {index_out=sample(N, N, replace = FALSE)}
  Y=Y_all[index_out]; X=X_all[index_out,]; M=M_all[index_out,]; S=S_all[index_out,]; A=A_all[index_out]; X_poly=X_expanded[index_out,]; S_poly=S_expanded[index_out,]
  index_out_sdr=index_out
  index_out=1:N
  n=N/Kmb
  
  X_poly=scale(X_poly); S_poly=scale(S_poly)
  # Existing mediator-containing columns, including precomputed M-X interactions.
  # Pure X columns remain in X_poly; no new powers or products are constructed.
  S_poly_log=scale(S_sdr_log_expanded[index_out_sdr,,drop=FALSE])
  M_poly=S_poly_log[,M_vars,drop=FALSE]
  
  d1_poly=ncol(X_poly)
  d2=ncol(M)
  d3_poly=ncol(S_poly)
  
  A_fac=factor(A, levels=c(0,1), labels=c("Control", "Treatment"))
  
  sel.rho_pi_a=c();sel.rho_pi_b=c();sel.rho_q=c(); sel.rho_p=c();lambda_mu=c();lambda_tau=c();lambda_tau1=c();
  
  mtry_pi=c();mtry_p=c();mtry_mu=c();mtry_tau1=c();mtry_tau=c();
  minsize_pi=c();minsize_p=c(); minsize_mu=c(); minsize_tau1=c();minsize_tau=c();
  
  pred_pi_a=matrix(0,nrow=num_pi_est,ncol=N)
  pred_pi_b=matrix(0,nrow=num_pi_est,ncol=N)
  pred_q=matrix(0,nrow=num_q_est,ncol=N)
  pred_p=matrix(0,nrow=num_p_est,ncol=N)
  pred_f1=matrix(0,nrow=1,ncol=N)
  pred_f0=matrix(0,nrow=1,ncol=N)
  pred_mu=matrix(0,nrow=num_mu_est,ncol=N)
  pred_tau=matrix(0,nrow=num_mu_est,ncol=N)
  pred_tau1=matrix(0,nrow=num_mu_est,ncol=N)
  
  pred_theta=matrix(0,nrow=num_mu_est,ncol=N)
  pred_11theta=matrix(0,nrow=num_mu_est,ncol=N)
  pred_00theta=matrix(0,nrow=num_mu_est,ncol=N)
  
  pred_pi_00=matrix(0,nrow=2,ncol=N)
  pred_mu_00=matrix(0,nrow=2,ncol=N)
  pred_mu_11=matrix(0,nrow=2,ncol=N)
  
  log_appender(appender_file(log_file))
  log_formatter(formatter_paste)
  log_threshold(INFO)
  
  set.seed(seed+t)
    
  #------------MQR1-----------
  for (k in 1:Kmb){
    index=(1:n)+(k-1)*n
    index_nk=(1:N)[-index]
    index_nk1=index_nk[1:floor(length(index_nk)/2)]
    # index_nk2=index_nk[-(1:floor(length(index_nk)/2))]
    #
    # index_train1=index_nk1[1:floor(length(index_nk1)/2)]
    # index_train2=index_nk1[-(1:floor(length(index_nk1)/2))]
    # index_train3=index_nk2[1:floor(length(index_nk2)/2)]
    # index_train4=index_nk2[-(1:floor(length(index_nk2)/2))]
    # index_train5=index_train4[A[index_train4]==1]
    # index_train4=index_train4[A[index_train4]==0]
  
    index_train1=index_nk[1:floor(length(index_nk)/5)]
    index_nk_2345=index_nk[-(1:floor(length(index_nk)/5))]
    index_train23=index_nk_2345[1:floor(length(index_nk_2345)/2)]
    index_train45=index_nk_2345[-(1:floor(length(index_nk_2345)/2))]
    index_train2=index_train23[1:floor(length(index_train23)/2)]
    index_train3=index_train23[-(1:floor(length(index_train23)/2))]
    index_train4=index_train45[1:floor(length(index_train45)/2)]
    index_train5=index_train45[-(1:floor(length(index_train45)/2))]
    rm(index_nk_2345,index_train23, index_train45)
  
    #pi_a
    # sel.rho_pi_a=c()
    rho.range0=sqrt(log(d1_poly)/length(index_train1))*rho.range0_default
    boolFalse=FALSE #repeat until no errors occur
    while(boolFalse==FALSE){
      tryCatch({
        fit_Log1_1=glm.regu.cv(x=X_poly[index_train1,],y=A[index_train1],fold=5,loss="cal",rho.seq=rho.range0, ann.init=FALSE,tune.cut=FALSE)
        boolFalse=TRUE
      },error=function(e){
      },finally={})
    }
    sel.rho_pi_a=c(sel.rho_pi_a,fit_Log1_1$sel.rho[sw1]/sqrt(log(d1_poly)/length(index_train1)))
    alpha_a=fit_Log1_1$sel.bet[,sw1]
    pred_pi_a[1,index]=pmax(pmin(plogis(alpha_a[1]+X_poly[index,]%*%alpha_a[-1]),1),0)
    # if(max(pred_pi_a[1,index])>0.95 || min(pred_pi_a[1,index])<0.05) {
    #   alpha_a=fit_Log1_1$sel.bet[,2]
    #   pred_pi_a[1,index]=pmax(pmin(plogis(alpha_a[1]+X_poly[index,]%*%alpha_a[-1]),1),0)
    # }
    # mean((abs(pred_pi_a[1,index]-pi[index])/pi[index]))
  
    #pi_b
    # sel.rho_pi_b=c()
    rho.range0=sqrt(log(d1_poly)/length(index_train1))*rho.range0b_default
    boolFalse=FALSE #repeat until no errors occur
    while(boolFalse==FALSE){
      tryCatch({
        fit_Log1_1=glm.regu.cv(x=-X_poly[index_train1,],y=1-A[index_train1],fold=5,loss="cal",rho.seq=rho.range0, ann.init=FALSE,tune.cut=FALSE)
        boolFalse=TRUE
      },error=function(e){
      },finally={})
    }
    sel.rho_pi_b=c(sel.rho_pi_b,fit_Log1_1$sel.rho[sw1]/sqrt(log(d1_poly)/length(index_train2)))
    alpha_b=fit_Log1_1$sel.bet[,sw1]
    alpha_b[1]=-alpha_b[1]
    pred_pi_b[1,index]=pmax(pmin(plogis(alpha_b[1]+X_poly[index,]%*%alpha_b[-1]),1),0)
    # if(max(pred_pi_b[1,index])>0.95 || min(pred_pi_b[1,index])<0.05){
    #   alpha_b=fit_Log1_1$sel.bet[,2]
    #   alpha_b[1]=-alpha_b[1]
    #   pred_pi_b[1,index]=pmax(pmin(plogis(alpha_b[1]+X_poly[index,]%*%alpha_b[-1]),1),0)
    # }
  
    #q
    # sel.rho_q=c()
    col_variances <- apply(S_poly[index_train2,], 2, var, na.rm = TRUE)
    zero_var_cols <- which(col_variances < .Machine$double.eps)
    S_index2_filtered <- S_poly[index_train2,col_variances > .Machine$double.eps, drop = FALSE]
    rho.range1=sqrt(log(d3_poly)/length(index_train2))*rho.range1_default
    boolFalse=FALSE #repeat until no errors occur
    while(boolFalse==FALSE){
      tryCatch({
        fit_Log2_1=glm.regu.cv(x=-S_index2_filtered, y=A[index_train2], iw=as.numeric(A[index_train2]*(1+exp(-alpha_a[1]-X_poly[index_train2,]%*%alpha_a[-1]))+(1-A[index_train2])*(1+exp(alpha_b[1]+X_poly[index_train2,]%*%alpha_b[-1])))
                               , fold=5,loss="cal",rho.seq=rho.range1, ann.init=FALSE,tune.cut=FALSE)
        boolFalse=TRUE
      },error=function(e){
      },finally={})
    }
    sel.rho_q=c(sel.rho_q,fit_Log2_1$sel.rho[sw1]/sqrt(log(d3_poly)/length(index_train2)))
    beta_q1=fit_Log2_1$sel.bet[,sw1]
    beta_q1[1]=-beta_q1[1]
    pred_q[1,index]=exp(beta_q1[1]+S_poly[index,col_variances > .Machine$double.eps, drop = FALSE]%*%beta_q1[-1])
    # if(max(pred_q[1,index])>20 || min(pred_q[1,index])<0.05){
    #   beta_q1=fit_Log2_1$sel.bet[,2]
    #   beta_q1[1]=-beta_q1[1]
    #   pred_q[1,index]=exp(beta_q1[1]+S_poly[index,col_variances > .Machine$double.eps, drop = FALSE]%*%beta_q1[-1])
    # }
  
    #mu
    # lambda_mu=c()
    fit_lasso_mu=cv.glmnet(x=S_poly[index_train3,],y=Y[index_train3],weights=A[index_train3]*exp(beta_q1[1]+S_poly[index_train3,col_variances > .Machine$double.eps, drop = FALSE]%*%beta_q1[-1])*(1+exp(-alpha_a[1]-X_poly[index_train3,]%*%alpha_a[-1]))#/sum(A[index_train3]*exp(beta_q1[1]+S_poly[index_train3,col_variances > .Machine$double.eps, drop = FALSE]%*%beta_q1[-1])*(1+exp(-alpha_a[1]-X_poly[index_train3,]%*%alpha_a[-1])))
                           ,family="gaussian",alpha=1,standardize=FALSE,parallel=FALSE,nfolds=5)
    pred_mu[1,index]=predict(fit_lasso_mu,newx=S_poly[index,],s="lambda.min",standardize=FALSE)
    lambda_mu=c(lambda_mu,fit_lasso_mu$lambda.min/sqrt(log(d3_poly)/length(index_train3)))
  
    #tau1
    index_train4_A0=index_train4[which(A[index_train4]==0)]
    find_mu=which(fit_lasso_mu$lambda==fit_lasso_mu$lambda.min)
    a0_mu=fit_lasso_mu$glmnet.fit$a0[find_mu]
    alpha_mu_X=fit_lasso_mu$glmnet.fit$beta[1:d1_poly,find_mu]
    alpha_mu_M=fit_lasso_mu$glmnet.fit$beta[(d1_poly+1):nrow(fit_lasso_mu$glmnet.fit$beta),find_mu]
    pred_mu_new1=S_poly[index_train4,(d1_poly+1):(d3_poly)]%*%alpha_mu_M+a0_mu
    #pred_mu_new1=predict(fit_lasso_mu,newx=S_poly[index_train4,],s="lambda.min")
    if (sd(pred_mu_new1)<=1e-5 || length(pred_mu_new1)<=8){pred_tau1[1,index]=X_poly[index,]%*%alpha_mu_X+
      rep(mean(pred_mu_new1),length(index)); lambda_tau1=c(lambda_tau1,-1)} else{
        fit_lasso_tau1=tryCatch({fit=cv.glmnet(x=X_poly[index_train4,],y=pred_mu_new1, weights=((1-A[index_train4])*exp(alpha_b[1]+X_poly[index_train4,]%*%alpha_b[-1]))#/sum((1-A[index_train4])*exp(alpha_b[1]+X_poly[index_train4,]%*%alpha_b[-1]))
                                               ,family="gaussian",alpha=1,standardize=FALSE,parallel=FALSE,nfolds=5)
        list(fit = fit, error = NULL)
        }, error = function(e) {
          list(fit = NULL, error = e)
        })
        if (!is.null(fit_lasso_tau1$error)) {
          pred_tau1[1, index]=X_poly[index,]%*%alpha_mu_X +
            rep(mean(pred_mu_new1), length(index))
          lambda_tau1 = c(lambda_tau1, -1)
        } else {
          pred_tau1[1,index]=X_poly[index,]%*%alpha_mu_X+predict(fit_lasso_tau1$fit,newx=X_poly[index,],s="lambda.min")
          lambda_tau1=c(lambda_tau1,fit_lasso_tau1$fit$lambda.min/sqrt(log(d1_poly)/length(index_train4)))
        }
      }
  
    #tau
    pred_mu_new2=S_poly[index_train5,(d1_poly+1):(d3_poly)]%*%alpha_mu_M+a0_mu
    #pred_mu_new2=predict(fit_lasso_mu,newx=S_poly[index_train5,], s="lambda.min")
    if (sd(pred_mu_new1)<=1e-5 || length(pred_mu_new1)<=8 ){pred_tau1_new1=X_poly[index_train5,]%*%alpha_mu_X+
      rep(mean(pred_mu_new2),length(index_train5))} else{
        pred_tau1_new1=X_poly[index_train5,]%*%alpha_mu_X+
          predict(fit_lasso_tau1$fit,newx=X_poly[index_train5,],s="lambda.min",standardize=FALSE)
      }
    y_pseudo=pred_tau1_new1+exp(beta_q1[1]+S_poly[index_train5,col_variances > .Machine$double.eps, drop = FALSE]%*%beta_q1[-1])*(Y[index_train5]-pred_mu_new2-X_poly[index_train5,]%*%alpha_mu_X)
    if (sd(y_pseudo)<=1e-5 || length(y_pseudo)<=8 ){pred_tau[1,index]=rep(mean(y_pseudo),length(index));lambda_tau=c(lambda_tau,0)} else{
      fit_lasso_tau=cv.glmnet(x=X_poly[index_train5,],y=y_pseudo, weights=A[index_train5]*exp(-alpha_a[1]-X_poly[index_train5,]%*%alpha_a[-1])#/sum(A[index_train5]*exp(-alpha_a[1]-X_poly[index_train5,]%*%alpha_a[-1]))
                              ,family="gaussian",alpha=1,standardize=FALSE,parallel=FALSE,nfolds=5)
      pred_tau[1,index]=predict(fit_lasso_tau,newx=X_poly[index,],s="lambda.min",standardize=FALSE)
      lambda_tau=c(lambda_tau,fit_lasso_tau$lambda.min/sqrt(log(d1_poly)/length(index_train4)))}
  
  }
  lock_obj=lock(log_file, timeout=5)
  unlock(lock_obj)
  log_info(paste("Time=", t  ,"MQR1 accomplished"))
  set.seed(seed+t)
  # 
  #------------MQR2------------------
  for (k in 1:Kmb){
    index=(1:n)+(k-1)*n
    index_nk=(1:N)[-index]
    index_nk_A=index_nk[A[index_nk]==1]
    index_nk_A0=index_nk[A[index_nk]==0]

    #pi_a
    rho.range0=sqrt(log(d1_poly)/length(index_nk))*rho.range02_default
    boolFalse=FALSE #repeat until no errors occur
    while(boolFalse==FALSE){
      tryCatch({
        fit_Log1_1=glm.regu.cv(x=X_poly[index_nk,],y=A[index_nk],fold=5,loss="cal",rho.seq=rho.range0, ann.init=FALSE,tune.cut=FALSE)
        boolFalse=TRUE
      },error=function(e){
      },finally={})
    }
    sel.rho_pi_a=c(sel.rho_pi_a,fit_Log1_1$sel.rho[sw2]/sqrt(log(d1_poly)/length(index_nk)))
    alpha_a=fit_Log1_1$sel.bet[,sw2]
    pred_pi_a[2,index]=pmax(pmin(plogis(alpha_a[1]+X_poly[index,]%*%alpha_a[-1]),1),0) #pi_a(X)
    # if(max(pred_pi_a[2,index])>0.95 || min(pred_pi_a[2,index])<0.05){
    #   alpha_a=fit_Log1_1$sel.bet[,2]
    #   pred_pi_a[2,index]=pmax(pmin(plogis(alpha_a[1]+X_poly[index,]%*%alpha_a[-1]),1),0)
    # }

    #pi_b
    rho.range0=sqrt(log(d1_poly)/length(index_nk))*rho.range0b_default
    boolFalse=FALSE #repeat until no errors occur
    while(boolFalse==FALSE){
      tryCatch({
        fit_Log1_1=glm.regu.cv(x=-X_poly[index_nk,],y=1-A[index_nk],fold=5,loss="cal",rho.seq=rho.range0, ann.init=FALSE,tune.cut=FALSE)
        boolFalse=TRUE
      },error=function(e){
      },finally={})
    }
    sel.rho_pi_b=c(sel.rho_pi_b,fit_Log1_1$sel.rho[sw2]/sqrt(log(d1_poly)/length(index_nk)))
    alpha_b=fit_Log1_1$sel.bet[,sw2]
    alpha_b[1]=-alpha_b[1]
    pred_pi_b[2,index]=pmax(pmin(plogis(alpha_b[1]+X_poly[index,]%*%alpha_b[-1]),1),0) #pi_b(X)
    # if(max(pred_pi_b[2,index])>0.95 || min(pred_pi_b[2,index])<0.05){
    #   alpha_b=fit_Log1_1$sel.bet[,2]
    #   alpha_b[1]=-alpha_b[1]
    #   pred_pi_b[2,index]=pmax(pmin(plogis(alpha_b[1]+X_poly[index,]%*%alpha_b[-1]),1),0)
    # }

    #q
    S_index2=S_poly[index_nk,]
    rho.range1=sqrt(log(d3_poly)/length(index_nk))*rho.range12_default
    col_variances <- apply(S_index2, 2, var, na.rm = TRUE)
    zero_var_cols <- which(col_variances < .Machine$double.eps)
    S_index2_filtered <- S_index2[, col_variances > .Machine$double.eps, drop = FALSE]
    boolFalse=FALSE #repeat until no errors occur
    while(boolFalse==FALSE){
      tryCatch({
        fit_Log2_1=glm.regu.cv(x=-S_index2_filtered, y=A[index_nk], iw=as.numeric(A[index_nk]/plogis(alpha_a[1]+X_poly[index_nk,]%*%alpha_a[-1])+(1-A[index_nk])/(1-plogis(alpha_b[1]+X_poly[index_nk,]%*%alpha_b[-1]))),fold=5,loss="cal",rho.seq=rho.range1, ann.init=FALSE,tune.cut=FALSE)
        boolFalse=TRUE
      },error=function(e){
      },finally={})
    }
    sel.rho_q=c(sel.rho_q,fit_Log2_1$sel.rho[sw2]/sqrt(log(d3_poly)/length(index_nk)))
    beta_q=fit_Log2_1$sel.bet[,sw2]
    beta_q[1]=-beta_q[1]
    pred_q[2,index]=exp(beta_q[1]+S_poly[index,col_variances > .Machine$double.eps, drop = FALSE]%*%beta_q[-1])  #q(S)
    # if(max(pred_q[2,index])>20 || min(pred_q[2,index])<0.05){
    #   beta_q=fit_Log2_1$sel.bet[,2]
    #   beta_q[1]=-beta_q[1]
    #   pred_q[2,index]=exp(beta_q[1]+S_poly[index,col_variances > .Machine$double.eps, drop = FALSE]%*%beta_q[-1])
    # }

    #mu
    fit_lasso_mu=cv.glmnet(x=S_poly[index_nk,],y=Y[index_nk],weights=A[index_nk]*exp(beta_q[1]+S_poly[index_nk,col_variances > .Machine$double.eps, drop = FALSE]%*%beta_q[-1])*(1+exp(-alpha_a[1]-X_poly[index_nk,]%*%alpha_a[-1])),nfolds=10, family="gaussian",alpha=1,standardize=FALSE,parallel=FALSE)
    pred_mu[2,index]=predict(fit_lasso_mu,newx=S_poly[index,],s="lambda.min",standardize=FALSE)
    lambda_mu=c(lambda_mu,fit_lasso_mu$lambda.min/sqrt(log(d3_poly)/length(index_nk)))

    #tau1
    find_mu=which(fit_lasso_mu$lambda==fit_lasso_mu$lambda.min)
    a0_mu=fit_lasso_mu$glmnet.fit$a0[find_mu]
    alpha_mu_X=fit_lasso_mu$glmnet.fit$beta[1:d1_poly,find_mu]
    alpha_mu_M=fit_lasso_mu$glmnet.fit$beta[(d1_poly+1):(d3_poly),find_mu]
    pred_mu_new1=a0_mu+S_poly[index_nk,(d1_poly+1):(d3_poly)]%*%alpha_mu_M
    #pred_mu_new1=predict(fit_lasso_mu,newx=S_poly[index_nk,],s="lambda.min")
    if (sd(pred_mu_new1)<=1e-5){pred_tau1[2,index]=X_poly[index,]%*%alpha_mu_X+
      rep(mean(pred_mu_new1),length(index));lambda_tau1=c(lambda_tau1,-1)} else{
        fit_lasso_tau1=cv.glmnet(x=X_poly[index_nk,],y=pred_mu_new1, weights=(1-A[index_nk])*exp(alpha_b[1]+X_poly[index_nk,]%*%alpha_b[-1]),nfolds=10
                                 ,family="gaussian",alpha=1,standardize=FALSE,parallel=FALSE)
        pred_tau1[2,index]=X_poly[index,]%*%alpha_mu_X+
          predict(fit_lasso_tau1,newx=X_poly[index,],s="lambda.min",standardize=FALSE)
        lambda_tau1=c(lambda_tau1,fit_lasso_tau1$lambda.min/sqrt(log(d1_poly)/length(index_nk)))
      }

    #tau
    pred_mu_new1=a0_mu+S_poly[index_nk,(d1_poly+1):(d3_poly)]%*%alpha_mu_M
    if (sd(pred_mu_new1)<=1e-5){pred_tau1_new1=X_poly[index_nk,]%*%alpha_mu_X+
      rep(mean(pred_mu_new1),length(index_nk))} else{
        pred_tau1_new1=X_poly[index_nk,]%*%alpha_mu_X+
          predict(fit_lasso_tau1,newx=X_poly[index_nk,],s="lambda.min",standardize=FALSE)
      }
    y_pseudo=pred_tau1_new1+exp(beta_q[1]+S_poly[index_nk,col_variances > .Machine$double.eps, drop = FALSE]%*%beta_q[-1])*(Y[index_nk]-pred_mu_new1-X_poly[index_nk,]%*%alpha_mu_X)
    if (sd(y_pseudo)<=1e-5){pred_tau[2,index]=rep(mean(y_pseudo),length(index));lambda_tau=c(lambda_tau,-1)} else{
      fit_lasso_tau=cv.glmnet(x=X_poly[index_nk,],y=y_pseudo, weights=exp(-alpha_a[1]-X_poly[index_nk,]%*%alpha_a[-1]),nfolds=10, family="gaussian",alpha=1,standardize=FALSE,parallel=FALSE)
      pred_tau[2,index]=predict(fit_lasso_tau,newx=X_poly[index,],s="lambda.min",standardize=FALSE)
      lambda_tau=c(lambda_tau,fit_lasso_tau$lambda.min/sqrt(log(d1_poly)/length(index_nk)))
    }

  }
  lock_obj=lock(log_file, timeout=5)
  log_info(paste("Time=", t  ,"MQR2 accomplished"))
  unlock(lock_obj)
  set.seed(seed+t)

  #----------------QR1-----------
  for (k in 1:Kmb) {
    index=(1:n)+(k-1)*n
    index_nk=(1:N)[-index]
    index_nk3=index_nk[floor(length(index_nk)/3*2):length(index_nk)]
    index_nk12=index_nk[-(floor(length(index_nk)/3*2):length(index_nk))]
    index_nk12_A=index_nk12[A[index_nk12]==1]
    index_nk12_A0=index_nk12[A[index_nk12]==0]
    index_nk1=index_nk12[1:floor(length(index_nk12)/2)]
    index_nk2=index_nk12[-(1:floor(length(index_nk12)/2))]
    index_nk1_A=index_nk1[A[index_nk1]==1]
    index_nk2_A0=index_nk2[A[index_nk2]==0]
    index_nk3_A=index_nk3[A[index_nk3]==1]
  
    #pi using Random forest
    lasso_pi <- lasso_var_selection(X_poly[index_nk12,],A[index_nk12],task=1)
    selected_vars_pi <- lasso_pi$selected_vars
    pi_rf=caret::train(x=X_poly[index_nk12,selected_vars_pi,drop=FALSE], y=A_fac[index_nk12], method= "ranger", trControl=ctrl_pi, tuneGrid=generate_tune_grid(length(index_nk12),length(selected_vars_pi),1), metric="ROC", importance="impurity", num.trees=generate_num_tree(length(index_nk12)))#,class.weights = c("Control" = 0.5, "Treatment" = 0.5)
    prediction_pi=predict(pi_rf$finalModel, data=X_poly[index,selected_vars_pi,drop=FALSE])
    pred_pi_a[3,index]=pmax(pmin(prediction_pi$predictions[,2], 0.99), 0.01)
    pred_pi_b[3,index]=pred_pi_a[3,index]
    mtry_pi=c(mtry_pi,pi_rf$bestTune$mtry)
    minsize_pi=c(minsize_pi,pi_rf$bestTune$min.node.size)
  
    # f semiparametrically
    tm_mediator_screen <- proc.time()[["elapsed"]]
    selected_M <- select_different_M_components_ks(
      X_poly[index_nk12, ],
      M_poly[index_nk12, ],
      A[index_nk12],
      max_selected = 8,
      min_selected = 1
    )
    mediator_screen_sec <- proc.time()[["elapsed"]] - tm_mediator_screen
  
    x_vars <- multi_y_var_selection(
      X_poly[index_nk12, ],
      M_poly[index_nk12, selected_M, drop = FALSE]
    )$union_vars
    q_basis_poly <- sdr_existing_mediator_basis(
      S_poly_log, M_vars, colnames(M_poly)[selected_M]
    )
    tm_basis_selection <- proc.time()[["elapsed"]]
    basis_sel <- select_el_drm_interactions(
      X = X_poly[index_nk12, , drop = FALSE],
      initial_x_vars = x_vars,
      M = q_basis_poly[index_nk12, , drop = FALSE],
      linear_names = colnames(M_poly)[selected_M],
      A = A[index_nk12],
      stability_reps = basis_stability_reps,
      subsample_frac = 0.8,
      cv_folds = 5,
      alpha = 0.5,
      stability_cutoff = 0.8,
      seed = seed + 20000 * t + 200 * k
    )
    basis_selection_sec <- proc.time()[["elapsed"]] - tm_basis_selection
  
    selected_basis <- basis_sel$selected_basis
    x_vars <- basis_sel$selected_x
  
    log_info(paste(
      "QR1 basis fold", k,
      "; selected_M =", paste(selected_M, collapse = ","),
      "; selected_basis =", paste(selected_basis, collapse = ","),
      "; successful_reps =", basis_sel$successful_reps,
      "; stable_interactions =", nrow(basis_sel$selected_interactions),
      "; selected_X =", paste(x_vars, collapse = ","),
      "; mediator_screen_sec =", round(mediator_screen_sec, 3),
      "; basis_selection_sec =", round(basis_selection_sec, 3),
      "; quadratic_freq =",
      paste(
        names(basis_sel$quadratic_frequency),
        round(basis_sel$quadratic_frequency, 2),
        sep = ":",
        collapse = ","
      )
    ))
  
    tm_el_drm <- proc.time()[["elapsed"]]
    f_el_drm_sdr1 <- el_drm_fit_density_ratio(
      X_poly[index_nk12, x_vars, drop = FALSE],
      q_basis_poly[index_nk12, , drop = FALSE],
      A[index_nk12],
      basis_names = selected_basis,
      ridge_grid = c(
        0, 0.01, 0.03, 0.05, 0.1,
        0.2, 0.3, 0.5, 1
      ),
      ridge_cv_folds = 5,
      ridge_cv_seed = seed + 10000 * t + 100 * k + 2
    )
    el_drm_sec <- proc.time()[["elapsed"]] - tm_el_drm
    log_ratio_index <- el_drm_predict_log_ratio(
      f_el_drm_sdr1,
      X_poly[index, x_vars, drop = FALSE],
      q_basis_poly[index, , drop = FALSE]
    )
  
    ref_log_q_sdr1 <- f_el_drm_sdr1$full_train_log_ratio[
      A[index_nk12] == 1
    ]
  
    pred_q[3,index] <- normalized_truncate_q(
      log_q = log_ratio_index,
      ref_log_q = ref_log_q_sdr1
    )
    # pred_q[3,index]=exp(pmax(pmin(log_ratio_index, log(100)), log(0.01)))
    # pred_q[3,index]=exp(log_ratio_index)
  
    log_info(paste(
      "QR1 fold", k,
      "ridge_lambda =", f_el_drm_sdr1$ridge_lambda,
      "x_dim =", f_el_drm_sdr1$x_dim,
      "q_dim =", f_el_drm_sdr1$q_dim,
      "theta_dim =", length(f_el_drm_sdr1$theta),
      "el_drm_sec =", round(el_drm_sec, 3)
    ))
  
    #mu
    lasso_mu<-lasso_var_selection(S_poly_log[index_nk1_A,,drop=FALSE], Y[index_nk1_A], task = 0)
    selected_vars_mu<-lasso_mu$selected_vars
    mu_rf=caret::train(x=S_poly_log[index_nk1_A,,drop=FALSE][,selected_vars_mu,drop=FALSE], y=Y[index_nk1_A], method= "ranger", tuneGrid=generate_tune_grid(length(index_nk1_A),length(selected_vars_mu),0), trControl=ctrl_mu, importance="impurity", num.trees=generate_num_tree(length(index_nk1_A)) )
    pred_mu[3,index]=predict(mu_rf$finalModel, data=S_poly_log[index,,drop=FALSE][,selected_vars_mu,drop=FALSE])$predictions
    mtry_mu=c(mtry_mu,mu_rf$bestTune$mtry)
    minsize_mu=c(minsize_mu,mu_rf$bestTune$min.node.size)
    # summary(pred_mu[3,index]-mu[index])
  
    #tau1
    prediction_mu_1=predict(mu_rf$finalModel, data=S_poly_log[index_nk2_A0,,drop=FALSE][,selected_vars_mu,drop=FALSE])
    lasso_tau1<-lasso_var_selection(X_poly[index_nk2_A0,], prediction_mu_1$predictions, task = 0)
    selected_vars_tau1<-lasso_tau1$selected_vars
    tau1_rf=caret::train(x=X_poly[index_nk2_A0,selected_vars_tau1,drop=FALSE], y=prediction_mu_1$predictions, method= "ranger", tuneGrid=generate_tune_grid(length(index_nk2_A0),length(selected_vars_tau1),0), trControl=ctrl_mu, importance="impurity", num.trees=generate_num_tree(length(index_nk2_A0)))
    mtry_tau1=c(mtry_tau1,tau1_rf$bestTune$mtry)
    minsize_tau1=c(minsize_tau1,tau1_rf$bestTune$min.node.size)
    pred_tau1[3,index]=predict(tau1_rf$finalModel, data=X_poly[index,selected_vars_tau1,drop=FALSE])$predictions
    # summary(pred_tau1[3,index])
  
    #tau
    # prediction_pi_1=pmax(pmin(predict(pi_rf$finalModel, data=X_poly[index_nk3_A,selected_vars_pi,drop=FALSE])$predictions[, 2], 0.9999), 0.0001)
    # prediction_p_1=pmax(pmin(predict(p_rf$finalModel, data=S_poly_log[index_nk3_A,,drop=FALSE][,selected_vars_p,drop=FALSE])$predictions[, 2], 0.9999), 0.0001)
    # prediction_q_1=(1-prediction_p_1)/prediction_p_1*prediction_pi_1/(1-prediction_pi_1)
  
    log_ratio_index_nk3_A <- el_drm_predict_log_ratio(f_el_drm_sdr1,
                                                      X_poly[index_nk3_A, x_vars, drop = FALSE],
                                                      q_basis_poly[index_nk3_A, , drop = FALSE])
    prediction_q_1 <- normalized_truncate_q(
      log_q = log_ratio_index_nk3_A,
      ref_log_q = ref_log_q_sdr1
    )
    # prediction_q_1=exp(pmax(pmin(log_ratio_index_nk3_A, log(100)), log(0.01)))
    # prediction_q_1=exp(log_ratio_index_nk3_A)
  
    prediction_mu_2=predict(mu_rf$finalModel, data=S_poly_log[index_nk3_A,,drop=FALSE][,selected_vars_mu,drop=FALSE])
    prediction_tau1_1=predict(tau1_rf$finalModel, data=X_poly[index_nk3_A,selected_vars_tau1,drop=FALSE])
    pred_pseudoY=prediction_q_1*(Y[index_nk3_A]-prediction_mu_2$predictions)+prediction_tau1_1$predictions
    lasso_tau<-lasso_var_selection(X_poly[index_nk3_A,], pred_pseudoY, task = 0)
    selected_vars_tau<-lasso_tau$selected_vars
    tau_rf=caret::train(x=X_poly[index_nk3_A,selected_vars_tau,drop=FALSE], y=pred_pseudoY, method= "ranger", tuneGrid=generate_tune_grid(length(index_nk3_A),length(selected_vars_tau),0), trControl=ctrl_mu, importance="impurity", num.trees=generate_num_tree(length(index_nk3_A)) )
    mtry_tau=c(mtry_tau,tau_rf$bestTune$mtry)
    minsize_tau=c(minsize_tau,tau_rf$bestTune$min.node.size)
  
    pred_tau[3,index]=predict(tau_rf$finalModel, data=X_poly[index,selected_vars_tau,drop=FALSE])$predictions
    # p_rf=butcher(p_rf); rm(p_rf);rm(prediction_pi_1,prediction_p_1);
    #q1_crulsif=butcher(q1_crulsif);q2_crulsif=butcher(q2_crulsif);rm(q1_crulsif,q2_crulsif)
    pi_rf=butcher(pi_rf); rm(pi_rf,prediction_pi)
    mu_rf=butcher(mu_rf);tau1_rf=butcher(tau1_rf);tau_rf=butcher(tau_rf); rm(mu_rf,tau1_rf,tau_rf,prediction_mu_1,prediction_mu_2,prediction_tau1_1)
    rm(prediction_q_1)
    gc()
    # pred_theta[3,index]=A[index]*pred_q[3,index]/(pred_pi_a[3,index])*(Y[index]-pred_mu[3,index])+(1-A[index])/(1-pred_pi_b[3,index])*(pred_mu[3,index]-pred_tau1[3,index])+A[index]/(pred_pi_a[3,index])*(pred_tau1[3,index]-pred_tau[3,index])+pred_tau[3,index];summary(pred_theta[3,index]);
  }
  lock_obj=lock(log_file, timeout=5)
  unlock(lock_obj)
  log_info(paste("Time=", t, "QR1 accomplished"))
  set.seed(seed+t)

  # --------------------QR2---------
  for (k in 1:Kmb) {
    index=(1:n)+(k-1)*n
    index_nk=(1:N)[-index]
    index_nk_A=index_nk[A[index_nk]==1]
    index_nk_A0=index_nk[A[index_nk]==0]

    #pi using Random forest
    lasso_pi <- lasso_var_selection(X_poly[index_nk,],A[index_nk],task=1)
    selected_vars_pi <- lasso_pi$selected_vars
    pi_rf=caret::train(x=X_poly[index_nk,selected_vars_pi,drop=FALSE], y=A_fac[index_nk], method= "ranger", trControl=ctrl_pi, tuneGrid=generate_tune_grid(length(index_nk),length(selected_vars_pi),1), metric="ROC", importance="impurity", num.trees=generate_num_tree(length(index_nk)))
    prediction_pi=predict(pi_rf$finalModel, data=X_poly[index,selected_vars_pi,drop=FALSE])
    pred_pi_a[4,index]=pmax(pmin(prediction_pi$predictions[, 2],0.99), 0.01)
    pred_pi_b[4,index]=pred_pi_a[4,index]
    mtry_pi=c(mtry_pi,pi_rf$bestTune$mtry)
    minsize_pi=c(minsize_pi,pi_rf$bestTune$min.node.size)

    # q using Random forest
    tm_mediator_screen <- proc.time()[["elapsed"]]
    selected_M <- select_different_M_components_ks(
      X_poly[index_nk, ],
      M_poly[index_nk, ],
      A[index_nk],
      max_selected = 8,
      min_selected = 1
    )
    mediator_screen_sec <- proc.time()[["elapsed"]] - tm_mediator_screen
    x_vars <- multi_y_var_selection(X_poly[index_nk, ], M_poly[index_nk, selected_M, drop = FALSE])$union_vars
    q_basis_poly <- sdr_existing_mediator_basis(
      S_poly_log, M_vars, colnames(M_poly)[selected_M]
    )
    tm_basis_selection <- proc.time()[["elapsed"]]
    basis_sel <- select_el_drm_interactions(
      X = X_poly[index_nk, , drop = FALSE],
      initial_x_vars = x_vars,
      M = q_basis_poly[index_nk, , drop = FALSE],
      linear_names = colnames(M_poly)[selected_M],
      A = A[index_nk],
      stability_reps = basis_stability_reps,
      subsample_frac = 0.8,
      cv_folds = 5,
      alpha = 0.5,
      stability_cutoff = 0.8,
      seed = seed + 20000 * t + 200 * k
    )
    basis_selection_sec <- proc.time()[["elapsed"]] - tm_basis_selection

    selected_basis <- basis_sel$selected_basis
    x_vars <- basis_sel$selected_x

    log_info(paste(
      "QR2 basis fold", k,
      "; selected_M =", paste(selected_M, collapse = ","),
      "; selected_basis =", paste(selected_basis, collapse = ","),
      "; successful_reps =", basis_sel$successful_reps,
      "; stable_interactions =", nrow(basis_sel$selected_interactions),
      "; selected_X =", paste(x_vars, collapse = ","),
      "; mediator_screen_sec =", round(mediator_screen_sec, 3),
      "; basis_selection_sec =", round(basis_selection_sec, 3),
      "; quadratic_freq =",
      paste(
        names(basis_sel$quadratic_frequency),
        round(basis_sel$quadratic_frequency, 2),
        sep = ":",
        collapse = ","
      )
    ))

    tm_el_drm <- proc.time()[["elapsed"]]
    f_el_drm <- el_drm_fit_density_ratio(
      X_poly[index_nk, x_vars, drop = FALSE],
      q_basis_poly[index_nk, , drop = FALSE],
      A[index_nk],
      basis_names = selected_basis,
      ridge_grid = c(
        0, 0.01, 0.03, 0.05, 0.1,
        0.2, 0.3, 0.5, 1
      ),
      ridge_cv_folds = 5,
      ridge_cv_seed = seed + 10000 * t + 100 * k + 2
    )
    el_drm_sec <- proc.time()[["elapsed"]] - tm_el_drm
    log_ratio_index <- el_drm_predict_log_ratio(f_el_drm,
                                                X_poly[index, x_vars, drop = FALSE],
                                                q_basis_poly[index, , drop = FALSE])
    ref_log_q_sdr2_full <- f_el_drm$full_train_log_ratio[
      A[index_nk] == 1
    ]

    pred_q[4,index] <- normalized_truncate_q(
      log_q = log_ratio_index,
      ref_log_q = ref_log_q_sdr2_full
    )
    # pred_q[4,index]=exp(pmax(pmin(log_ratio_index, log(100)), log(0.01)))
    # pred_q[4,index]=exp(log_ratio_index)
    log_info(paste(
      "QR2 fold", k,
      "ridge_lambda =", f_el_drm$ridge_lambda,
      "x_dim =", f_el_drm$x_dim,
      "q_dim =", f_el_drm$q_dim,
      "theta_dim =", length(f_el_drm$theta),
      "el_drm_sec =", round(el_drm_sec, 3)
    ))

    #mu
    lasso_mu<-lasso_var_selection(S_poly_log[index_nk_A,,drop=FALSE], Y[index_nk_A], task = 0)
    selected_vars_mu<-lasso_mu$selected_vars
    mu_rf=caret::train(x=S_poly_log[index_nk_A,,drop=FALSE][,selected_vars_mu,drop=FALSE], y=Y[index_nk_A], method= "ranger", tuneGrid=generate_tune_grid(length(index_nk_A),length(selected_vars_mu),0), trControl=ctrl_mu, importance="impurity", num.trees=generate_num_tree(length(index_nk_A)) )
    pred_mu[4,index]=predict(mu_rf$finalModel, data=S_poly_log[index,,drop=FALSE][,selected_vars_mu,drop=FALSE])$predictions
    mtry_mu=c(mtry_mu,mu_rf$bestTune$mtry)
    minsize_mu=c(minsize_mu,mu_rf$bestTune$min.node.size)
    # summary(pred_mu[3,index]-mu[index])

    #tau1
    prediction_mu_1=predict(mu_rf$finalModel, data=S_poly_log[index_nk_A0,,drop=FALSE][,selected_vars_mu,drop=FALSE])
    lasso_tau1<-lasso_var_selection(X_poly[index_nk_A0,], prediction_mu_1$predictions, task = 0)
    selected_vars_tau1<-lasso_tau1$selected_vars
    tau1_rf=caret::train(x=X_poly[index_nk_A0,selected_vars_tau1,drop=FALSE], y=prediction_mu_1$predictions, method= "ranger", tuneGrid=generate_tune_grid(length(index_nk_A0),length(selected_vars_tau1),0), trControl=ctrl_mu, importance="impurity", num.trees=generate_num_tree(length(index_nk_A0)) )
    pred_tau1[4,index]=predict(tau1_rf$finalModel, data=X_poly[index,selected_vars_tau1,drop=FALSE])$predictions
    mtry_tau1=c(mtry_tau1,tau1_rf$bestTune$mtry)
    minsize_tau1=c(minsize_tau1,tau1_rf$bestTune$min.node.size)
    # summary(pred_tau1[3,index])

    #tau
    prediction_q_1 <- normalized_truncate_q(
      log_q = ref_log_q_sdr2_full,
      ref_log_q = ref_log_q_sdr2_full
    )

    prediction_mu_2_preds <- predict_caret_regression_final(
      mu_rf,
      newdata = S_poly_log[index_nk_A,,drop=FALSE][, selected_vars_mu, drop = FALSE]
    )
    prediction_tau1_1=predict(tau1_rf$finalModel, data=X_poly[index_nk_A,selected_vars_tau1,drop=FALSE])
    pred_pseudoY=prediction_q_1*(Y[index_nk_A]-prediction_mu_2_preds)+prediction_tau1_1$predictions
    lasso_tau<-lasso_var_selection(X_poly[index_nk_A,], pred_pseudoY, task = 0)
    selected_vars_tau<-lasso_tau$selected_vars
    tau_rf=caret::train(x=X_poly[index_nk_A,selected_vars_tau,drop=FALSE], y=pred_pseudoY, method= "ranger", tuneGrid=generate_tune_grid(length(index_nk_A),length(selected_vars_tau),0), trControl=ctrl_mu, importance="impurity", num.trees=generate_num_tree(length(index_nk_A)) )
    pred_tau[4,index]=predict(tau_rf$finalModel, data=X_poly[index,selected_vars_tau,drop=FALSE])$predictions
    mtry_tau=c(mtry_tau,tau_rf$bestTune$mtry)
    minsize_tau=c(minsize_tau,tau_rf$bestTune$min.node.size)
    # pi_rf=butcher(pi_rf); rm(pi_rf,prediction_pi);p_rf=butcher(p_rf); rm(p_rf);rm(prediction_pi_1,prediction_p_1);
    #q1_crulsif=butcher(q1_crulsif);q2_crulsif=butcher(q2_crulsif);rm(q1_crulsif,q2_crulsif);rm(prediction_q_1)
    mu_rf=butcher(mu_rf);tau1_rf=butcher(tau1_rf);tau_rf=butcher(tau_rf); rm(mu_rf,tau1_rf,tau_rf,prediction_mu_1,prediction_mu_2_preds,prediction_tau1_1)
    gc()
  }
  lock_obj=lock(log_file, timeout=5)
  unlock(lock_obj)
  log_info(paste("Time=", t, "QR2 accomplished"))
  set.seed(seed+t)


  #--------------O-TR'---------------
  for (k in 1:Kmb){
    index=(1:n)+(k-1)*n
    index_nk=(1:N)[-index]
    index_nk1=index_nk[1:floor(length(index_nk)/2)]
    index_nk2=index_nk[-(1:floor(length(index_nk)/2))]

    #pi
    # sel.rho_pi_a=c()
    rho.range0=sqrt(log(d1_poly)/length(index_nk))*rho.range0_default
    boolFalse=FALSE #repeat until no errors occur
    while(boolFalse==FALSE){
      tryCatch({
        fit_Log1_1=glm.regu.cv(x=X_poly[index_nk,],y=A[index_nk],fold=5,loss="ml",rho.seq=rho.range0, ann.init=FALSE,tune.cut=FALSE)
        boolFalse=TRUE
      },error=function(e){
      },finally={})
    }
    sel.rho_pi_a=c(sel.rho_pi_a,fit_Log1_1$sel.rho[sw2]/sqrt(log(d1_poly)/length(index_nk)))
    sel.rho_pi_b=c(sel.rho_pi_b,fit_Log1_1$sel.rho[sw2]/sqrt(log(d1_poly)/length(index_nk)))
    alpha_a_miles=fit_Log1_1$sel.bet[,sw2]
    pred_pi_a[5,index]=pmax(pmin(plogis(alpha_a_miles[1]+X_poly[index,]%*%alpha_a_miles[-1]),1),0) #pi_a(X)
    # if((max(pred_pi_a[5,index])>0.95) || (min(pred_pi_a[5,index])<0.05) ){
    #   alpha_a_miles=fit_Log1_1$sel.bet[,2]
    #   pred_pi_a[5,index]=pmax(pmin(plogis(alpha_a_miles[1]+X_poly[index,]%*%alpha_a_miles[-1]),1),0)
    # }
    pred_pi_b[5,index]=pred_pi_a[5,index]
    pred_pi_a[8,index]=pred_pi_a[5,index]
    pred_pi_b[8,index]=pred_pi_a[5,index]

    #p
    # sel.rho_p=c()
    rho.range0=sqrt(log(d3_poly)/length(index_nk))*rho.range1_default
    boolFalse=FALSE #repeat until no errors occur
    while(boolFalse==FALSE){
      tryCatch({
        fit_Log1_1=glm.regu.cv(x=S_poly[index_nk,], y=A[index_nk],fold=5,loss="ml",rho.seq=rho.range0, ann.init=FALSE,tune.cut=FALSE)
        boolFalse=TRUE
      },error=function(e){
      },finally={})
    }
    sel.rho_p=c(sel.rho_p,fit_Log1_1$sel.rho[sw2]/sqrt(log(d3_poly)/length(index_nk)))
    alpha_p=fit_Log1_1$sel.bet[,sw2]
    pred_p[3,index]=pmax(pmin(plogis(alpha_p[1]+S_poly[index,]%*%alpha_p[-1]),1),0) #pi_a(X)
    # if((max(pred_p[3,index])>0.95) || (min(pred_p[3,index])<0.05) ){
    #   alpha_p=fit_Log1_1$sel.bet[,2]
    #   pred_p[3,index]=pmax(pmin(plogis(alpha_p[1]+S_poly[index,]%*%alpha_p[-1]),1),0)
    # }
    pred_q[5,index]=(1-pred_p[3,index])/pred_p[3,index]*pred_pi_a[5,index]/(1-pred_pi_a[5,index])
    pred_p[5,index]=pred_p[3,index]
    pred_q[8,index]=pred_q[5,index]

    #mu-1
    # lambda_mu=c()
    fit_lasso_mu=cv.glmnet(x=S_poly[index_nk1[A[index_nk1]==1],],y=Y[index_nk1[A[index_nk1]==1]],family="gaussian",alpha=1,standardize=FALSE)
    pred_mu[5,index]=predict(fit_lasso_mu,newx=S_poly[index,],s="lambda.min",standardize=FALSE,parallel=FALSE)
    lambda_mu=c(lambda_mu,fit_lasso_mu$lambda.min/sqrt(log(d3_poly)/length(index_nk1)))

    #tau1-1
    pred_mu_new1_miles=predict(fit_lasso_mu,newx=S_poly[index_nk2[A[index_nk2]==0],],s="lambda.min",standardize=FALSE)
    if (sd(pred_mu_new1_miles)<=1e-5){pred_tau1[5,index]=rep(mean(pred_mu_new1_miles),length(index));lambda_tau1=c(lambda_tau1,0)} else{

      fit_lasso_tau1 = tryCatch({
        fit = cv.glmnet(
          x = X_poly[index_nk2[A[index_nk2]==0],],
          y = pred_mu_new1_miles,
          family = "gaussian",
          alpha = 1,
          standardize = FALSE,parallel=FALSE
        )
        list(fit = fit, error = NULL)
      }, error = function(e) {
        list(fit = NULL, error = e)
      })
      if (!is.null(fit_lasso_tau1$error)) {
        pred_tau1[5, index]=rep(mean(pred_mu_new1_miles), length(index))
        lambda_tau1 = c(lambda_tau1, -1)
      } else {
        pred_tau1[5, index]=predict(fit_lasso_tau1$fit, newx = X_poly[index, ], s = "lambda.min",standardize=FALSE)
        lambda_tau1 = c(lambda_tau1, fit_lasso_tau1$fit$lambda.min / sqrt(log(d1_poly)/length(index_nk2)))
      }

    }

   #  #mu-2
   #  # lambda_mu=c()
   # fit_lasso_mu=cv.glmnet(x=S_poly[index_nk[A[index_nk]==1],],y=Y[index_nk[A[index_nk]==1]],family="gaussian",alpha=1,standardize=FALSE,nfolds=10)
   #  pred_mu[8,index]=predict(fit_lasso_mu,newx=S_poly[index,],s="lambda.min",standardize=FALSE,parallel=FALSE)
   #  lambda_mu=c(lambda_mu,fit_lasso_mu$lambda.min/sqrt(log(d3_poly)/length(index_nk)))
   #
   #  #tau1-2
   #  pred_mu_new1_miles=predict(fit_lasso_mu,newx=S_poly[index_nk[A[index_nk]==0],],s="lambda.min",standardize=FALSE)
   #  if (sd(pred_mu_new1_miles)<=1e-5){pred_tau1[8,index]=rep(mean(pred_mu_new1_miles),length(index));lambda_tau1=c(lambda_tau1,0)} else{
   #    fit_lasso_tau1 = tryCatch({
   #      fit = cv.glmnet(
   #        x = X_poly[index_nk[A[index_nk]==0],],
   #        y = pred_mu_new1_miles,
   #        family = "gaussian",
   #        alpha = 1,
   #        standardize = FALSE,parallel=FALSE,nfolds=10
   #      )
   #      list(fit = fit, error = NULL)
   #    }, error = function(e) {
   #      list(fit = NULL, error = e)
   #    })
   #    if (!is.null(fit_lasso_tau1$error)) {
   #      pred_tau1[8, index]=#X_poly[index, ] %*% alpha_mu_X +
   #        rep(mean(pred_mu_new1_miles), length(index))
   #      lambda_tau1 = c(lambda_tau1, -1)
   #    } else {
   #      pred_tau1[8, index] = #X_poly[index, ] %*% alpha_mu_X +
   #        predict(fit_lasso_tau1$fit, newx = X_poly[index, ], s = "lambda.min",standardize=FALSE)
   #      lambda_tau1 = c(lambda_tau1, fit_lasso_tau1$fit$lambda.min / sqrt(log(d1_poly)/length(index_nk)))
   #    }
   #    #fit_lasso_tau1=cv.glmnet(x=X_poly[index_nk[A[index_nk]==0],],y=pred_mu_new1_miles,family="gaussian",alpha=1,standardize=FALSE)
   #    #pred_tau1[8,index]=predict(fit_lasso_tau1,newx=X_poly[index,],s="lambda.min",standardize=FALSE)
   #    #lambda_tau1=c(lambda_tau1,fit_lasso_tau1$lambda.min/sqrt(log(d1_poly)/length(index_nk)))
   #  }

  }
  lock_obj=lock(log_file, timeout=5)
  unlock(lock_obj)
  log_info(paste("Time=", t  ,"O-TR' accomplished"))
  set.seed(seed+t)

  #-------------TR & O-TR------------------
  for (k in 1:Kmb){
    index=(1:n)+(k-1)*n
    index_nk=(1:N)[-index]
    index_nk1=index_nk[1:floor(length(index_nk)/2)]
    index_nk2=index_nk[-(1:floor(length(index_nk)/2))]
    index_nk1_A=index_nk1[A[index_nk1]==1]
    index_nk2_A0=index_nk2[A[index_nk2]==0]
    index_nk_A=index_nk[A[index_nk]==1]
    index_nk_A0=index_nk[A[index_nk]==0]
    n_nk_A=length(index_nk_A)
    n_nk_A0=length(index_nk_A0)

    #pi using Random forest
    lasso_pi <- lasso_var_selection(X[index_nk,],A[index_nk],task=1)
    selected_vars_pi<-lasso_pi$selected_vars
    pi_rf_TR=caret::train(x=X[index_nk,selected_vars_pi,drop=FALSE], y=A_fac[index_nk], method= "ranger", trControl=ctrl_pi, tuneGrid=generate_tune_grid(length(index_nk),length(selected_vars_pi),1), metric="ROC", importance="impurity", num.trees=generate_num_tree(length(index_nk)) )
    prediction_pi_TR=predict(pi_rf_TR$finalModel, data=X[index,selected_vars_pi,drop=FALSE])
    pred_pi_a[6,index]=pmax(pmin(prediction_pi_TR$predictions[, 2],0.99),0.01)
    pred_pi_b[6,index]=pred_pi_a[6,index]
    pred_pi_a[7,index]=pred_pi_a[6,index]
    pred_pi_b[7,index]=pred_pi_a[6,index]
    pred_pi_a[9,index]=pred_pi_a[6,index]
    pred_pi_b[9,index]=pred_pi_a[6,index]
    mtry_pi=c(mtry_pi,pi_rf_TR$bestTune$mtry)
    minsize_pi=c(minsize_pi,pi_rf_TR$bestTune$min.node.size)
    pi_rf_TR=butcher(pi_rf_TR);rm(pi_rf_TR,prediction_pi_TR);gc()

    #p using Random forest-O-TR
    lasso_p <- lasso_var_selection(cbind(X[index_nk,],M[index_nk,]),A[index_nk],task=1)
    selected_vars_p <- lasso_p$selected_vars
    # selected_M=colnames(M)[select_different_M_components_ks(M[index_nk, ], A[index_nk])]
    # selected_vars_p <- c(colnames(X[index_nk,]),selected_M)
    # selected_A1=multi_y_var_selection(X[index_nk_A,], M[index_nk_A, selected_M])$union_vars
    # selected_A0=multi_y_var_selection(X[index_nk_A0,], M[index_nk_A0, selected_M])$union_vars
    # selected_vars_p <- unique(c(selected_A1,selected_A0,selected_M))
    # p_rf_TR=caret::train(x=cbind(X[index_nk,],M[index_nk,])[,selected_vars_p,drop=FALSE], y=A_fac[index_nk], method= "ranger", trControl=ctrl_pi, tuneGrid=generate_tune_grid(length(index_nk),length(selected_vars_p),1), metric="ROC", importance="impurity", num.trees=generate_num_tree(length(index_nk)) )

      p_rf_TR = caret::train(
        x = cbind(X[index_nk,], M[index_nk,])[, selected_vars_p, drop = FALSE],
        y = A_fac[index_nk],
        method = "ranger",
        trControl = ctrl_pi,
        tuneGrid = generate_tune_grid_otr_p(
          length(index_nk),
          length(selected_vars_p)
      ),
      metric = "ROC",
      importance = "impurity",
      num.trees = generate_num_tree(length(index_nk))
    )

    pred_p[4,index]=pmax(pmin(predict(p_rf_TR$finalModel, data=cbind(X[index,],M[index,])[,selected_vars_p,drop=FALSE])$predictions[, 2],0.99),0.01)
    pred_q[7,index]=(1-pred_p[4,index])/pred_p[4,index]*pred_pi_a[7,index]/(1- pred_pi_a[7,index])
    pred_p[6,index]=pred_p[4,index]
    pred_q[9,index]=pred_q[7,index]
    mtry_p=c(mtry_p,p_rf_TR$bestTune$mtry)
    minsize_p=c(minsize_p,p_rf_TR$bestTune$min.node.size)
    p_rf_TR=butcher(p_rf_TR);rm(p_rf_TR)

    #mu for O-TR
    lasso_mu<-lasso_var_selection(cbind(X[index_nk1_A,], M[index_nk1_A,]), Y[index_nk1_A], task = 0)
    selected_vars_mu<-lasso_mu$selected_vars
    mu_rf_TR=caret::train(x=cbind(X[index_nk1_A,],M[index_nk1_A,])[,selected_vars_mu,drop=FALSE], y=Y[index_nk1_A], method= "ranger", trControl=ctrl_mu, importance="impurity", num.trees=generate_num_tree(length(index_nk1_A)) , tuneGrid=generate_tune_grid(length(index_nk1_A),length(selected_vars_mu),0))
    mtry_mu=c(mtry_mu,mu_rf_TR$bestTune$mtry)
    minsize_mu=c(minsize_mu,mu_rf_TR$bestTune$min.node.size)
    pred_mu[7,index]=predict(mu_rf_TR$finalModel, data=cbind(X[index,],M[index,])[,selected_vars_mu,drop=FALSE])$predictions

    #tau1 for O-TR
    prediction_mu_new1_TR=predict(mu_rf_TR$finalModel, data=cbind(X[index_nk2_A0,],M[index_nk2_A0,])[,selected_vars_mu,drop=FALSE])
    lasso_tau1<-lasso_var_selection(X[index_nk2_A0,], prediction_mu_new1_TR$predictions, task = 0)
    selected_vars_tau1<-lasso_tau1$selected_vars
    tau1_rf_TR=caret::train(x=X[index_nk2_A0,selected_vars_tau1,drop=FALSE], y=prediction_mu_new1_TR$predictions, method= "ranger", trControl=ctrl_mu, tuneGrid=generate_tune_grid(length(index_nk2_A0),length(selected_vars_tau1),0), importance="impurity", num.trees=generate_num_tree(length(index_nk2_A0)) )
    pred_tau1[7,index]=predict(tau1_rf_TR$finalModel, data=X[index, selected_vars_tau1,drop=FALSE])$predictions
    mtry_tau1=c(mtry_tau1,tau1_rf_TR$bestTune$mtry)
    minsize_tau1=c(minsize_tau1,tau1_rf_TR$bestTune$min.node.size)

    #mu for TR/O-TR2
    lasso_mu<-lasso_var_selection(cbind(X[index_nk_A,], M[index_nk_A,]), Y[index_nk_A], task = 0)
    selected_vars_mu<-lasso_mu$selected_vars
    mu_rf_TR=caret::train(x=cbind(X[index_nk_A,],M[index_nk_A,])[,selected_vars_mu,drop=FALSE], y=Y[index_nk_A], method= "ranger", trControl=ctrl_mu, importance="impurity", num.trees=generate_num_tree(length(index_nk_A)), tuneGrid=generate_tune_grid(length(index_nk_A),length(selected_vars_mu),0))
    pred_mu[6,index]=predict(mu_rf_TR$finalModel, data=cbind(X[index,],M[index,])[,selected_vars_mu,drop=FALSE])$predictions
    pred_mu[9,index]=pred_mu[6,index]
    mtry_mu=c(mtry_mu,mu_rf_TR$bestTune$mtry)
    minsize_mu=c(minsize_mu,mu_rf_TR$bestTune$min.node.size)
    # summary(pred_mu[5,index]-mu[index])

    #tau1 for O-TR2
    prediction_mu_new1_TR=predict(mu_rf_TR$finalModel, data=cbind(X[index_nk_A0, ],M[index_nk_A0,])[,selected_vars_mu,drop=FALSE])
    lasso_tau1<-lasso_var_selection(X[index_nk_A0,], prediction_mu_new1_TR$predictions, task = 0)
    selected_vars_tau1<-lasso_tau1$selected_vars
    tau1_rf_TR=caret::train(x=as.matrix(X[index_nk_A0,selected_vars_tau1,drop=FALSE]), y=prediction_mu_new1_TR$predictions, method= "ranger", trControl=ctrl_mu, tuneGrid=generate_tune_grid(length(index_nk_A0),length(selected_vars_tau1),0), importance="impurity", num.trees=generate_num_tree(length(index_nk_A0)) )
    pred_tau1[9,index]=predict(tau1_rf_TR$finalModel, data=X[index,selected_vars_tau1,drop=FALSE])$predictions
    mtry_tau1=c(mtry_tau1,tau1_rf_TR$bestTune$mtry)
    minsize_tau1=c(minsize_tau1,tau1_rf_TR$bestTune$min.node.size)

    # Gaussian F0 for tau1 only; q uses the separate legacy KNN module.
    selected_M <- 1:ncol(M)
    ds=length(selected_M)
    x_vars <- multi_y_var_selection(X[index_nk, ], M[index_nk, selected_M, drop = FALSE])$union_vars
    selected_A1=x_vars#multi_y_var_selection(X[index_nk_A,], M[index_nk_A, selected_M])$union_vars
    selected_A0=x_vars#multi_y_var_selection(X[index_nk_A0,], M[index_nk_A0, selected_M])$union_vars

    # A=0: independently tune mtry and nodesize
    rf_mean_A0 <- fit_mv_rf_cf(
      X_train = X[index_nk_A0, selected_A0, drop = FALSE],
      M_train = M[index_nk_A0, selected_M, drop = FALSE],
      K = 3,
      ntree_tune = 50,
      ntree_final = 100,
      seed = seed
    )

    mean_rf_A0 <- rf_mean_A0$fit
    mean_train_A0 <- rf_mean_A0$oob_pred
    residuals_A0 <- M[index_nk_A0, selected_M, drop = FALSE] - mean_train_A0

    log_info(paste(
      "TR tau1 mean RF fold", k,
      "A0: mtry =", rf_mean_A0$best_mtry,
      "nodesize =", rf_mean_A0$best_nodesize
    ))
    residuals_A0 <- sweep(residuals_A0, 2, colMeans(residuals_A0), "-")
    cov_A0 <- compute_shrinkage_cov(residuals_A0, lambda = "auto")
    sig_M_A0 <- stabilize_cov_eigen(cov_A0$matrix)
    mean_pred <- predict(mean_rf_A0, newdata = data.frame(X[index, selected_A0, drop=FALSE]))
    mean_pred_A0 <- sapply(seq_len(ds), function(i) mean_pred$regrOutput[[i]][[1]])
    colnames(mean_pred_A0) <- colnames(M)[selected_M]

    # q: 175_2 KNN density, 2.5%-97.5% truncation, test-fold treated normalization.
    pred_q[6, index] <- estimate_tr_q_legacy(
      X, M, A, Y, train_index = index_nk, query_index = index,
      selected_vars_mu = selected_vars_mu
    )
    if (!all(is.finite(pred_q[6, index])) || any(pred_q[6, index] <= 0)) {
      stop("Invalid legacy TR density ratio.")
    }

    # tau_1: retain the jasa_sparse6 Gaussian F0(M|X) integration
    m_vars <- match(intersect(selected_vars_mu, colnames(M)), colnames(M))
    m_names <- colnames(M)[m_vars]
    ds <- length(m_vars)

    if (!all(m_names %in% colnames(mean_pred_A0))) {
      stop("m_vars used by mu must be contained in selected_M used for density estimation.")
    }

    mean_pred_tau1 <- mean_pred_A0[, m_names, drop=FALSE]
    sig_M_tau1 <- sig_M_A0[m_names, m_names, drop=FALSE]

    M_samples <- matrix(0, nrow=length(index)*n_mc, ncol=ds)

    for (i in seq_len(length(index))) {
      M_samples[((i-1)*n_mc+1):(i*n_mc), ] <-
        rmvnorm_eigen_custom(
          n_mc,
          mean = mean_pred_tau1[i, ],
          sigma = sig_M_tau1
        )
    }

    colnames(M_samples) <- m_names

    pred_data <- do.call(
      rbind,
      lapply(seq_len(length(index)), function(i) X[index, , drop=FALSE][rep(i, n_mc), ])
    )

    pred_data <- cbind(pred_data, M_samples)

    mu_samples <- predict(
      mu_rf_TR$finalModel,
      data=pred_data[, selected_vars_mu, drop=FALSE]
    )$predictions

    pred_tau1[6,index] <- tapply(
      mu_samples,
      rep(seq_len(length(index)), each=n_mc),
      mean
    )

  }
  lock_obj=lock(log_file, timeout=5)
  unlock(lock_obj)
  log_info(paste("Time=", t, "TR & O-TR accomplished"))
  set.seed(seed+t)

  # --------------theta00-nonpa-----------------------
  for (k in 1:Kmb){
    index=(1:n)+(k-1)*n
    index_nk=(1:N)[-index]
    index_nk_A0=index_nk[which(A[index_nk]==0)]
    index_nk_A=index_nk[which(A[index_nk]==1)]
    
    #pi-nonp
    lasso_pi <- var_selection(X[index_nk,],A[index_nk],task=1)
    selected_vars_pi <- lasso_pi$selected_vars
    pi_rf=caret::train(x=X[index_nk,selected_vars_pi,drop=FALSE], y=A_fac[index_nk], method= "ranger", trControl=ctrl_pi, metric="ROC", importance="impurity", tuneGrid = generate_tune_grid(length(index_nk),length(selected_vars_pi),1), num.trees=generate_num_tree(length(index_nk)))
    pred_pi_00[2,index]=pmax(pmin(predict(pi_rf$finalModel, data=X[index,selected_vars_pi,drop=FALSE])$predictions[,2],0.99),0.01) #pi_a(X)
    
    #mu0-nonp
    lasso_mu<-var_selection(cbind(X[index_nk_A0,], M[index_nk_A0,]), Y[index_nk_A0], task = 0)
    selected_vars_mu<-lasso_mu$selected_vars
    mu_rf=caret::train(x=cbind(X[index_nk_A0,],M[index_nk_A0,])[,selected_vars_mu,drop=FALSE], y=Y[index_nk_A0], method= "ranger", trControl=ctrl_mu, importance="impurity",  tuneGrid=generate_tune_grid(length(index_nk_A0),length(selected_vars_mu),0), num.trees=generate_num_tree(length(index_nk_A0)) )
    pred_mu_00[2,index]=predict(mu_rf$finalModel, data=cbind(X[index,],M[index,])[,selected_vars_mu,drop=FALSE])$predictions
    
    #mu1-nonp
    lasso_mu<-var_selection(cbind(X[index_nk_A,], M[index_nk_A,]), Y[index_nk_A], task = 0)
    selected_vars_mu<-lasso_mu$selected_vars
    mu_rf=caret::train(x=cbind(X[index_nk_A,],M[index_nk_A,])[,selected_vars_mu,drop=FALSE], y=Y[index_nk_A], method= "ranger", trControl=ctrl_mu, importance="impurity",  tuneGrid=generate_tune_grid(length(index_nk_A),length(selected_vars_mu),0), num.trees=generate_num_tree(length(index_nk_A)) )
    pred_mu_11[2,index]=predict(mu_rf$finalModel, data=cbind(X[index,],M[index,])[,selected_vars_mu,drop=FALSE])$predictions
  }
  lock_obj=lock(log_file, timeout=5)
  unlock(lock_obj)
  log_info(paste("Time=", t  ,"00_nonp accomplished"))
  set.seed(seed+t)
  
  #-------------------MSE check-------------
  mse_pi=rowMeans((pred_pi_a-matrix(rep(A,each=num_pi_est),ncol=N))^2)[c(1,2,5,8,3,4,6,7,9)]
  mse_mu=rowMeans((pred_mu[,A==1]-matrix(rep(Y[A==1],each=num_mu_est),nrow=num_mu_est))^2)[c(1,2,5,8,3,4,6,7,9)]
  pred_tau[5:9,]=pred_tau1[5:9,]
  
  #estimation of theta & asymptotic variance estimator
  pred_theta[1,]=A*pred_q[1,]/(pred_pi_a[1,])*(Y-pred_mu[1,])+(1-A)/(1-pred_pi_b[1,])*(pred_mu[1,]-pred_tau1[1,])+A/(pred_pi_a[1,])*(pred_tau1[1,]-pred_tau[1,])+pred_tau[1,];
  pred_11theta[1,]=A/(pred_pi_a[1,])*(Y-pred_mu[1,])+pred_mu[1,];#MQR1
  pred_theta[2,]=A*pred_q[2,]/(pred_pi_a[2,])*(Y-pred_mu[2,])+(1-A)/(1-pred_pi_b[2,])*(pred_mu[2,]-pred_tau1[2,])+A/(pred_pi_a[2,])*(pred_tau1[2,]-pred_tau[2,])+pred_tau[2,];
  pred_11theta[2,]=A/(pred_pi_a[2,])*(Y-pred_mu[2,])+pred_mu[2,];#MQR2
  pred_theta[3,]=A/(1-pred_pi_a[5,])*(1-pred_p[3,])/(pred_p[3,])*(Y-pred_mu[5,])+(1-A)/(1-pred_pi_a[5,])*(pred_mu[5,]-pred_tau1[5,])+pred_tau1[5,];
  pred_11theta[3,]=A/(pred_pi_a[5,])*(Y-pred_mu[5,])+pred_mu[5,];#O-TR' 
  pred_theta[4,]=A/(1-pred_pi_a[8,])*(1-pred_p[5,])/(pred_p[5,])*(Y-pred_mu[8,])+(1-A)/(1-pred_pi_a[8,])*(pred_mu[8,]-pred_tau1[8,])+pred_tau1[8,];
  pred_11theta[4,]=A/(pred_pi_a[8,])*(Y-pred_mu[8,])+pred_mu[8,];#MTR2 (internal comparator)
  
  pred_theta[5,]=A/(pred_pi_a[3,])*pred_q[3,]*(Y-pred_mu[3,])+(1-A)/(1-pred_pi_b[3,])*(pred_mu[3,]-pred_tau1[3,])+A/(pred_pi_a[3,])*(pred_tau1[3,]-pred_tau[3,])+pred_tau[3,];
  pred_11theta[5,]=A/(pred_pi_a[3,])*(Y-pred_mu[3,])+pred_mu[3,];#QR1
  
  pred_theta[6,]=A/(pred_pi_a[4,])*pred_q[4,]*(Y-pred_mu[4,])+(1-A)/(1-pred_pi_b[4,])*(pred_mu[4,]-pred_tau1[4,])+A/(pred_pi_a[4,])*(pred_tau1[4,]-pred_tau[4,])+pred_tau[4,];
  pred_11theta[6,]=A/(pred_pi_a[4,])*(Y-pred_mu[4,])+pred_mu[4,];#QR2
  
  pred_theta[7,]=A*pred_q[6,]/(pred_pi_a[6,])*(Y-pred_mu[6,])+(1-A)/(1-pred_pi_a[6,])*(pred_mu[6,]-pred_tau1[6,])+pred_tau1[6,];
  pred_11theta[7,]=A/(pred_pi_a[6,])*(Y-pred_mu[6,])+pred_mu[6,];#TR
  pred_theta[8,]=A/(1-pred_pi_a[7,])*(1-pred_p[4,])/(pred_p[4,])*(Y-pred_mu[7,])+(1-A)/(1-pred_pi_a[7,])*(pred_mu[7,]-pred_tau1[7,])+pred_tau1[7,];
  pred_11theta[8,]=A/(pred_pi_a[7,])*(Y-pred_mu[7,])+pred_mu[7,];#O-TR
  pred_theta[9,]=A/(1-pred_pi_a[9,])*(1-pred_p[6,])/(pred_p[6,])*(Y-pred_mu[9,])+(1-A)/(1-pred_pi_a[9,])*(pred_mu[9,]-pred_tau1[9,])+pred_tau1[9,];
  pred_11theta[9,]=A/(pred_pi_a[9,])*(Y-pred_mu[9,])+pred_mu[9,];#O-TR2
  
  pred_11theta=do.call(rbind, replicate(10,A/(pred_pi_00[2,])*(Y-pred_mu_11[2,])+pred_mu_11[2,], simplify = FALSE))
  pred_00theta=do.call(rbind, replicate(10,(1-A)/(1-pred_pi_00[2,])*(Y-pred_mu_00[2,])+pred_mu_00[2,], simplify = FALSE))
  
  theta_hat=rowMeans(pred_theta);
  theta11_hat=rowMeans(pred_11theta);
  theta00_hat=rowMeans(pred_00theta);
  sigma_hat=rowMeans((pred_theta-do.call(cbind, replicate(N, theta_hat, simplify = FALSE)))^2)
  sigma_DE=rowMeans((pred_theta-pred_00theta-do.call(cbind, replicate(N, theta_hat-theta00_hat, simplify = FALSE)))^2)
  sigma_IDE=rowMeans((pred_11theta-pred_theta-do.call(cbind, replicate(N, theta11_hat-theta_hat, simplify = FALSE)))^2)
  sigma_TE=rowMeans((pred_11theta-pred_00theta-do.call(cbind, replicate(N, theta11_hat-theta00_hat, simplify = FALSE)))^2)
  
  cat(theta_hat,"\n",sqrt(sigma_hat/N)*qnorm(1-0.05/2),"\n")
  cat(theta11_hat-theta00_hat,"\n",sqrt(sigma_TE/N)*qnorm(1-0.05/2),"\n")
  cat(theta11_hat-theta_hat,"\n",sqrt(sigma_IDE/N)*qnorm(1-0.05/2),"\n")
  cat(theta_hat-theta00_hat,"\n",sqrt(sigma_DE/N)*qnorm(1-0.05/2),"\n")
  
  result=c(pred_theta=theta_hat, pred_theta0=theta00_hat, pred_theta1=theta11_hat,pred_sig2=sigma_hat, pred_sig2_IDE=sigma_IDE, pred_sig2_DE=sigma_DE,pred_sig2_TE=sigma_TE,mse_pi=mse_pi,mse_mu=mse_mu)
  temp_file <- tempfile(pattern = "temp_", fileext = ".rds")
  saveRDS(result, file = temp_file)
  file.rename(temp_file, paste0("175_3_", current_time, "_", t , ".rds")) 
  lock_obj=lock(log_file, timeout=5)
  unlock(lock_obj)
  log_info(paste("Time=", t ,"saved"))                                   
  c(pred_theta=theta_hat,pred_theta0=theta00_hat,pred_theta1=theta11_hat,pred_sig2=sigma_hat,pred_sig2_IDE=sigma_IDE,pred_sig2_DE=sigma_DE,pred_sig2_TE=sigma_TE,mse_pi=mse_pi,mse_mu=mse_mu) 
}                    
pred_theta=results_par[,1:num_totest]
pred_theta0=results_par[,(num_totest+1):(2*num_totest)]
pred_theta1=results_par[,(2*num_totest+1):(3*num_totest)]
pred_sig2=results_par[,(3*num_totest+1):(4*num_totest)]
pred_sig2_IDE=results_par[,(4*num_totest+1):(5*num_totest)]
pred_sig2_DE=results_par[,(5*num_totest+1):(6*num_totest)]
pred_sig2_TE=results_par[,(6*num_totest+1):(7*num_totest)]

# Keep all historical fitting calls and internal legacy slots; report the seven compared methods.
method_columns <- c(1L, 2L, 3L, 5L, 6L, 7L, 8L)
method_names <- c("MQR1", "MQR2", "O-TR'", "QR1", "QR2", "TR", "O-TR")
results_par_legacy <- results_par
output_matrices <- c("pred_theta", "pred_theta0", "pred_theta1", "pred_sig2",
                     "pred_sig2_IDE", "pred_sig2_DE", "pred_sig2_TE")
for (object_name in output_matrices) {
  compact <- get(object_name)[, method_columns, drop = FALSE]
  colnames(compact) <- method_names
  assign(object_name, compact)
}
results_par <- do.call(cbind, lapply(output_matrices, function(object_name) get(object_name)))
colnames(results_par) <- unlist(lapply(output_matrices, function(metric) {
  paste(metric, method_names, sep = ":")
}))
num_totest <- length(method_names)
v=tmax

# Definition 3.3: aggregate each scalar estimand over repeated sample splits.
# pred_sig2_* are influence-function variances, not squared standard errors.
theta=apply(pred_theta[1:v,,drop=FALSE],2,median)
theta0=apply(pred_theta0[1:v,,drop=FALSE],2,median)
theta1=apply(pred_theta1[1:v,,drop=FALSE],2,median)
DE_split=pred_theta[1:v,,drop=FALSE]-pred_theta0[1:v,,drop=FALSE]
IDE_split=pred_theta1[1:v,,drop=FALSE]-pred_theta[1:v,,drop=FALSE]
TE_split=pred_theta1[1:v,,drop=FALSE]-pred_theta0[1:v,,drop=FALSE]
DE=apply(DE_split,2,median)
IDE=apply(IDE_split,2,median)
TE=apply(TE_split,2,median)

# Equation (3.14), applied separately to each scalar effect and method.
sig2=apply(pred_sig2[1:v,,drop=FALSE]+sweep(pred_theta[1:v,,drop=FALSE],2,theta,"-")^2,2,median)
sig2_IDE=apply(pred_sig2_IDE[1:v,,drop=FALSE]+sweep(IDE_split,2,IDE,"-")^2,2,median)
sig2_DE=apply(pred_sig2_DE[1:v,,drop=FALSE]+sweep(DE_split,2,DE,"-")^2,2,median)
sig2_TE=apply(pred_sig2_TE[1:v,,drop=FALSE]+sweep(TE_split,2,TE,"-")^2,2,median)

# Two-sided Wald tests of a zero effect, using the aggregated estimate and SE.
p_theta=2*pnorm(-abs(theta/sqrt(sig2/N)))
p_ide=2*pnorm(-abs(IDE/sqrt(sig2_IDE/N)))
p_de=2*pnorm(-abs(DE/sqrt(sig2_DE/N)))
p_te=2*pnorm(-abs(TE/sqrt(sig2_TE/N)))

proportion=IDE/TE
mean(proportion)

result_summary <- data.frame(
  method = method_names,
  theta = theta, theta_length = 2*qnorm(.975)*sqrt(sig2/N),
  NDE = DE, NDE_length = 2*qnorm(.975)*sqrt(sig2_DE/N), NDE_p = p_de,
  NIE = IDE, NIE_length = 2*qnorm(.975)*sqrt(sig2_IDE/N), NIE_p = p_ide,
  proportion = proportion,
  ATE = TE, ATE_length = 2*qnorm(.975)*sqrt(sig2_TE/N), ATE_p = p_te,
  row.names = NULL
)
print(result_summary)

file_name <- paste0("real_2_",  current_time, "_", seed, "_", N, "_", tmax, ".RData"); save(log_file, current_time,N_all,N,tmax,seed,method_names,method_columns,result_summary,results_par, pred_theta, pred_theta0, pred_theta1, pred_sig2, pred_sig2_DE,pred_sig2_IDE, pred_sig2_TE, proportion, sig2, sig2_IDE,sig2_DE,sig2_TE,theta, theta0,theta1,DE,IDE, TE, p_theta,p_ide,p_de,p_te, file = file_name)
# cat("\n theta=\n",paste(capture.output(theta), collapse = "\n"), file = log_file, append = TRUE)
# cat("\n proportion=\n",paste(capture.output(proportion), collapse = "\n"), file = log_file, append = TRUE)
# cat("\n sig2=\n",paste(capture.output(sig2), collapse = "\n"), "\n", file = log_file, append = TRUE)
cat("\n IDE=\n",paste(capture.output(IDE), collapse = "\n"), "\n", file = log_file, append = TRUE)
cat("\n AL/2=\n",paste(capture.output(sqrt(sig2_IDE/N)*qnorm(1-0.05/2)), collapse = "\n"), "\n", file = log_file, append = TRUE)
cat("\n DE=\n",paste(capture.output(DE), collapse = "\n"), "\n", file = log_file, append = TRUE)
cat("\n AL/2=\n",paste(capture.output(sqrt(sig2_DE/N)*qnorm(1-0.05/2)), collapse = "\n"), "\n", file = log_file, append = TRUE)
cat("\n proportion=\n",paste(capture.output(proportion), collapse = "\n"), "\n", file = log_file, append = TRUE)
cat(paste("Important data is saved in", file_name,"\n"),file = log_file, append = TRUE)
cat("\n IDE=\n",paste(capture.output(IDE), collapse = "\n"), "\n")
cat("\n AL/2=\n",paste(capture.output(sqrt(sig2_IDE/N)*qnorm(1-0.05/2)), collapse = "\n"), "\n")
cat("\n DE=\n",paste(capture.output(DE), collapse = "\n"), "\n")
cat("\n AL/2=\n",paste(capture.output(sqrt(sig2_DE/N)*qnorm(1-0.05/2)), collapse = "\n"), "\n")
cat("\n proportion=\n",paste(capture.output(proportion), collapse = "\n"), "\n")
save.image(file = paste0("real_2_workspace_",current_time,".RData"))
stopCluster(cl)
