# simu_4_N1800: reconstructed from the matched result's saved source.
# All reported methods are enabled; deprecated fits remain internal for RNG continuity.
# Only switches, replication-range plumbing, RNG preservation and output metadata are normalized.
#离散变量M
rm(list=ls())
require(glmnet);require(itertools);require(foreach);require(doParallel);require(RCAL);require(MASS);require(matrixStats);require(ranger);require(caret);require(cobalt);require(densratio);require(butcher);require(randomForestSRC);require(logger);require(filelock);require(sn);
require(truncnorm);require(mvtnorm)
ptm=proc.time()
current_time <- paste0(format(Sys.time(), "%Y%m%d_%H_%M_%S"), "_simu_4_N1800")

N=1800;d1=100L;d2=1L #sample size & dimensions
n_mc=50L;
Setting="a" #choose Setting a or b

tmax=100;Kmb=5;value_trunc=2;n=floor(N/Kmb)
# The final N=1000 replacement batch used seed 123; other sample sizes used 1234.
seed=if (N==1000) 123 else 1234;limit=30
set.seed(seed)
#13是有依赖1200100, 否则是15, 无截断无依赖10
#--------data generating--------------
# func_link=function(x){0.2+0.6*(x>0)}
func_link=function(x){(abs(x)+sin(x)+0.1)/(abs(x)+sin(x)+0.7)}
#func_link=function(x){0.2+0.6/(1+abs(x)^3)} #non-logistic function
#func_link=function(x){pnorm(x)-1/2}
beta_pi_i=1
beta_pi_X=c(rep(1,3),rep(0,d1-3))

beta_mu_i=1
beta_mu_X=c(rep(0.5,3),rep(0,d1-3))
beta_mu_M=0

beta_q_X=c(rep(1,3),rep(0,d1-3))

#X_all follows the normal distribution. 
# X_all=matrix(rtruncnorm(N*d1*tmax,-value_trunc, value_trunc,mean=0),nrow=N*tmax)
Sigma=matrix(0,nrow=d1,ncol=d1)
for (i in 1:d1){
  for (j in 1:d1){
   Sigma[i,j]=0.3^abs(i-j)
  }
}
# X_all <- rtmvnorm(n = N * tmax,  mean = rep(0, d1), sigma = Sigma, lower = rep(-2, d1), upper = rep(2, d1))
X_all=rmvnorm(n=N*tmax, mean=rep(0,d1), sigma=Sigma)

# X_all=matrix(rnorm(N*d1*tmax,mean=0, sd=1),nrow=N*tmax)
# X_all=matrix(rbeta(N*d1*tmax,1,2),nrow=N*tmax)-1/3
#X_all <- matrix(rsn(N*d1*tmax, xi = 0, omega = 1.3, alpha = 5),nrow=N*tmax)-0-1.3*5/sqrt(1 + 5^2)*sqrt(2 / pi)

#pi and A
pi_all=func_link(beta_pi_i+X_all%*%beta_pi_X)
A_all=rbinom(N*tmax,1,pi_all)

#没用上的数据
delta5_all_i=-1*1
delta5_all_X=-1*beta_mu_X
delta5_all_M=-1*beta_mu_M

#Phi
# Phi=matrix(0,nrow=d1,ncol=d2)
# for (i in 1:d1){
#   for (j in 1:d2-5){
#     if (i==j){
#       Phi[i,j]=1 
#     } else if (abs(i-j)==1){
#       Phi[i,j]=0.5
#     } else if (abs(i-j)==2){
#       Phi[i,j]=0
#     }
#   }
# }
# 
# #Sigma & U
# Sigma=diag(d2)
# Sigma[1,2]=0.3; Sigma[2,1]=0.3
# #Sigma[1,3]=0.5; Sigma[3,1]=0.5
# U_all=rmvnorm(n=N*tmax, mean=rep(0,d2), sigma=Sigma)

#M_all
# scale_XM=0.6
# scale_Me=1
# delta1_all_matrix=matrix(nrow=N*tmax,ncol=d2)
# delta2_all_matrix=matrix(nrow=N*tmax,ncol=d2)
# for (k in 1:d2){delta1_all_matrix[,k]=delta1_all[k];delta2_all_matrix[,k]=delta2_all[k]}
# M_all=scale_XM*(X_all)%*%Phi+scale_Me*U_all+delta2_all_matrix+A_all*delta1_all_matrix
# U_tot_all=scale_Me*U_all+delta2_all_matrix+A_all*delta1_all_matrix

#q
# Inverse=diag(d2)/(scale_Me)^2
# Inverse[1:3,1:3]=solve(scale_Me^2*Sigma[1:3,1:3])
# beta_q_i=as.numeric((delta2_all+1/2*delta1_all)%*%(Inverse%*%delta1_all))
# beta_q_X=scale_XM*Phi%*%(Inverse%*%delta1_all)
# beta_q_M=-Inverse%*%delta1_all
# beta_q_S=c(beta_q_i,beta_q_X,beta_q_M)
# q_all=exp(beta_q_i+X_all%*%beta_q_X+M_all%*%beta_q_M)

#M
delta1=0.3
# p_X=(0.6-delta1)/3*((X_all[,1]>0)+(X_all[,2]>0)+(X_all[,3]>0))
p_X=(0.6-delta1)*plogis(1+X_all%*%beta_q_X)

p_all=0.2+delta1*(1-A_all)*(0.5+0.5*(X_all[,2]>0))+p_X
p0_all=0.2+delta1*(0.5+0.5*(X_all[,2]>0))+p_X

p1_all=0.2+p_X
M_all=rbinom(N*tmax,1,p_all)
f1_all=M_all*p1_all+(1-M_all)*(1-p1_all)
f0_all=M_all*p0_all+(1-M_all)*(1-p0_all)
q_all=f0_all/f1_all

scale_mix=1
scale_sub=1

#mu_all
mu_all=beta_mu_i+scale_sub*(1-M_all)*(rowSums(X_all[,1:3]*X_all[,c(3,2,1)]))+scale_mix*(X_all[,3]*(X_all[,1]>=0)+X_all[,2]*(X_all[,2]>=0)+X_all[,1]*(X_all[,3]>=0))+abs(X_all-1)%*%beta_mu_X
#2*(U_tot_all[,1]>delta2_all_matrix[,1])*beta_mu_i+M_all%*%beta_mu_M+scale_mix*(M_all[,3]*abs(X_all[,1])*beta_mu_X[1]+M_all[,1]*abs(X_all[,2])*beta_mu_X[2]+M_all[,2]*abs(X_all[,3])*beta_mu_X[3])+scale_sub*((U_tot_all[,1]>delta2_all_matrix[,1]))*(abs(X_all)%*%beta_mu_X)#+(X_all[,1]>=0)*(X_all[,1]<=1/2)

# E(M | A=0, X); f0_all is the probability mass at the observed M.
bsm_all=p0_all
#tau1_all & tau_all
tau_all=beta_mu_i+scale_sub*(1-bsm_all)*(rowSums(X_all[,1:3]*X_all[,c(3,2,1)]))+scale_mix*(X_all[,3]*(X_all[,1]>=0)+X_all[,2]*(X_all[,2]>=0)+X_all[,1]*(X_all[,3]>=0))+abs(X_all-1)%*%beta_mu_X
tau1_all=tau_all

#theta
theta=mean(tau_all)

#Y
epsilon=rnorm(N*tmax,mean=0)
Y_all=mu_all+(1-A_all)*(delta5_all_i+X_all%*%delta5_all_X+M_all*delta5_all_M)+epsilon

#range of tuning parameters
if (d1>30){
  #HD
  rho.range0_default=c(seq(from=0.01,to=2,length.out=46),50,100,1000,10000)
  rho.range1_default=c(seq(from=0.01,to=2,length.out=46),50,100,1000,10000)
  rho.range2_default=c(seq(from=0.01,to=11,length.out=46),50,100,1000,10000)
  rho.range0_default_CF=c(seq(from=0.01,to=2,length.out=46),50,100,1000,10000)
  rho.range1_default_CF=c(seq(from=0.01,to=4,length.out=46),50,100,1000,10000)
  rho.range2_default_CF=c(seq(from=0.01,to=10,length.out=46),50,100,1000,10000)
} else {
  #LD
  rho.range0_default=c(seq(from=0.01,to=4,length.out=100),50,100,1000,10000)
  rho.range1_default=c(seq(from=0.01,to=5,length.out=100),50,100,1000,10000)
  rho.range2_default=c(seq(from=0.01,to=20,length.out=60),50,100,1000,10000)
  rho.range0_default_CF=c(seq(from=0.01,to=5,length.out=60),50,100,1000,10000)
  rho.range1_default_CF=c(seq(from=0.01,to=8,length.out=60),50,100,1000,10000)
  rho.range2_default_CF=c(seq(from=0.01,to=20,length.out=60),50,100,1000,10000)
}
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


standardize_data <- function(X_train, X_test = NULL) {
  mean_X <- colMeans(X_train)
  sd_X <- apply(X_train, 2, sd)
  X_train_scaled <- sweep(sweep(X_train, 2, mean_X, "-"), 2, sd_X, "/")
  if (!is.null(X_test)) {
    X_test_scaled <- sweep(sweep(X_test, 2, mean_X, "-"), 2, sd_X, "/")
    return(list(train = X_train_scaled, test = X_test_scaled))
  } else {
    return(list(train = X_train_scaled))
  }
}
# Lasso 变量筛选函数（自动标准化）
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



# Binary M: M^2=M, so only the linear basis is identifiable.
el_drm_raw_q <- function(M) {
  out <- as.matrix(M)
  colnames(out) <- "L_M1"
  out
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

# Select X predictors of the single mediator; M itself is always retained.
select_x_for_m <- function(X, M) {
  fit <- ranger(y ~ ., data = data.frame(X = X, y = as.numeric(M)),
                importance = "impurity", num.trees = 100, num.threads = 1)
  imp <- importance(fit)
  selected <- names(imp)[imp > mean(imp) + sd(imp)]
  selected <- sub("^X\\.", "", selected)
  if (length(selected) == 0L) selected <- colnames(X)
  selected
}

# Binary-M counterpart of sparse7's fit_mv_rf_cf: inner CV plus OOB predictions.
fit_binary_rf_cf <- function(X_train, M_train, K = 3,
                             mtry_grid = NULL, nodesize_grid = c(3, 5, 10, 50),
                             ntree_tune = 100, ntree_final = 50, seed = 123) {
  X_train <- as.matrix(X_train)
  y <- factor(as.numeric(M_train), levels = c(0, 1))
  n <- nrow(X_train)
  p <- ncol(X_train)
  if (is.null(mtry_grid)) {
    mtry_grid <- unique(pmax(1L, pmin(p, round(c(sqrt(p), p/10, p/5, p/3)))))
  }
  nodesize_grid <- unique(pmax(1L, pmin(n - 1L, nodesize_grid)))
  grid <- expand.grid(mtry = mtry_grid, nodesize = nodesize_grid)
  old_seed <- .Random.seed
  on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv))
  fold_id <- make_stratified_fold_id(as.numeric(M_train), K, seed)
  cv_loss <- matrix(NA_real_, nrow(grid), max(fold_id))
  for (v in seq_len(ncol(cv_loss))) {
    train <- which(fold_id != v)
    valid <- which(fold_id == v)
    for (g in seq_len(nrow(grid))) {
      fit <- ranger(x = X_train[train, , drop = FALSE], y = y[train],
                    probability = TRUE, num.trees = ntree_tune,
                    mtry = grid$mtry[g], min.node.size = grid$nodesize[g],
                    num.threads = 1, seed = seed + 1000 * v)
      prob <- predict(fit, data = X_train[valid, , drop = FALSE])$predictions[, "1"]
      cv_loss[g, v] <- mean((prob - as.numeric(M_train)[valid])^2)
    }
  }
  best <- which.min(rowMeans(cv_loss))
  fit <- ranger(x = X_train, y = y, probability = TRUE,
                num.trees = ntree_final, mtry = grid$mtry[best],
                min.node.size = grid$nodesize[best],
                num.threads = 1, seed = seed + 99999)
  list(fit = fit, oob_pred = fit$predictions[, "1"],
       best_mtry = grid$mtry[best], best_nodesize = grid$nodesize[best],
       cv_loss = rowMeans(cv_loss), tune_grid = grid)
}

#-------------experimental log----------------------
log_dir <- "ld1"
if (!dir.exists(log_dir)) {
  dir.create(log_dir)
}
log_file<-file.path(log_dir, paste0("l_ld1", current_time, ".txt"))
res_dir <- "res_ld1"
if (!dir.exists(res_dir)) {
  dir.create(res_dir)
}

log_appender(appender_file(log_file))
log_formatter(formatter_paste)
log_threshold(INFO)  # 设置日志级别为 INFO
log_info("Logger initialized and ready to record logs.")
log_message<-paste("N=", N, "; d1=", d1, "; d2=", d2, "; tmax=", tmax, ";  scale-sub=", scale_sub, "; scale_mix=", scale_mix)
write(log_message, file=log_file, append = TRUE)  # 写入日志
cat("\n beta_pi_X=", capture.output(beta_pi_X),file = log_file, append = TRUE)
cat("\n beta_mu_i=", capture.output(beta_mu_i),file = log_file, append = TRUE)
cat("\n beta_mu_X=", capture.output(beta_mu_X),file = log_file, append = TRUE)
cat("\n beta_mu_M=", capture.output(beta_mu_M),file = log_file, append = TRUE)
cat("\n summary(pi_all)=",paste(capture.output(summary(pi_all)), collapse = "\n"), file = log_file, append = TRUE)
cat("\n\n summary(p_all)=",paste(capture.output(summary(p_all)), collapse = "\n"), "\n", file = log_file, append = TRUE)
cat("\n\n summary(q_all)=",paste(capture.output(summary(q_all)), collapse = "\n"), "\n", file = log_file, append = TRUE)
cat("Simulation begun\n")

#num_cores <- detectCores()
#print(paste("Available cores:", num_cores))
if (!is.null(getOption("cl"))) { stopCluster(getOption("cl")) }
cl <- makeCluster(min(limit,tmax))  # 使用核心数减 1，保留一个核心给操作系统
registerDoParallel(cl)

#----------Computing estimators-------
num_mu_est=10;num_pi_est=10;num_q_est=10;num_p_est=6;num_totest=num_pi_est
grid=data.frame(est_mu=1:num_mu_est,est_pi=1:num_pi_est)
t=1;k=1

replication_ids <- 1:100
results_par=foreach (t=replication_ids, .combine=rbind, .packages=c("ranger", "MASS", "matrixStats","glmnet","RCAL", "mvtnorm", "caret", "densratio","randomForestSRC","butcher","logger","filelock")) %dopar% {
# for (t in 1:tmax){
  sel.rho_pi_a=c();sel.rho_pi_b=c();sel.rho_q=c(); sel.rho_p=c()
  lambda_mu=c();lambda_tau=c();lambda_tau1=c();
  
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
  
  X=X_all[(t-1)*N+1:N,];colnames(X)=paste0("X", 1:ncol(X))
  M <- as.matrix(M_all[(t-1)*N+1:N])
  colnames(M) <- "M1"
  S=cbind(X, M);colnames(S)[ncol(S)] <- "M1"
  
  A=A_all[(t-1)*N+1:N];Y=Y_all[(t-1)*N+1:N]; A_fac=factor(A, levels=c(0,1), labels=c("Control", "Treatment")); M_fac=factor(M, levels=c(0,1), labels=c("Control", "Treatment")) 
  pi=pi_all[(t-1)*N+1:N];
  # p=p_all[(t-1)*N+1:N];
  q=q_all[(t-1)*N+1:N];
  f0=f0_all[(t-1)*N+1:N];
  f1=f1_all[(t-1)*N+1:N];
  tau=tau_all[(t-1)*N+1:N];
  tau1=tau1_all[(t-1)*N+1:N];
  mu=mu_all[(t-1)*N+1:N];
  
  #polynomial transform
  #X_poly=poly(X,degree=2);M_poly=poly(M,degree=2);
  #S_poly=poly(cbind(X,M),degree=2);
  X_poly=X;S_poly=cbind(X,M);M_poly=M;  
  
  #Oracle
  pred_pi_a[num_pi_est,]=pi; pred_pi_b[num_pi_est,]=pi
  pred_q[num_q_est,]=q;
  pred_mu[num_mu_est,]=mu;
  pred_tau[num_mu_est,]=tau;
  
  log_appender(appender_file(log_file))
  log_formatter(formatter_paste)
  log_threshold(INFO)

  #estimation of pi, p, q, mu, tau1, tau
  
  # Retain the original pre-method RNG consumption and method call order.
  archive_qr1_entry_rng <- get(".Random.seed", envir = .GlobalEnv)
  # ------------MQR1-----------
  set.seed(seed)
  for (k in 1:Kmb){
    index=(1:n)+(k-1)*n
    index_nk=(1:N)[-index]
    index_nk1=index_nk[1:floor(length(index_nk)/2)]
    index_nk2=index_nk[-(1:floor(length(index_nk)/2))]
  
    index_train1=index_nk[1:floor(length(index_nk)/5)]
    index_nk_2345=index_nk[-(1:floor(length(index_nk)/5))]
    index_train23=index_nk_2345[1:floor(length(index_nk_2345)/2)]
    index_train45=index_nk_2345[-(1:floor(length(index_nk_2345)/2))]
    index_train2=index_train23[1:floor(length(index_train23)/2)]
    index_train3=index_train23[-(1:floor(length(index_train23)/2))]
    index_train4=index_train45[1:floor(length(index_train45)/2)]
    index_train5=index_train45[-(1:floor(length(index_train45)/2))]
    rm(index_nk_2345,index_train23, index_train45)
  
    # index_train12=index_nk[1:floor(length(index_nk)/2)]
    # index_train345=index_nk[-(1:floor(length(index_nk)/2))]
    # index_train1=index_train12[1:floor(length(index_train12)/2)]
    # index_train2=index_train12[-(1:floor(length(index_train12)/2))]
    # index_train3=index_train345[1:floor(length(index_train345)/2)]
    # index_train5=index_train345[-(1:floor(length(index_train345)/2))]
    # index_train4=index_train345
  
    # index_train1=index_nk[1:floor(length(index_nk)/4)]
    # index_train2=index_nk[(floor(length(index_nk)/4)+1):floor(length(index_nk)/2)]
    # index_train3=index_nk[(floor(length(index_nk)/2)+1):(floor(length(index_nk)/4*3))]
    # index_train4=index_nk[(floor(length(index_nk)/4*3)+1):(length(index_nk))]
    # index_train5=index_train4
  
  
    #pi_a
    # sel.rho_pi_a=c()
    rho.range0=sqrt(log(d1)/length(index_train1))*rho.range0_default
    boolFalse=FALSE #repeat until no errors occur
    while(boolFalse==FALSE){
      tryCatch({
        fit_Log1_1=glm.regu.cv(x=X_poly[index_train1,],y=A[index_train1],fold=5,loss="cal",rho.seq=rho.range0, ann.init=FALSE,tune.cut=FALSE)
        boolFalse=TRUE
      },error=function(e){
      },finally={})
    }
    sel.rho_pi_a=c(sel.rho_pi_a,fit_Log1_1$sel.rho[1]/sqrt(log(d1)/length(index_train1)))
    alpha_a=fit_Log1_1$sel.bet[,1]
    pred_pi_a[1,index]=plogis(alpha_a[1]+X_poly[index,]%*%alpha_a[-1]) #pi_a(X)
    # mean((abs(pred_pi_a[1,index]-pi[index])/pi[index]))
  
    #pi_b
    # sel.rho_pi_b=c()
    rho.range0=sqrt(log(d1)/length(index_train1))*rho.range0_default
    boolFalse=FALSE #repeat until no errors occur
    while(boolFalse==FALSE){
      tryCatch({
        fit_Log1_1=glm.regu.cv(x=-X_poly[index_train1,],y=1-A[index_train1],fold=5,loss="cal",rho.seq=rho.range0, ann.init=FALSE,tune.cut=FALSE)
        boolFalse=TRUE
      },error=function(e){
      },finally={})
    }
    sel.rho_pi_b=c(sel.rho_pi_b,fit_Log1_1$sel.rho[1]/sqrt(log(d1)/length(index_train2)))
    alpha_b=fit_Log1_1$sel.bet[,1]
    alpha_b[1]=-alpha_b[1]
    pred_pi_b[1,index]=plogis(alpha_b[1]+X_poly[index,]%*%alpha_b[-1]) #pi_b(X)
    # mean((abs(pred_pi_b[1,index]-pi[index])/pi[index]))
  
    #For test of imputation error of q
    # grad_pi_a=colSums(as.numeric(1-A[index_train2]/plogis(alpha_a[1]+X[index_train2,]%*%alpha_a[-1]))*X[index_train1,])/length(index_train2)
    # grad_pi_b=colSums(as.numeric(1-(1-A[index_train2])/(1-plogis(alpha_b[1]+X[index_train2,]%*%alpha_b[-1])))*X[index_train1,])/length(index_train2)
    # iw=as.numeric(A[index_train2]/plogis(alpha_a[1]+X[index_train2,]%*%alpha_a[-1])-(1-A[index_train2])/(1-plogis(alpha_b[1]+X[index_train2,]%*%alpha_b[-1])))
    # mean(iw)
    # mean(colSums(iw*X[index_train2,])/length(index_train2))
    # iw_oracle=as.numeric(A[index_train2]/pi[index_train2]-(1-A[index_train2])/(1-pi[index_train2]))
    # mean(iw_oracle)
    # mean(colSums(iw_oracle*X[index_train2,])/length(index_train2))
  
    #q
    # sel.rho_q=c()
    S_index2=S_poly[index_train2,]
    rho.range1=sqrt(log(d1+d2)/length(index_train2))*rho.range1_default
    boolFalse=FALSE #repeat until no errors occur
    while(boolFalse==FALSE){
      tryCatch({
        fit_Log2_1=glm.regu.cv(x=scale(-S_index2), y=A[index_train2], iw=as.numeric(A[index_train2]*(1+exp(-alpha_a[1]-X_poly[index_train2,]%*%alpha_a[-1]))+(1-A[index_train2])*(1+exp(alpha_b[1]+X_poly[index_train2,]%*%alpha_b[-1]))),fold=5,loss="cal",rho.seq=rho.range1, ann.init=TRUE,tune.cut=FALSE)
        boolFalse=TRUE
      },error=function(e){
      },finally={})
    }
    sel.rho_q=c(sel.rho_q,fit_Log2_1$sel.rho[1]/sqrt(log(d1+d2)/length(index_train2)))
    beta_q=fit_Log2_1$sel.bet[,1]
    mean_S_index2=colMeans(-S_index2);sd_S_index2=colSds(-S_index2);beta_q_after=rep(0,length(beta_q));beta_q_after[-1]=beta_q[-1]/sd_S_index2;beta_q_after[1]=beta_q[1]-sum((mean_S_index2/sd_S_index2)%*%beta_q[-1]);beta_q_after[1]=-beta_q_after[1]
    beta_q=beta_q_after
    pred_q[1,index]=exp(beta_q[1]+S_poly[index,]%*%beta_q[-1])  #q(S)
  
    # index_A=index[A[index]==1]
    # index_A=index_A[q[index_A]>=1]
    # cbind(mean(abs(pred_q[1,index_A]-q[index_A])),max(abs(pred_q[1,index_A]-q[index_A])),max(abs(pred_q[1,index_A]-q[index_A])/q[index_A]));
    # cbind(q[index_A[abs(pred_q[1,index_A]-q[index_A])/q[index_A]==max(abs(pred_q[1,index_A]-q[index_A])/q[index_A])]],pred_q[1,index_A[abs(pred_q[1,index_A]-q[index_A])/q[index_A]==max(abs(pred_q[1,index_A]-q[index_A])/q[index_A])]])
    # sum((beta_q-beta_q_S)^2)
  
    # max(abs(beta_q[1]+S[index, , drop = FALSE]%*%beta_q[-1]-beta_q_i-X[index,]%*%beta_q_X-M[index,]%*%beta_q_M))
    # grad_q=colSums(as.numeric(A[index_train2]/func_link(alpha_a[1]+X[index_train2,]%*%alpha_a[-1])*exp(beta_q[1]+S[index_train2, , drop = FALSE]%*%beta_q[-1])-(1-A[index_train2])/(1-func_link(alpha_b[1]+X[index_train2,]%*%alpha_b[-1])))*S[index_train1, , drop = FALSE])/length(index_train2)
    # grad_q_phi=colSums(as.numeric(A[index_train2]/func_link(alpha_a[1]+X[index_train2,]%*%alpha_a[-1]))*X[index_train2,]%*%Phi-as.numeric((1-A[index_train2])/(1-func_link(alpha_b[1]+X[index_train2,]%*%alpha_b[-1])))*M[index_train2,])/length(index_train2)
  
    #mu
    # lambda_mu=c()
    fit_lasso_mu=cv.glmnet(x=S_poly[index_train3,],y=Y[index_train3],weights=A[index_train3]*exp(beta_q[1]+S_poly[index_train3,]%*%beta_q[-1])*(1+exp(-alpha_a[1]-X_poly[index_train3,]%*%alpha_a[-1])),family="gaussian",alpha=1)
    pred_mu[1,index]=predict(fit_lasso_mu,newx=S_poly[index,],s="lambda.min")
    lambda_mu=c(lambda_mu,fit_lasso_mu$lambda.min/sqrt(log(d1+d2)/length(index_train3)))
    # beta_mu=as.vector(coef(fit_lasso_mu$glmnet.fit, s=fit_lasso_mu$lambda.min))
    # summary(pred_mu[1,index]-mu[index])
  
    #tau1
    index_train4_A0=index_train4[which(A[index_train4]==0)]
    find_mu=which(fit_lasso_mu$lambda==fit_lasso_mu$lambda.min)
    a0_mu=fit_lasso_mu$glmnet.fit$a0[find_mu]
    alpha_mu_X=fit_lasso_mu$glmnet.fit$beta[1:d1,find_mu]
    alpha_mu_M=fit_lasso_mu$glmnet.fit$beta[(d1+1):(d1+d2),find_mu]
    pred_mu_new1=S_poly[index_train4_A0,(d1+1):(d1+d2),drop=FALSE]%*%alpha_mu_M+a0_mu
    #pred_mu_new1=predict(fit_lasso_mu,newx=S_poly[index_train4,],s="lambda.min")
    if (sd(pred_mu_new1)==0 || length(pred_mu_new1)<8){pred_tau1[1,index]=X_poly[index,]%*%alpha_mu_X+
      rep(mean(pred_mu_new1),length(index)); lambda_tau1=c(lambda_tau1,-1)} else{
        fit_lasso_tau1=cv.glmnet(x=X_poly[index_train4_A0,],y=pred_mu_new1, weights=exp(alpha_b[1]+X_poly[index_train4_A0,]%*%alpha_b[-1]),family="gaussian",alpha=1)
        pred_tau1[1,index]=X_poly[index,]%*%alpha_mu_X+
          predict(fit_lasso_tau1,newx=X_poly[index,],s="lambda.min");
        lambda_tau1=c(lambda_tau1,fit_lasso_tau1$lambda.min/sqrt(log(d1)/length(index_train4)))}
    # beta_tau1=as.vector(coefficients=coef(fit_lasso_tau1$glmnet.fit, s=fit_lasso_tau1$lambda.min))
    # pred_tau1_test=pred_tau1[1,index]
  
  
    #tau
    pred_mu_new2=S_poly[index_train5,(d1+1):(d1+d2),drop=FALSE]%*%alpha_mu_M+a0_mu
    #pred_mu_new2=predict(fit_lasso_mu,newx=S_poly[index_train5,], s="lambda.min")
    if (sd(pred_mu_new1)==0 || length(pred_mu_new1)<8){pred_tau1_new1=X_poly[index_train5,]%*%alpha_mu_X+
      rep(mean(pred_mu_new1),length(index_train5))} else{      pred_tau1_new1=X_poly[index_train5,]%*%alpha_mu_X+
        predict(fit_lasso_tau1,newx=X_poly[index_train5,],s="lambda.min")
      }
    y_pseudo=pred_tau1_new1+exp(beta_q[1]+S_poly[index_train5,]%*%beta_q[-1])*(Y[index_train5]-pred_mu_new2-X_poly[index_train5,]%*%alpha_mu_X)
    if (sd(y_pseudo)==0){pred_tau[1,index]=rep(mean(y_pseudo),length(index));lambda_tau=c(lambda_tau,0)} else{
      fit_lasso_tau=cv.glmnet(x=X_poly[index_train5,],y=y_pseudo, weights=A[index_train5]*exp(-alpha_a[1]-X_poly[index_train5,]%*%alpha_a[-1]),family="gaussian",alpha=1)
      pred_tau[1,index]=predict(fit_lasso_tau,newx=X_poly[index,],s="lambda.min")
      lambda_tau=c(lambda_tau,fit_lasso_tau$lambda.min/sqrt(log(d1)/length(index_train4)))}
  }
  lock_obj=lock(log_file, timeout=5)
  unlock(lock_obj)
  log_info(paste("Time=", t, "MQR1 accomplished"))
  set.seed(seed)

  #------------MQR2------------------
  set.seed(seed)
  for (k in 1:Kmb){
    index=(1:n)+(k-1)*n
    index_nk=(1:N)[-index]
    index_nk_A=index_nk[A[index_nk]==1]
    index_nk_A0=index_nk[A[index_nk]==0]
  
    #pi_a
    rho.range0=sqrt(log(d1)/length(index_nk))*rho.range0_default
    boolFalse=FALSE #repeat until no errors occur
    while(boolFalse==FALSE){
      tryCatch({
        fit_Log1_1=glm.regu.cv(x=X_poly[index_nk,],y=A[index_nk],fold=5,loss="cal",rho.seq=rho.range0, ann.init=FALSE,tune.cut=FALSE)
        boolFalse=TRUE
      },error=function(e){
      },finally={})
    }
    sel.rho_pi_a=c(sel.rho_pi_a,fit_Log1_1$sel.rho[1]/sqrt(log(d1)/length(index_nk)))
    alpha_a=fit_Log1_1$sel.bet[,1]
    pred_pi_a[2,index]=plogis(alpha_a[1]+X_poly[index,]%*%alpha_a[-1]) #pi_a(X)
    # mean((abs(pred_pi_a[2,index]-pi[index])/pi[index]))
  
    #pi_b
    rho.range0=sqrt(log(d1)/length(index_nk))*rho.range0_default
    boolFalse=FALSE #repeat until no errors occur
    while(boolFalse==FALSE){
      tryCatch({
        fit_Log1_1=glm.regu.cv(x=-X_poly[index_nk,],y=1-A[index_nk],fold=5,loss="cal",rho.seq=rho.range0, ann.init=FALSE,tune.cut=FALSE)
        boolFalse=TRUE
      },error=function(e){
      },finally={})
    }
    sel.rho_pi_b=c(sel.rho_pi_b,fit_Log1_1$sel.rho[1]/sqrt(log(d1)/length(index_nk)))
    alpha_b=fit_Log1_1$sel.bet[,1]
    alpha_b[1]=-alpha_b[1]
    pred_pi_b[2,index]=plogis(alpha_b[1]+X_poly[index,]%*%alpha_b[-1]) #pi_b(X)
    # mean((abs(pred_pi_b[2,index]-pi[index])/pi[index]))
  
    #q
    S_index2=S_poly[index_nk,]
    rho.range1=sqrt(log(d1+d2)/length(index_nk))*rho.range1_default
    boolFalse=FALSE #repeat until no errors occur
    while(boolFalse==FALSE){
      tryCatch({
        fit_Log2_1=glm.regu.cv(x=scale(-S_index2), y=A[index_nk], iw=as.numeric(A[index_nk]/plogis(alpha_a[1]+X_poly[index_nk,]%*%alpha_a[-1])+(1-A[index_nk])/(1-plogis(alpha_b[1]+X_poly[index_nk,]%*%alpha_b[-1]))),fold=5,loss="cal",rho.seq=rho.range1, ann.init=TRUE,tune.cut=FALSE)
        # fit_Log2_1=glm.regu.cv(x=-S_poly[index_nk,], y=A[index_nk], iw=as.numeric(A[index_nk]/plogis(alpha_a[1]+X_poly[index_nk,]%*%alpha_a[-1])+(1-A[index_nk])/(1-plogis(alpha_b[1]+X_poly[index_nk,]%*%alpha_b[-1]))),fold=5,loss="cal",rho.seq=rho.range1, ann.init=FALSE,tune.cut=FALSE)
        boolFalse=TRUE
      },error=function(e){
      },finally={})
    }
    sel.rho_q=c(sel.rho_q,fit_Log2_1$sel.rho[1]/sqrt(log(d1+d2)/length(index_nk)))
    beta_q=fit_Log2_1$sel.bet[,1]
    mean_S_index2=colMeans(-S_index2);sd_S_index2=colSds(-S_index2);beta_q_after=rep(0,length(beta_q));beta_q_after[-1]=beta_q[-1]/sd_S_index2;beta_q_after[1]=beta_q[1]-sum((mean_S_index2/sd_S_index2)%*%beta_q[-1]);beta_q_after[1]=-beta_q_after[1]
    beta_q[1]=-beta_q[1]
    beta_q=beta_q_after
    pred_q[2,index]=exp(beta_q[1]+S_poly[index,]%*%beta_q[-1])  #q(S)
  
    # index_A=index[A[index]==1]
    # index_A=index_A[q[index_A]>=1]
    # cbind(mean(abs(pred_q[2,index_A]-q[index_A])),max(abs(pred_q[2,index_A]-q[index_A])),max(abs(pred_q[2,index_A]-q[index_A])/q[index_A]));
    # cbind(q[index_A[abs(pred_q[2,index_A]-q[index_A])/q[index_A]==max(abs(pred_q[2,index_A]-q[index_A])/q[index_A])]],pred_q[2,index_A[abs(pred_q[2,index_A]-q[index_A])/q[index_A]==max(abs(pred_q[2,index_A]-q[index_A])/q[index_A])]])
    # sum((beta_q-beta_q_S)^2)
  
    #mu
    fit_lasso_mu=cv.glmnet(x=S_poly[index_nk_A,],y=Y[index_nk_A],weights=exp(beta_q[1]+S_poly[index_nk_A,]%*%beta_q[-1])*(1+exp(-alpha_a[1]-X_poly[index_nk_A,]%*%alpha_a[-1])),family="gaussian",alpha=1)
    pred_mu[2,index]=predict(fit_lasso_mu,newx=S_poly[index,],s="lambda.min")
    lambda_mu=c(lambda_mu,fit_lasso_mu$lambda.min/sqrt(log(d1+d2)/length(index_nk)))
    # beta_mu=as.vector(coefficients=coef(fit_lasso_mu$glmnet.fit, s=fit_lasso_mu$lambda.min))
    # summary(pred_mu[2,index]-mu[index])
    # summary(A[index_A]*q[index_A]/(pi[index_A])*(Y[index_A]-mu[index_A]))
    # summary(A[index_A]*q[index_A]/(pi[index_A])*(Y[index_A]-mu[index_A])-A[index_A]*pred_q[2,index_A]/(pred_pi_a[2,index_A])*(Y[index_A]-pred_mu[2,index_A]))
    # summary(A[index_A]*pred_q[2,index_A]/(pred_pi_a[2,index_A])*(Y[index_A]-pred_mu[2,index_A]))
  
    #tau1
    find_mu=which(fit_lasso_mu$lambda==fit_lasso_mu$lambda.min)
    a0_mu=fit_lasso_mu$glmnet.fit$a0[find_mu]
    alpha_mu_X=fit_lasso_mu$glmnet.fit$beta[1:d1,find_mu]
    alpha_mu_M=fit_lasso_mu$glmnet.fit$beta[(d1+1):(d1+d2),find_mu]
    pred_mu_new1=a0_mu+S_poly[index_nk_A0,(d1+1):(d1+d2),drop=FALSE]%*%alpha_mu_M
    #pred_mu_new1=predict(fit_lasso_mu,newx=S_poly[index_nk,],s="lambda.min")
    if (sd(pred_mu_new1)==0){pred_tau1[2,index]=X_poly[index,]%*%alpha_mu_X+
      rep(mean(pred_mu_new1),length(index));lambda_tau1=c(lambda_tau1,-1)} else{
        fit_lasso_tau1=cv.glmnet(x=X_poly[index_nk_A0,],y=pred_mu_new1, weights=exp(alpha_b[1]+X_poly[index_nk_A0,]%*%alpha_b[-1]),family="gaussian",alpha=1)
        pred_tau1[2,index]=X_poly[index,]%*%alpha_mu_X+
          predict(fit_lasso_tau1,newx=X_poly[index,],s="lambda.min")
        lambda_tau1=c(lambda_tau1,fit_lasso_tau1$lambda.min/sqrt(log(d1)/length(index_nk)))
      }
    # beta_tau1=as.vector(coefficients=coef(fit_lasso_tau1$glmnet.fit, s=fit_lasso_tau1$lambda.min))
    # pred_tau1_test=pred_tau1[2,index]
  
    #tau
    pred_mu_new1=a0_mu+S_poly[index_nk_A,(d1+1):(d1+d2),drop=FALSE]%*%alpha_mu_M
    if (sd(pred_mu_new1)==0){pred_tau1_new1=X_poly[index_nk_A,]%*%alpha_mu_X+
      rep(mean(pred_mu_new1),length(index_nk_A))} else{
        pred_tau1_new1=X_poly[index_nk_A,]%*%alpha_mu_X+
          predict(fit_lasso_tau1,newx=X_poly[index_nk_A,],s="lambda.min")
      }
    y_pseudo=pred_tau1_new1+exp(beta_q[1]+S_poly[index_nk_A,]%*%beta_q[-1])*(Y[index_nk_A]-pred_mu_new1-X_poly[index_nk_A,]%*%alpha_mu_X)
    if (sd(y_pseudo)==0){pred_tau[2,index]=rep(mean(y_pseudo),length(index));lambda_tau=c(lambda_tau,-1)} else{
      fit_lasso_tau=cv.glmnet(x=X_poly[index_nk_A,],y=y_pseudo, weights=exp(-alpha_a[1]-X_poly[index_nk_A,]%*%alpha_a[-1]),family="gaussian",alpha=1)
      pred_tau[2,index]=predict(fit_lasso_tau,newx=X_poly[index,],s="lambda.min")
      lambda_tau=c(lambda_tau,fit_lasso_tau$lambda.min/sqrt(log(d1)/length(index_nk)))
    }
  }
  lock_obj=lock(log_file, timeout=5)
  unlock(lock_obj)
  log_info(paste("Time=", t, "MQR2 accomplished"))
  set.seed(seed)

  # Preserve the QR1 entry state of its published QR-only batch.
  assign(".Random.seed", archive_qr1_entry_rng, envir = .GlobalEnv)
  #----------------QR1-----------
  set.seed(seed)
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
    lasso_pi <- lasso_var_selection(X[index_nk12,],A[index_nk12],task=1)
    selected_vars_pi <- lasso_pi$selected_vars
    pi_rf=caret::train(x=X[index_nk12,selected_vars_pi,drop=FALSE], y=A_fac[index_nk12], method= "ranger", trControl=ctrl_pi, tuneGrid=generate_tune_grid(length(index_nk12),length(selected_vars_pi),1), metric="ROC", importance="none", num.trees=generate_num_tree(length(index_nk12)))#,class.weights = c("Control" = 0.5, "Treatment" = 0.5)
    prediction_pi=predict(pi_rf$finalModel, data=X[index,selected_vars_pi,drop=FALSE])
    pred_pi_a[3,index]=pmax(pmin(prediction_pi$predictions[,2], 0.99), 0.01)
    pred_pi_b[3,index]=pred_pi_a[3,index]
    mtry_pi=c(mtry_pi,pi_rf$bestTune$mtry)
    minsize_pi=c(minsize_pi,pi_rf$bestTune$min.node.size)
  
    # f semiparametrically
    # M is binary and one-dimensional: no mediator or quadratic-basis selection.
    x_vars <- select_x_for_m(X[index_nk12, , drop = FALSE], M[index_nk12, , drop = FALSE])
    selected_basis <- "L_M1"
    log_info(paste("QR1 basis fold", k, "; selected_basis =", selected_basis))

    tm_el_drm <- proc.time()[["elapsed"]]
    f_el_drm_sdr1 <- el_drm_fit_density_ratio(
      X[index_nk12, x_vars, drop = FALSE],
      M[index_nk12, , drop = FALSE],
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
      X[index, x_vars, drop = FALSE],
      M[index, , drop = FALSE]
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
    lasso_mu<-lasso_var_selection(S[index_nk1_A, , drop = FALSE], Y[index_nk1_A], task = 0)
    selected_vars_mu<-lasso_mu$selected_vars
    mu_rf=caret::train(x=S[index_nk1_A, , drop = FALSE][,selected_vars_mu,drop=FALSE], y=Y[index_nk1_A], method= "ranger", tuneGrid=generate_tune_grid(length(index_nk1_A),length(selected_vars_mu),0), trControl=ctrl_mu, importance="none", num.trees=generate_num_tree(length(index_nk1_A)) )
    pred_mu[3,index]=predict(mu_rf$finalModel, data=S[index, , drop = FALSE][,selected_vars_mu,drop=FALSE])$predictions
    mtry_mu=c(mtry_mu,mu_rf$bestTune$mtry)
    minsize_mu=c(minsize_mu,mu_rf$bestTune$min.node.size)
    # summary(pred_mu[3,index]-mu[index])
  
    #tau1
    prediction_mu_1=predict(mu_rf$finalModel, data=S[index_nk2_A0, , drop = FALSE][,selected_vars_mu,drop=FALSE])
    lasso_tau1<-lasso_var_selection(X[index_nk2_A0,], prediction_mu_1$predictions, task = 0)
    selected_vars_tau1<-lasso_tau1$selected_vars
    tau1_rf=caret::train(x=X[index_nk2_A0,selected_vars_tau1,drop=FALSE], y=prediction_mu_1$predictions, method= "ranger", tuneGrid=generate_tune_grid(length(index_nk2_A0),length(selected_vars_tau1),0), trControl=ctrl_mu, importance="none", num.trees=generate_num_tree(length(index_nk2_A0)))
    mtry_tau1=c(mtry_tau1,tau1_rf$bestTune$mtry)
    minsize_tau1=c(minsize_tau1,tau1_rf$bestTune$min.node.size)
    pred_tau1[3,index]=predict(tau1_rf$finalModel, data=X[index,selected_vars_tau1,drop=FALSE])$predictions
    # summary(pred_tau1[3,index])
  
    #tau
    # prediction_pi_1=pmax(pmin(predict(pi_rf$finalModel, data=X[index_nk3_A,selected_vars_pi,drop=FALSE])$predictions[, 2], 0.99), 0.01)
    # prediction_p_1=pmax(pmin(predict(p_rf$finalModel, data=S[index_nk3_A, , drop = FALSE][,selected_vars_p,drop=FALSE])$predictions[, 2], 0.99), 0.01)
    # prediction_q_1=(1-prediction_p_1)/prediction_p_1*prediction_pi_1/(1-prediction_pi_1)
  
    log_ratio_index_nk3_A <- el_drm_predict_log_ratio(f_el_drm_sdr1,
                                                      X[index_nk3_A, x_vars, drop = FALSE],
                                                      M[index_nk3_A, , drop = FALSE])
    prediction_q_1 <- normalized_truncate_q(
      log_q = log_ratio_index_nk3_A,
      ref_log_q = ref_log_q_sdr1
    )
    # prediction_q_1=exp(pmax(pmin(log_ratio_index_nk3_A, log(100)), log(0.01)))
    # prediction_q_1=exp(log_ratio_index_nk3_A)
  
    prediction_mu_2=predict(mu_rf$finalModel, data=S[index_nk3_A, , drop = FALSE][,selected_vars_mu,drop=FALSE])
    prediction_tau1_1=predict(tau1_rf$finalModel, data=X[index_nk3_A,selected_vars_tau1,drop=FALSE])
    pred_pseudoY=prediction_q_1*(Y[index_nk3_A]-prediction_mu_2$predictions)+prediction_tau1_1$predictions
    lasso_tau<-lasso_var_selection(X[index_nk3_A,], pred_pseudoY, task = 0)
    selected_vars_tau<-lasso_tau$selected_vars
    tau_rf=caret::train(x=X[index_nk3_A,selected_vars_tau,drop=FALSE], y=pred_pseudoY, method= "ranger", tuneGrid=generate_tune_grid(length(index_nk3_A),length(selected_vars_tau),0), trControl=ctrl_mu, importance="none", num.trees=generate_num_tree(length(index_nk3_A)) )
    mtry_tau=c(mtry_tau,tau_rf$bestTune$mtry)
    minsize_tau=c(minsize_tau,tau_rf$bestTune$min.node.size)
  
    pred_tau[3,index]=predict(tau_rf$finalModel, data=X[index,selected_vars_tau,drop=FALSE])$predictions
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
  set.seed(seed)
  
  # --------------------QR2---------
  set.seed(seed)
  for (k in 1:Kmb) {
    index=(1:n)+(k-1)*n
    index_nk=(1:N)[-index]
    index_nk_A=index_nk[A[index_nk]==1]
    index_nk_A0=index_nk[A[index_nk]==0]
  
    #pi using Random forest
    lasso_pi <- lasso_var_selection(X[index_nk,],A[index_nk],task=1)
    selected_vars_pi <- lasso_pi$selected_vars
    pi_rf=caret::train(x=X[index_nk,selected_vars_pi,drop=FALSE], y=A_fac[index_nk], method= "ranger", trControl=ctrl_pi, tuneGrid=generate_tune_grid(length(index_nk),length(selected_vars_pi),1), metric="ROC", importance="none", num.trees=generate_num_tree(length(index_nk)))
    prediction_pi=predict(pi_rf$finalModel, data=X[index,selected_vars_pi,drop=FALSE])
    pred_pi_a[4,index]=pmax(pmin(prediction_pi$predictions[, 2],0.99), 0.01)
    pred_pi_b[4,index]=pred_pi_a[4,index]
    mtry_pi=c(mtry_pi,pi_rf$bestTune$mtry)
    minsize_pi=c(minsize_pi,pi_rf$bestTune$min.node.size)
  
    # q using EL-DRM
    # M is binary and one-dimensional: no mediator or quadratic-basis selection.
    x_vars <- select_x_for_m(X[index_nk, , drop = FALSE], M[index_nk, , drop = FALSE])
    selected_basis <- "L_M1"
    log_info(paste("QR2 basis fold", k, "; selected_basis =", selected_basis))

    tm_el_drm <- proc.time()[["elapsed"]]
    f_el_drm <- el_drm_fit_density_ratio(
      X[index_nk, x_vars, drop = FALSE],
      M[index_nk, , drop = FALSE],
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
                                                X[index, x_vars, drop = FALSE],
                                                M[index, , drop = FALSE])
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
    lasso_mu<-lasso_var_selection(S[index_nk_A, , drop = FALSE], Y[index_nk_A], task = 0)
    selected_vars_mu<-lasso_mu$selected_vars
    mu_rf=caret::train(x=S[index_nk_A, , drop = FALSE][,selected_vars_mu,drop=FALSE], y=Y[index_nk_A], method= "ranger", tuneGrid=generate_tune_grid(length(index_nk_A),length(selected_vars_mu),0), trControl=ctrl_mu, importance="none", num.trees=generate_num_tree(length(index_nk_A)) )
    pred_mu[4,index]=predict(mu_rf$finalModel, data=S[index, , drop = FALSE][,selected_vars_mu,drop=FALSE])$predictions
    mtry_mu=c(mtry_mu,mu_rf$bestTune$mtry)
    minsize_mu=c(minsize_mu,mu_rf$bestTune$min.node.size)
    # summary(pred_mu[3,index]-mu[index])
  
    #tau1
    prediction_mu_1=predict(mu_rf$finalModel, data=S[index_nk_A0, , drop = FALSE][,selected_vars_mu,drop=FALSE])
    lasso_tau1<-lasso_var_selection(X[index_nk_A0,], prediction_mu_1$predictions, task = 0)
    selected_vars_tau1<-lasso_tau1$selected_vars
    tau1_rf=caret::train(x=X[index_nk_A0,selected_vars_tau1,drop=FALSE], y=prediction_mu_1$predictions, method= "ranger", tuneGrid=generate_tune_grid(length(index_nk_A0),length(selected_vars_tau1),0), trControl=ctrl_mu, importance="none", num.trees=generate_num_tree(length(index_nk_A0)) )
    pred_tau1[4,index]=predict(tau1_rf$finalModel, data=X[index,selected_vars_tau1,drop=FALSE])$predictions
    mtry_tau1=c(mtry_tau1,tau1_rf$bestTune$mtry)
    minsize_tau1=c(minsize_tau1,tau1_rf$bestTune$min.node.size)
    # summary(pred_tau1[3,index])
  
    #tau: use full outer-training-fold predictions for the QR2 pseudo-outcome.
    prediction_q_1 <- normalized_truncate_q(
      log_q = ref_log_q_sdr2_full,
      ref_log_q = ref_log_q_sdr2_full
    )
  
    prediction_mu_2_preds <- predict_caret_regression_final(
      mu_rf,
      newdata = S[index_nk_A, , drop = FALSE][, selected_vars_mu, drop = FALSE]
    )
    prediction_tau1_1=predict(tau1_rf$finalModel, data=X[index_nk_A,selected_vars_tau1,drop=FALSE])
    pred_pseudoY=prediction_q_1*(Y[index_nk_A]-prediction_mu_2_preds)+prediction_tau1_1$predictions
    lasso_tau<-lasso_var_selection(X[index_nk_A,], pred_pseudoY, task = 0)
    selected_vars_tau<-lasso_tau$selected_vars
    tau_rf=caret::train(x=X[index_nk_A,selected_vars_tau,drop=FALSE], y=pred_pseudoY, method= "ranger", tuneGrid=generate_tune_grid(length(index_nk_A),length(selected_vars_tau),0), trControl=ctrl_mu, importance="none", num.trees=generate_num_tree(length(index_nk_A)) )
    pred_tau[4,index]=predict(tau_rf$finalModel, data=X[index,selected_vars_tau,drop=FALSE])$predictions
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
  set.seed(seed)
  
  #--------------O-TR'---------------
  set.seed(seed)
  for (k in 1:Kmb){
    index=(1:n)+(k-1)*n
    index_nk=(1:N)[-index]
    index_nk1=index_nk[1:floor(length(index_nk)/2)]
    index_nk2=index_nk[-(1:floor(length(index_nk)/2))]
  
    #pi
    # sel.rho_pi_a=c()
    rho.range0=sqrt(log(d1)/length(index_nk))*rho.range0_default
    boolFalse=FALSE #repeat until no errors occur
    while(boolFalse==FALSE){
      tryCatch({
        fit_Log1_1=glm.regu.cv(x=X_poly[index_nk,],y=A[index_nk],fold=5,loss="cal",rho.seq=rho.range0, ann.init=FALSE,tune.cut=FALSE)
        boolFalse=TRUE
      },error=function(e){
      },finally={})
    }
    sel.rho_pi_a=c(sel.rho_pi_a,fit_Log1_1$sel.rho[1]/sqrt(log(d1)/length(index_nk)))
    sel.rho_pi_b=c(sel.rho_pi_b,fit_Log1_1$sel.rho[1]/sqrt(log(d1)/length(index_nk)))
    alpha_a_miles=fit_Log1_1$sel.bet[,1]
    pred_pi_a[5,index]=plogis(alpha_a_miles[1]+X_poly[index,]%*%alpha_a_miles[-1]) #pi_a(X)
    pred_pi_b[5,index]=pred_pi_a[5,index]
    pred_pi_a[8,index]=pred_pi_a[5,index]
    pred_pi_b[8,index]=pred_pi_a[5,index]
  
    #p
    # sel.rho_p=c()
    rho.range0=sqrt(log(d1+d2)/length(index_nk))*rho.range0_default
    boolFalse=FALSE #repeat until no errors occur
    while(boolFalse==FALSE){
      tryCatch({
        fit_Log1_1=glm.regu.cv(x=S_poly[index_nk,], y=A[index_nk],fold=5,loss="cal",rho.seq=rho.range0, ann.init=FALSE,tune.cut=FALSE)
        boolFalse=TRUE
      },error=function(e){
      },finally={})
    }
    sel.rho_p=c(sel.rho_p,fit_Log1_1$sel.rho[1]/sqrt(log(d1+d2)/length(index_nk)))
    alpha_p=fit_Log1_1$sel.bet[,1]
    pred_p[3,index]=plogis(alpha_p[1]+S_poly[index,]%*%alpha_p[-1]) #pi_a(X)
    pred_q[5,index]=(1-pred_p[3,index])/pred_p[3,index]*pred_pi_a[5,index]/(1-pred_pi_a[5,index])
    pred_p[5,index]=pred_p[3,index]
    pred_q[8,index]=pred_q[5,index]
  
    #mu-1
    # lambda_mu=c()
    fit_lasso_mu=cv.glmnet(x=S_poly[index_nk1[A[index_nk1]==1],],y=Y[index_nk1[A[index_nk1]==1]],family="gaussian",alpha=1)
    pred_mu[5,index]=predict(fit_lasso_mu,newx=S_poly[index,],s="lambda.min")
    lambda_mu=c(lambda_mu,fit_lasso_mu$lambda.min/sqrt(log(d1+d2)/length(index_nk1)))
    # beta_mu=as.vector(coefficients=coef(fit_lasso_mu$glmnet.fit, s=fit_lasso_mu$lambda.min))
    # summary(pred_mu[4,index]-mu[index])
  
    #tau1-1
    pred_mu_new1_miles=predict(fit_lasso_mu,newx=S_poly[index_nk2[A[index_nk2]==0],],s="lambda.min")
    if (sd(pred_mu_new1_miles)==0){pred_tau1[5,index]=rep(mean(pred_mu_new1_miles),length(index));lambda_tau1=c(lambda_tau1,0)} else{
      fit_lasso_tau1=cv.glmnet(x=X_poly[index_nk2[A[index_nk2]==0],],y=pred_mu_new1_miles,family="gaussian",alpha=1)
      pred_tau1[5,index]=predict(fit_lasso_tau1,newx=X_poly[index,],s="lambda.min")
      lambda_tau1=c(lambda_tau1,fit_lasso_tau1$lambda.min/sqrt(log(d1)/length(index_nk2)))
    }
  
    #mu-2
    # lambda_mu=c()
    fit_lasso_mu=cv.glmnet(x=S_poly[index_nk[A[index_nk]==1],],y=Y[index_nk[A[index_nk]==1]],family="gaussian",alpha=1)
    pred_mu[8,index]=predict(fit_lasso_mu,newx=S_poly[index,],s="lambda.min")
    lambda_mu=c(lambda_mu,fit_lasso_mu$lambda.min/sqrt(log(d1+d2)/length(index_nk)))
  
    #tau1-2
    pred_mu_new1_miles=predict(fit_lasso_mu,newx=S_poly[index_nk[A[index_nk]==0],],s="lambda.min")
    if (sd(pred_mu_new1_miles)==0){pred_tau1[8,index]=rep(mean(pred_mu_new1_miles),length(index));lambda_tau1=c(lambda_tau1,0)} else{
      fit_lasso_tau1=cv.glmnet(x=X_poly[index_nk[A[index_nk]==0],],y=pred_mu_new1_miles,family="gaussian",alpha=1)
      pred_tau1[8,index]=predict(fit_lasso_tau1,newx=X_poly[index,],s="lambda.min")
      lambda_tau1=c(lambda_tau1,fit_lasso_tau1$lambda.min/sqrt(log(d1)/length(index_nk)))
    }
  
  }
  lock_obj=lock(log_file, timeout=5)
  unlock(lock_obj)
  log_info(paste("Time=", t, "O-TR' accomplished"))
  set.seed(seed)
  

  #-------------TR & O-TR------------------
  set.seed(seed)
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
    pi_rf_TR=caret::train(x=X[index_nk,selected_vars_pi,drop=FALSE], y=A_fac[index_nk], method= "ranger", trControl=ctrl_pi, tuneGrid=generate_tune_grid(length(index_nk),length(selected_vars_pi),1), metric="ROC", importance="none", num.trees=generate_num_tree(length(index_nk)) )
    prediction_pi_TR=predict(pi_rf_TR$finalModel, data=X[index,selected_vars_pi,drop=FALSE])
    pred_pi_a[6,index]=pmax(pmin(prediction_pi_TR$predictions[, 2],0.99), 0.01)
    pred_pi_b[6,index]=pred_pi_a[6,index]
    pred_pi_a[7,index]=pred_pi_a[6,index]
    pred_pi_b[7,index]=pred_pi_a[6,index]
    pred_pi_a[9,index]=pred_pi_a[6,index]
    pred_pi_b[9,index]=pred_pi_a[6,index]
    mtry_pi=c(mtry_pi,pi_rf_TR$bestTune$mtry)
    minsize_pi=c(minsize_pi,pi_rf_TR$bestTune$min.node.size)
    pi_rf_TR=butcher(pi_rf_TR);rm(pi_rf_TR,prediction_pi_TR);gc()
  
    #p using Random forest-O-TR
    lasso_p <- lasso_var_selection(S[index_nk, , drop = FALSE],A[index_nk],task=1)
    selected_vars_p <- lasso_p$selected_vars

        p_rf_TR = caret::train(
      x = S[index_nk, , drop = FALSE][, selected_vars_p, drop = FALSE],
      y = A_fac[index_nk],
      method = "ranger",
      trControl = ctrl_pi,
      tuneGrid = generate_tune_grid_otr_p(
        length(index_nk),
        length(selected_vars_p)
      ),
      metric = "ROC",
      importance = "none",
      num.trees = generate_num_tree(length(index_nk))
    )
      
    pred_p[4,index]=pmax(pmin(predict(p_rf_TR$finalModel, data=S[index, , drop = FALSE][,selected_vars_p,drop=FALSE])$predictions[, 2],0.99), 0.01)
    pred_q[7,index]=(1-pred_p[4,index])/pred_p[4,index]*pred_pi_a[7,index]/(1- pred_pi_a[7,index])
    pred_p[6,index]=pred_p[4,index]
    pred_q[9,index]=pred_q[7,index]
    mtry_p=c(mtry_p,p_rf_TR$bestTune$mtry)
    minsize_p=c(minsize_p,p_rf_TR$bestTune$min.node.size)
    p_rf_TR=butcher(p_rf_TR);rm(p_rf_TR)
  
    #mu for O-TR
    lasso_mu<-lasso_var_selection(S[index_nk1_A, , drop = FALSE], Y[index_nk1_A], task = 0)
    selected_vars_mu<-lasso_mu$selected_vars
    mu_rf_TR=caret::train(x=S[index_nk1_A, , drop = FALSE][,selected_vars_mu,drop=FALSE], y=Y[index_nk1_A], method= "ranger", trControl=ctrl_mu, importance="none", num.trees=generate_num_tree(length(index_nk1_A)) , tuneGrid=generate_tune_grid(length(index_nk1_A),length(selected_vars_mu),0))
    mtry_mu=c(mtry_mu,mu_rf_TR$bestTune$mtry)
    minsize_mu=c(minsize_mu,mu_rf_TR$bestTune$min.node.size)
    pred_mu[7,index]=predict(mu_rf_TR$finalModel, data=S[index, , drop = FALSE][,selected_vars_mu,drop=FALSE])$predictions
  
    #tau1 for O-TR
    prediction_mu_new1_TR=predict(mu_rf_TR$finalModel, data=S[index_nk2_A0, , drop = FALSE][,selected_vars_mu,drop=FALSE])
    lasso_tau1<-lasso_var_selection(X[index_nk2_A0,], prediction_mu_new1_TR$predictions, task = 0)
    selected_vars_tau1<-lasso_tau1$selected_vars
    tau1_rf_TR=caret::train(x=X[index_nk2_A0,selected_vars_tau1,drop=FALSE], y=prediction_mu_new1_TR$predictions, method= "ranger", trControl=ctrl_mu, tuneGrid=generate_tune_grid(length(index_nk2_A0),length(selected_vars_tau1),0), importance="none", num.trees=generate_num_tree(length(index_nk2_A0)) )
    pred_tau1[7,index]=predict(tau1_rf_TR$finalModel, data=X[index, selected_vars_tau1,drop=FALSE])$predictions
    mtry_tau1=c(mtry_tau1,tau1_rf_TR$bestTune$mtry)
    minsize_tau1=c(minsize_tau1,tau1_rf_TR$bestTune$min.node.size)
  
    #mu for TR/O-TR2
    lasso_mu<-lasso_var_selection(S[index_nk_A, , drop = FALSE], Y[index_nk_A], task = 0)
    selected_vars_mu<-lasso_mu$selected_vars
    mu_rf_TR=caret::train(x=S[index_nk_A, , drop = FALSE][,selected_vars_mu,drop=FALSE], y=Y[index_nk_A], method= "ranger", trControl=ctrl_mu, importance="none", num.trees=generate_num_tree(length(index_nk_A)), tuneGrid=generate_tune_grid(length(index_nk_A),length(selected_vars_mu),0))
    pred_mu[6,index]=predict(mu_rf_TR$finalModel, data=S[index, , drop = FALSE][,selected_vars_mu,drop=FALSE])$predictions
    pred_mu[9,index]=pred_mu[6,index]
    mtry_mu=c(mtry_mu,mu_rf_TR$bestTune$mtry)
    minsize_mu=c(minsize_mu,mu_rf_TR$bestTune$min.node.size)
    # summary(pred_mu[5,index]-mu[index])
  
    #tau1 for O-TR2
    prediction_mu_new1_TR=predict(mu_rf_TR$finalModel, data=S[index_nk_A0, , drop = FALSE][,selected_vars_mu,drop=FALSE])
    lasso_tau1<-lasso_var_selection(X[index_nk_A0,], prediction_mu_new1_TR$predictions, task = 0)
    selected_vars_tau1<-lasso_tau1$selected_vars
    tau1_rf_TR=caret::train(x=as.matrix(X[index_nk_A0,selected_vars_tau1,drop=FALSE]), y=prediction_mu_new1_TR$predictions, method= "ranger", trControl=ctrl_mu, tuneGrid=generate_tune_grid(length(index_nk_A0),length(selected_vars_tau1),0), importance="none", num.trees=generate_num_tree(length(index_nk_A0)) )
    pred_tau1[9,index]=predict(tau1_rf_TR$finalModel, data=X[index,selected_vars_tau1,drop=FALSE])$predictions
    mtry_tau1=c(mtry_tau1,tau1_rf_TR$bestTune$mtry)
    minsize_tau1=c(minsize_tau1,tau1_rf_TR$bestTune$min.node.size)
  
    # Bernoulli f_a(M|X), with the same inner-CV/OOB structure as sparse7.
    x_vars <- select_x_for_m(X[index_nk, , drop = FALSE], M[index_nk, , drop = FALSE])
    f1_fit <- fit_binary_rf_cf(
      X[index_nk_A, x_vars, drop = FALSE], M[index_nk_A, , drop = FALSE],
      K = 3, ntree_tune = 100, ntree_final = 50,
      seed = seed + 10000 * t + 100 * k + 1)
    f0_fit <- fit_binary_rf_cf(
      X[index_nk_A0, x_vars, drop = FALSE], M[index_nk_A0, , drop = FALSE],
      K = 3, ntree_tune = 100, ntree_final = 50,
      seed = seed + 10000 * t + 100 * k + 2)
    log_info(paste("TR Bernoulli RF fold", k,
                   "A1: mtry =", f1_fit$best_mtry, "nodesize =", f1_fit$best_nodesize,
                   "; A0: mtry =", f0_fit$best_mtry, "nodesize =", f0_fit$best_nodesize))

    prob1_ref <- pmax(pmin(f1_fit$oob_pred, 0.99), 0.01)
    prob0_ref <- pmax(pmin(predict(f0_fit$fit,
      data = X[index_nk_A, x_vars, drop = FALSE])$predictions[, "1"], 0.99), 0.01)
    log_q_train <- dbinom(M[index_nk_A, 1], 1, prob0_ref, log = TRUE) -
      dbinom(M[index_nk_A, 1], 1, prob1_ref, log = TRUE)
    prob1 <- pmax(pmin(predict(f1_fit$fit,
      data = X[index, x_vars, drop = FALSE])$predictions[, "1"], 0.99), 0.01)
    prob0 <- pmax(pmin(predict(f0_fit$fit,
      data = X[index, x_vars, drop = FALSE])$predictions[, "1"], 0.99), 0.01)
    log_q_test <- dbinom(M[index, 1], 1, prob0, log = TRUE) -
      dbinom(M[index, 1], 1, prob1, log = TRUE)
    pred_q[6,index] <- normalized_truncate_q_logstable(log_q_test, log_q_train)

    # Exact integration over {0,1}, using the same fitted F0 as the density ratio.
    data_m0 <- S[index, , drop = FALSE]
    data_m1 <- data_m0
    data_m0[, "M1"] <- 0
    data_m1[, "M1"] <- 1
    mu_m0 <- predict_caret_regression_final(mu_rf_TR, data_m0[, selected_vars_mu, drop = FALSE])
    mu_m1 <- predict_caret_regression_final(mu_rf_TR, data_m1[, selected_vars_mu, drop = FALSE])
    pred_tau1[6,index] <- (1 - prob0) * mu_m0 + prob0 * mu_m1
  }
  lock_obj=lock(log_file, timeout=5)
  unlock(lock_obj)
  log_info(paste("Time=", t, "TR & O-TR accomplished"))
  set.seed(seed)
  
  #-------------------MSE check-------------
  order_new=c(1,2,5,8,3,4,6,7,9,10)
  mse_pi=rowMeans((pred_pi_a-matrix(rep(pi,each=num_pi_est),ncol=N))^2)
  mse_q=rowMeans((pred_q[,which(A==1)]-matrix(rep(q[which(A==1)],each=num_q_est),nrow=num_q_est))^2)
  # mse_p=rowMeans((pred_p[,which(A==1)]-matrix(rep(p[which(A==1)],each=num_p_est),nrow=num_p_est))^2)
  mse_mu=rowMeans((pred_mu-matrix(rep(mu,each=num_mu_est),nrow=num_mu_est))^2)
  pred_tau[5:9,]=pred_tau1[5:9,]
  mse_tau=rowMeans((pred_tau-matrix(rep(tau,each=num_mu_est),nrow=num_mu_est))^2)    
  mse_tau1=rowMeans((pred_tau1-matrix(rep(tau,each=num_mu_est),nrow=num_mu_est))^2) 
  res_mu_tau=rowMeans((pred_q[,which(A==1)]*(Y[A==1]-pred_mu[,which(A==1)])+pred_tau1[,which(A==1)]-pred_tau[,which(A==1)]))
  mse_pi=mse_pi[order_new];mse_q=mse_q[order_new];mse_mu=mse_mu[order_new];mse_tau=mse_tau[order_new];mse_tau1=mse_tau1[order_new];res_mu_tau=res_mu_tau[order_new]
  
  #estimation of theta & asymptotic variance estimator
  pred_theta[1,]=A*pred_q[1,]/(pred_pi_a[1,])*(Y-pred_mu[1,])+(1-A)/(1-pred_pi_b[1,])*(pred_mu[1,]-pred_tau1[1,])+A/(pred_pi_a[1,])*(pred_tau1[1,]-pred_tau[1,])+pred_tau[1,];#MQR1
  pred_theta[2,]=A*pred_q[2,]/(pred_pi_a[2,])*(Y-pred_mu[2,])+(1-A)/(1-pred_pi_b[2,])*(pred_mu[2,]-pred_tau1[2,])+A/(pred_pi_a[2,])*(pred_tau1[2,]-pred_tau[2,])+pred_tau[2,];#MQR2
  pred_theta[3,]=A/(1-pred_pi_a[5,])*(1-pred_p[3,])/(pred_p[3,])*(Y-pred_mu[5,])+(1-A)/(1-pred_pi_a[5,])*(pred_mu[5,]-pred_tau1[5,])+pred_tau1[5,];#O-TR' 
  pred_theta[4,]=A/(1-pred_pi_a[8,])*(1-pred_p[5,])/(pred_p[5,])*(Y-pred_mu[8,])+(1-A)/(1-pred_pi_a[8,])*(pred_mu[8,]-pred_tau1[8,])+pred_tau1[8,];#MTR2 (internal comparator) 
  pred_theta[5,]=A/(pred_pi_a[3,])*pred_q[3,]*(Y-pred_mu[3,])+(1-A)/(1-pred_pi_b[3,])*(pred_mu[3,]-pred_tau1[3,])+A/(pred_pi_a[3,])*(pred_tau1[3,]-pred_tau[3,])+pred_tau[3,];#QR1
  pred_theta[6,]=A/(pred_pi_a[4,])*pred_q[4,]*(Y-pred_mu[4,])+(1-A)/(1-pred_pi_b[4,])*(pred_mu[4,]-pred_tau1[4,])+A/(pred_pi_a[4,])*(pred_tau1[4,]-pred_tau[4,])+pred_tau[4,];#QR2                                 
  pred_theta[7,]=A*pred_q[6,]/(pred_pi_a[6,])*(Y-pred_mu[6,])+(1-A)/(1-pred_pi_a[6,])*(pred_mu[6,]-pred_tau1[6,])+pred_tau1[6,];#TR
  pred_theta[8,]=A/(1-pred_pi_a[7,])*(1-pred_p[4,])/(pred_p[4,])*(Y-pred_mu[7,])+(1-A)/(1-pred_pi_a[7,])*(pred_mu[7,]-pred_tau1[7,])+pred_tau1[7,];#O-TR
  pred_theta[9,]=A/(1-pred_pi_a[9,])*(1-pred_p[6,])/(pred_p[6,])*(Y-pred_mu[9,])+(1-A)/(1-pred_pi_a[9,])*(pred_mu[9,]-pred_tau1[9,])+pred_tau1[9,];#O-TR2
  pred_theta[10,]=A*f0/f1/(pred_pi_a[10,])*(Y-pred_mu[10,])+(1-A)/(1-pred_pi_a[10,])*(pred_mu[10,]-pred_tau[10,])+pred_tau[10,];#Oracle
  
  theta_hat=rowMeans(pred_theta);
  sigma_hat=rowMeans((pred_theta-theta_hat)^2); 
  theta_hat-theta; 2*sqrt(sigma_hat/N)*qnorm(1-0.05/2)
                                     result=c(pred_theta=theta_hat,pred_sig2=sigma_hat,mse_tau,mse_mu,mse_q,mse_pi,mse_tau1,res_mu_tau,sel.rho_pi_a,sel.rho_pi_b,sel.rho_q,sel.rho_p,lambda_mu,lambda_tau1,lambda_tau,mtry_pi,mtry_p,mtry_mu,mtry_tau1,mtry_tau,minsize_pi,minsize_p,minsize_mu,minsize_tau1,minsize_tau)
  # result_list=list(result=result)                               
  temp_file <- tempfile(pattern = "temp_", fileext = ".rds")
  saveRDS(result, file = temp_file)
  file.rename(temp_file, file.path(res_dir, paste0("LD1_", current_time, "_", t, ".rds")))
  lock_obj=lock(log_file, timeout=5)
  unlock(lock_obj)
  log_info(paste("Time=", t ,"saved"))                                           
  
  c(pred_theta=theta_hat,pred_sig2=sigma_hat,mse_tau,mse_mu,mse_q,mse_pi,mse_tau1,res_mu_tau,sel.rho_pi_a,sel.rho_pi_b,sel.rho_q,sel.rho_p,lambda_mu,lambda_tau1,lambda_tau,mtry_pi,mtry_p,mtry_mu,mtry_tau1,mtry_tau,minsize_pi,minsize_p,minsize_mu,minsize_tau1,minsize_tau)
}

stopCluster(cl)
pred_theta=results_par[,1:num_totest]
bias=apply(pred_theta-theta,2,median)
pred_sig2=results_par[,(num_totest+1):(2*num_totest)];   
sig2=apply(pred_sig2,2,median)

RMSE=apply((pred_theta-theta)^2,2,mean)^0.5
AL=apply(2*sqrt(pred_sig2/N)*qnorm(1-0.05/2),2,median)
AC=colMeans(abs(pred_theta-theta)<=sqrt(pred_sig2/N)*qnorm(1-0.05/2))
ESD=apply(pred_theta,2,function(x){1.4826*median(abs(x-median(x)))})
ASD=apply(pred_sig2^0.5,2,median)/sqrt(N)                                    
mse_tau=apply(results_par[,(2*num_totest+1):(2*num_totest+num_mu_est)],2,median)
mse_mu=apply(results_par[,(2*num_totest+num_mu_est+1):(2*num_totest+2*num_mu_est)],2,median)
mse_q=apply(results_par[,(2*num_totest+2*num_mu_est+1):(2*num_totest+2*num_mu_est+num_q_est)],2,median)
mse_pi=apply(results_par[,(2*num_totest+2*num_mu_est+num_q_est+1):(2*num_totest+2*num_mu_est+num_q_est+num_pi_est)],2,median)
mse_tau1=apply(results_par[,(2*num_totest+2*num_mu_est+num_q_est+num_pi_est+1):(2*num_totest+3*num_mu_est+num_q_est+num_pi_est)],2,median)
res_mu_tau=apply(results_par[,(2*num_totest+3*num_mu_est+num_q_est+num_pi_est+1):(2*num_totest+4*num_mu_est+num_q_est+num_pi_est)],2,median)
# mse_tau1=apply(results_par[,(2*num_totest+2*num_mu_est+num_q_est+num_pi_est+47*Kmb+1):(2*num_totest+3*num_mu_est+num_q_est+num_pi_est+47*Kmb)],2,median)
# res_mu_tau=apply(results_par[,(2*num_totest+3*num_mu_est+num_q_est+num_pi_est+47*Kmb+1):(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+47*Kmb)],2,median)
# order_new=c(1,2,5,8,3,4,6,7,9,10)
# mse_pi=mse_pi[order_new];mse_q=mse_q[order_new];mse_mu=mse_mu[order_new];mse_tau=mse_tau[order_new];mse_tau1=mse_tau1[order_new];res_mu_tau=res_mu_tau[order_new]


file_name <- paste0("Result_", current_time, "_", N, "_", d1, d2, "_", tmax, ".RData"); 
cat("bias=\n",paste(capture.output(bias), collapse = "\n"), file = log_file, append = TRUE)
cat("\n RMSE=\n",paste(capture.output(RMSE), collapse = "\n"), file = log_file, append = TRUE)
cat("\n AC=\n",paste(capture.output(AC), collapse = "\n"), file = log_file, append = TRUE)     
cat("\n AL=\n",paste(capture.output(AL), collapse = "\n"), file = log_file, append = TRUE) 
log_info(paste("\n", "Important data is saved in", file_name,"\n"))                              
cat(theta,"\n",bias,"\n",RMSE,"\n",AC,"\n",sig2,"\n",sum(A_all==1),"\n", sum(A_all==0),"\n")  
save(seed,num_totest,d1,d2,N,scale_sub, scale_mix, func_link, results_par, theta, pred_theta, pred_sig2, bias,RMSE, AL,AC, ESD,ASD, mse_pi, mse_q, mse_mu, mse_tau, mse_tau1,res_mu_tau, file = file_name) 
#  range_lambda_pi_a, range_lambda_pi_b, range_lambda_q, range_lambda_p, range_lambda_mu, range_lambda_tau1, range_lambda_tau,
cat("数据已保存到文件:", file_name, "\n")                                    

# range_lambda_pi_a=range(results_par[,(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+1):(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+3*Kmb)])    
# range_lambda_pi_b=range(results_par[,(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+1+3*Kmb):(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+6*Kmb)])
# range_lambda_q=range(results_par[,(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+1+6*Kmb):(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+8*Kmb)])
# range_lambda_p=range(results_par[,(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+1+8*Kmb):(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+9*Kmb)])
# range_lambda_mu=range(results_par[,(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+1+9*Kmb):(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+12*Kmb)])
# range_lambda_tau1=range(results_par[,(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+1+12*Kmb):(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+15*Kmb)])
# range_lambda_tau=range(results_par[,(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+1+15*Kmb):(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+17*Kmb)])
# range_mtry_pi=range(results_par[,(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+17*Kmb+1):(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+20*Kmb)])
# range_mtry_p=range(results_par[,(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+20*Kmb+1):(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+23*Kmb)])
# range_mtry_mu=range(results_par[,(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+23*Kmb+1):(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+27*Kmb)])
# range_mtry_tau1=range(results_par[,(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+27*Kmb+1):(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+30*Kmb)])
# range_mtry_tau=range(results_par[,(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+30*Kmb+1):(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+32*Kmb)])
# range_minsize_pi=range(results_par[,(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+32*Kmb+1):(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+35*Kmb)])
# range_minsize_p=range(results_par[,(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+35*Kmb+1):(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+38*Kmb)])
# range_minsize_mu=range(results_par[,(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+38*Kmb+1):(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+42*Kmb)])
# range_minsize_tau1=range(results_par[,(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+42*Kmb+1):(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+45*Kmb)])
# range_minsize_tau=range(results_par[,(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+45*Kmb+1):(2*num_totest+4*num_mu_est+num_q_est+num_pi_est+47*Kmb)])

# est=c("MQR1","MQR2","MR1","MR2","QR1","QR2","TR","O-TR","O-TR2","Oracle")
# table_result=function(x){
#   paste(est[x],"&",format(round(bias[x],digits=3),nsmall=3,scientific=FALSE),"&",format(round(RMSE[x],digits=3),nsmall=3,scientific=FALSE),"&",
#         format(round(AL[x],digits=3),nsmall=3,scientific=FALSE),"&",format(round(AC[x],digits=2),nsmall=2,scientific=FALSE),"\\",sep="")
# }
# result=apply(t(1:num_totest),2,table_result)[c(10,7,8,9,5,6,3,4,1,2)]
# result=list(result=c("estimator&Bias&RMSE&AL&AC",result),est=est,
#             mse_prop_pi=round(mse_pi/as.numeric(var(pi_all)),digits=3),
#             mse_prop_q=round(mse_q/var(q_all[which(A_all==1)]),digits=3),
#             mse_prop_mu=round(mse_mu/var(mu_all[which(A_all==1)]),digits=3),
#             mse_prop_tau=round(mse_tau/var(tau_all[which(A_all==1)]),digits=3))
# result
# # result
# proc.time()-ptm
