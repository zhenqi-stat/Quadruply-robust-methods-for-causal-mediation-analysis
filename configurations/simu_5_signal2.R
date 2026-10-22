# simu_5_signal2: reconstructed from the matched result's saved source.
# All reported methods are enabled; deprecated fits remain internal for RNG continuity.
# Only switches, replication-range plumbing, RNG preservation and output metadata are normalized.
rm(list=ls())
Sys.unsetenv("OMP_NUM_THREADS")
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)
suppressPackageStartupMessages({
  require(glmnet);require(itertools);require(foreach);require(doParallel);require(RCAL);require(truncnorm);require(mvtnorm);require(corpcor);require(MASS);require(matrixStats);require(ranger);require(caret);require(cobalt);require(densratio);require(butcher);require(randomForestSRC);require(logger);require(filelock);require(sn);require(FNN)
})
ptm=proc.time()
simulation_id <- "decay"
current_time <- paste0(format(Sys.time(), "%Y%m%d_%H_%M_%S"), "_simu_5_signal2")
# current_time <- "20260625_19_20"

N=1800;d1=100L;d2=50L #sample size & dimensions
n_mc=50L;
Setting="a" #choose Setting a or b

tmax=100;Kmb=5;value_trunc=1.2;n=floor(N/Kmb)
seed=123;limit=30; signal=2
current_time <- paste0(format(Sys.time(), "%Y%m%d_%H_%M_%S"), "_simu_5_signal2")
set.seed(seed)
sdr_ridge_grid <- c(0, 0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1, 3, 10)
#--------data generating--------------
# if (Setting=="a")
# Setting a: non-logistic pi, nonlinear mu
# func_link=function(x){(abs(x)+sin(x)+0.1)/(abs(x)+sin(x)+0.7)}
func_link=function(x){plogis(x)}
# func_link=function(x){0.2+0.6/(1+abs(x)^3)} #non-logistic function
#func_link=function(x){pnorm(x)-1/2}
beta_pi_i=0.36
beta_pi_X = c(rep(0.55, 3), rep(0, 57), rep(0, d1 - 60))
beta_pi_h = c(-0.2, -0.2, 0.1, 0.1)


beta_mu_i=1
beta_mu_X=c(rep(0.2, 5),rep(0,d1-5))
# beta_mu_M = c(rep(0.60, 3), rep(0, d2 - 3))
beta_mu_M = c(1,1,1,rep(0, 7), rep(0, d2 - 10))
# beta_mu_M=c(rep(0.65,3),rep(0, 7),rep(0,d2-10))
beta_m0_h = c(-0.4, 0.4, 0.4, -0.4)
# beta_m0_h=c(0.15, -0.12, 0.10)
# Local nonlinear alternative for mu.  The nonlinear/interaction component
# vanishes at the N^{-mu_nonlinear_decay_rate} rate, so the linear working model becomes
# asymptotically correct while mu remains nonlinear at every finite N.
N_ref_mu_nonlinear <- 600L
mu_nonlinear_decay_rate <- 3
beta_mu_h_at_N_ref <- c(
  X1_sq = 1,
  X1_M1 = -1,
  X2_M2 = 1,
  M1_M2 = 1
)
beta_mu_h <- beta_mu_h_at_N_ref * signal
 
epsilon_sd = 1

#X_all follows the normal distribution. 
X_all=matrix(rtruncnorm(N*d1*tmax, -value_trunc, value_trunc, mean=0, sd=1),nrow=N*tmax)
# X_all=matrix(rnorm(N*d1*tmax,mean=0, sd=1),nrow=N*tmax)
#X_all=matrix(rbeta(N*d1*tmax,1,2),nrow=N*tmax)-1/3
#X_all <- matrix(rsn(N*d1*tmax, xi = 0, omega = 1.3, alpha = 5),nrow=N*tmax)-0-1.3*5/sqrt(1 + 5^2)*sqrt(2 / pi)

hX_all=cbind(
  X1_sq=X_all[,1]^2,
  X2_sq=X_all[,2]^2,
  X1_X2=X_all[,1]*X_all[,2],
  X3_X4=X_all[,3]*X_all[,4]
)

#pi and A
pi_all=func_link(beta_pi_i+X_all%*%beta_pi_X+hX_all%*%beta_pi_h)
A_all=rbinom(N*tmax,1,pi_all)

# delta_strength <- sqrt(1.08 / 10)
# delta1_all <- c(
#   rep(delta_strength, 10),
#   rep(0, d2 - 10)
# )
delta1_all=c(1, 0.4, 0.4, rep(0,2), rep(0,d2-5))
# delta1_all = c(0.22, 0.20, 0.20, rep(0, d2 - 3))
# delta1_all = c(0.30, 0.25, 0.25, rep(0, d2 - 3))
delta2_all = c(0.4, 0.4, 0.4, rep(0, d2 - 3))
# delta1_all = c(0.2, 0.2, 0.2, rep(0.0,7), rep(0, d2 - 10))
# delta2_all = c(rep(0.4,3),rep(0,d2-3))
delta2=delta2_all[1]
delta3_all=c(rep(0.3,10),rep(0,d2-10))
delta4_all=c(rep(0.3,10),rep(0,d2-10))
delta5_all_i=-1*1
delta5_all_X=0*beta_mu_X
delta5_all_M=0*beta_mu_M

#Phi
Phi=matrix(0,nrow=d1,ncol=d2)
for (i in 1:d1){
  for (j in 1:d2){
    if (i==j){
      Phi[i,j]=1 
    } else if (abs(i-j)==1){
      Phi[i,j]=0.5
    } else if (abs(i-j)==2){
      Phi[i,j]=0
    }
    # Phi[i,j]=0.35^(abs(i-j))
  }
}

# Dense common mediator means; preserve M1:M3 and the sum of M4:M10.
D_common <- matrix(0,nrow=d1,ncol=d2)
cycle_loadings <- c(0.10,0.10,0.10,0.10,-0.10,-0.10,-0.20)
for (i in seq_len(d1)) {
  D_common[i,4:10] <- cycle_loadings[((seq_len(7)+i-2L) %% 7L)+1L]
}
H_sign <- matrix(1,nrow=1,ncol=1)
while (nrow(H_sign)<max(d1,d2)) {
  H_sign <- rbind(cbind(H_sign,H_sign),cbind(H_sign,-H_sign))
}
D_common[,11:d2] <- 0.10*H_sign[seq_len(d1),2:(d2-9),drop=FALSE]
Phi <- Phi+2*D_common

#Sigma & U
Sigma=diag(d2)
# Sigma[1,2]=0.3; Sigma[2,1]=0.3
# Sigma[1,3]=0.5; Sigma[3,1]=0.5
for (i in 1:d2){
  for (j in 1:d2){
    if (abs(i-j)==1){
      Sigma[i,j]=0
    } else if (abs(i-j)==2){
      Sigma[i,j]=0
    }
  }
}
# Treatment-specific error variances, following jasa_nonexp4.
# Preserve the original correlation structure via diagonal rescaling.
s_q <- min(3L, d2)
weight0 <- c(rep(0.4, s_q), rep(1, d2-s_q))
weight1 <- c(rep(0.6, s_q), rep(1, d2-s_q))
Sigma0 <- diag(sqrt(weight0)) %*% Sigma %*% diag(sqrt(weight0))
Sigma1 <- diag(sqrt(weight1)) %*% Sigma %*% diag(sqrt(weight1))
U_all <- matrix(0, nrow=N*tmax, ncol=d2)
for (a in 0:1) {
  idx_a <- which(A_all == a)
  if (length(idx_a) > 0L) {
    U_all[idx_a,] <- rmvnorm(length(idx_a), mean=rep(0,d2),
                            sigma=if (a == 0) Sigma0 else Sigma1)
  }
}

#M_all
scale_XM=0.5
scale_Me=1
delta1_all_matrix=matrix(nrow=N*tmax,ncol=d2)
delta2_all_matrix=matrix(nrow=N*tmax,ncol=d2)
for (k in 1:d2){delta1_all_matrix[,k]=delta1_all[k];delta2_all_matrix[,k]=delta2_all[k]}
m0_all=scale_XM*(X_all)%*%Phi
m0_all[,1]=m0_all[,1]+beta_m0_h[1]*hX_all[,1]
m0_all[,2]=m0_all[,2]+beta_m0_h[2]*hX_all[,2]
m0_all[,3]=m0_all[,3]+beta_m0_h[3]*hX_all[,3]
m0_all[,4]=m0_all[,4]+beta_m0_h[4]*hX_all[,4]
M_all=m0_all+scale_Me*U_all+delta2_all_matrix+A_all*delta1_all_matrix
U_tot_all=scale_Me*U_all+delta2_all_matrix+A_all*delta1_all_matrix

# q = f(M | A=0,X) / f(M | A=1,X).
# Unequal covariances retain quadratic M terms and M-by-m0(X) interactions.
# Compute on the log scale to avoid dividing tiny high-dimensional densities.
U_1_all=U_all-(1-A_all)*(delta1_all_matrix)/scale_Me
U_0_all=U_all+A_all*(delta1_all_matrix)/scale_Me
log_f_1_all=dmvnorm(U_1_all, mean=rep(0,d2), sigma=Sigma1, log=TRUE)-
  d2*log(scale_Me)
log_f_0_all=dmvnorm(U_0_all, mean=rep(0,d2), sigma=Sigma0, log=TRUE)-
  d2*log(scale_Me)
log_q_all=log_f_0_all-log_f_1_all
q_all=exp(log_q_all)
f_1_all=exp(log_f_1_all)
f_0_all=exp(log_f_0_all)

# p = P(A=1 | M,X)
p_all=plogis(qlogis(pi_all)-log_q_all)

scale_mix=0
scale_sub=0
scale_mue=0
scale_exp=0
scale_xsq=0

#mu_all
mu_all=beta_mu_i+M_all%*%beta_mu_M+X_all%*%beta_mu_X+
  beta_mu_h["X1_sq"]*hX_all[,1]+
  beta_mu_h["X1_M1"]*X_all[,1]*M_all[,1]+
  beta_mu_h["X2_M2"]*X_all[,2]*M_all[,2]+
  beta_mu_h["M1_M2"]*M_all[,1]*M_all[,2]
# scale_mix* (U_tot_all[,s_nl:1]*X_all[,1:s_nl])%*%beta_mu_nl[1:s_nl]+scale_sub*((U_tot_all[,1]>=delta2_all[1]))*(X_all%*%beta_mu_X)+scale_xsq*rowSums((X_all[,1:3])^2*(U_tot_all[,1:3]-delta2_all_matrix[,1:3]))
#+exp(scale_exp*(U_tot_all)%*%delta4_all)+scale_mue*(U_tot_all)^2%*%delta3_all
#+scale_mix*((U_all[,d2]>0)*M_all[,d2-1]*X_all[,1]+(U_all[,d2-1]>0)*M_all[,d2-2]*X_all[,2]+(U_all[,d2-2]>0)*M_all[,d2]*X_all[,3])
#+(X_all%*%Phi*U_all)%*%beta_mu_mix+scale_abs*abs(U_all)%*%delta4_all

#tau_all
bsm_all=m0_all+delta2_all_matrix
#tau1_all & tau_all
tau_all=beta_mu_i+bsm_all%*%beta_mu_M+X_all%*%beta_mu_X+
  beta_mu_h["X1_sq"]*hX_all[,1]+
  beta_mu_h["X1_M1"]*X_all[,1]*bsm_all[,1]+
  beta_mu_h["X2_M2"]*X_all[,2]*bsm_all[,2]+
  beta_mu_h["M1_M2"]*(bsm_all[,1]*bsm_all[,2]+scale_Me^2*Sigma0[1,2])+scale_mue*as.vector((diag(Sigma)+delta2_all^2)%*%delta3_all)+scale_mix*(delta2)*(X_all%*%beta_mu_X)+scale_sub*(1-pnorm(0))*(X_all%*%beta_mu_X)
# +scale_mue*as.vector((diag(Sigma)+delta2_all^2)%*%delta3_all)+scale_mix*(delta2_all[s_nl:1]*X_all[,1:s_nl])%*%beta_mu_nl[1:s_nl]+scale_sub*(1-pnorm(0))*(X_all%*%beta_mu_X)
#+as.vector(exp(scale_exp*delta4_all%*%delta2_all+0.5*scale_exp^2*scale_Me^2*delta4_all%*%(Sigma%*%delta4_all)))
#+scale_mix/2*(delta2)*(X_all[,1]+X_all[,2]+X_all[,3])
#+scale_abs*sum(sqrt(2/base::pi)*delta4_all)
#+scale_sub*0*(3-pnorm(delta1_all[1])-pnorm(delta1_all[2])-pnorm(delta1_all[3]))*X_all%*%beta_mu_X
tau1_all=tau_all

#theta
theta=mean(tau_all)

#Y
epsilon=rnorm(N*tmax,mean=0,sd=epsilon_sd)
Y_all=mu_all+(1-A_all)*(delta5_all_i+X_all%*%delta5_all_X+M_all%*%delta5_all_M)+epsilon

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
ctrl_pi=trainControl(method="cv", number=5, summaryFunction=twoClassSummary, classProbs=TRUE)
ctrl_mu=trainControl(method="cv", number=5, summaryFunction=defaultSummary,classProbs=FALSE,
                     savePredictions=FALSE)

# 数据标准化函数（中心化 + 缩放）
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
    rf <- ranger(y ~ ., data = data.frame(X=X, y=y_current), importance = "impurity", num.trees=100)
    importance<-importance(rf)
    selected_vars <- names(importance)[importance > mean(importance)+sd(importance)]#
    # if(length(selected_vars)==0){
    #   rf <- ranger(y ~ ., data = data.frame(X=X, y=y_current), importance = "impurity", num.trees=50)
    #   importance<-importance(rf)
    #   selected_vars <- names(importance)[importance > mean(importance)+1/2*sd(importance)]}
    # if(length(selected_vars)==0){
    #   rf <- ranger(y ~ ., data = data.frame(X=X, y=y_current), importance = "impurity", num.trees=50)
    #   importance<-importance(rf)
    #   selected_vars <- names(importance)[importance > mean(importance)]}
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
    importance = FALSE,
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



el_drm_raw_q <- function(M) {
  M <- as.matrix(M)
  
  if (is.null(colnames(M))) {
    colnames(M) <- paste0("M", seq_len(ncol(M)))
  }
  
  d <- ncol(M)
  
  # First-order terms
  linear_terms <- M
  colnames(linear_terms) <- paste0("L_", colnames(M))
  
  # Squares
  square_terms <- M^2
  colnames(square_terms) <- paste0("Q_", colnames(M), "_sq")
  
  # Pairwise interactions
  if (d >= 2L) {
    pairs <- combn(seq_len(d), 2)
    
    interaction_terms <- sapply(seq_len(ncol(pairs)), function(k) {
      M[, pairs[1, k]] * M[, pairs[2, k]]
    })
    
    if (is.null(dim(interaction_terms))) {
      interaction_terms <- matrix(interaction_terms, ncol = 1)
    }
    
    colnames(interaction_terms) <- apply(pairs, 2, function(idx) {
      paste0(
        "Q_",
        colnames(M)[idx[1]],
        "_",
        colnames(M)[idx[2]]
      )
    })
    
    cbind(
      linear_terms,
      square_terms,
      interaction_terms
    )
  } else {
    cbind(
      linear_terms,
      square_terms
    )
  }
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

select_el_drm_basis <- function(
    X, M, A,
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
  
  linear_names <- grep(
    "^L_",
    colnames(full_basis),
    value = TRUE
  )
  
  quadratic_names <- grep(
    "^Q_",
    colnames(full_basis),
    value = TRUE
  )
  
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
        type.measure = "deviance"
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

#-------------experimental log----------------------
log_dir <- "logs"
if (!dir.exists(log_dir)) {
  dir.create(log_dir)
}
log_file<-file.path(log_dir, paste0("Log_", current_time, ".txt"))
res_dir <- "res_1"
if (!dir.exists(res_dir)) {
  dir.create(res_dir)
}
log_appender(appender_file(log_file))
log_formatter(formatter_paste)
log_threshold(INFO)  # 设置日志级别为 INFO
log_info("Logger initialized and ready to record logs.")
# Store the exact same parameter snapshot in the text log and RData.
# dput preserves vector names, matrix dimensions and numeric precision.
dgp_config <- list(
  simulation_id=simulation_id, run_id=current_time,
  N=N, d1=d1, d2=d2, tmax=tmax, seed=seed, Setting=Setting,
  Kmb=Kmb, n_mc=n_mc, value_trunc=value_trunc,
  beta_pi_i=beta_pi_i, beta_pi_X=beta_pi_X, beta_pi_h=beta_pi_h,
  func_link=paste(deparse(func_link), collapse="
"),
  hX_basis=c("X1^2", "X2^2", "X1*X2", "X3*X4"),
  beta_m0_h=beta_m0_h, scale_XM=scale_XM, scale_Me=scale_Me,
  delta1_all=delta1_all, delta2_all=delta2_all,
  q_definition="f(M|A=0,X)/f(M|A=1,X)",
  q_design="Gaussian, treatment-specific covariance; nonlinear m0(X)",
  s_q=s_q, weight0=weight0, weight1=weight1,
  Sigma=Sigma, Sigma0=Sigma0, Sigma1=Sigma1,
  covariance_definition="Cov(M|A=a,X)=scale_Me^2*Sigma_a; Sigma_a=D_a*Sigma*D_a, D_a=diag(sqrt(weight_a))",
  beta_mu_i=beta_mu_i, beta_mu_X=beta_mu_X, beta_mu_M=beta_mu_M,
  mu_basis=c("X1^2", "X1*M1", "X2*M2", "M1*M2"),
  N_ref_mu_nonlinear=N_ref_mu_nonlinear,
  mu_nonlinear_decay_rate=mu_nonlinear_decay_rate,
  beta_mu_h_at_N_ref=beta_mu_h_at_N_ref,
  mu_decay_multiplier=(N_ref_mu_nonlinear/N)^mu_nonlinear_decay_rate,
  beta_mu_h=beta_mu_h, epsilon_sd=epsilon_sd,
  mu_decay_definition="beta_mu_h=beta_mu_h_at_N_ref*(N_ref_mu_nonlinear/N)^mu_nonlinear_decay_rate",
  tau_covariance="scale_Me^2*Sigma0[1,2] for E(M1*M2|A=0,X)",
  delta5_all_i=delta5_all_i, delta5_all_X=delta5_all_X, delta5_all_M=delta5_all_M
)
cat("
DGP configuration:
",
    paste(capture.output(dput(dgp_config, control=c("keepNA", "keepInteger", "niceNames", "showAttributes", "digits17"))), collapse="
"), "
",
    file=log_file, append=TRUE)
for (label in c("pi_all", "p_all", "log_q_all")) {
  cat("
summary(", label, "):
",
      paste(capture.output(summary(get(label))), collapse="
"), "
",
      sep="", file=log_file, append=TRUE)
}
for (a in 0:1) {
  cat("
summary(log_q_all | A=", a, "), n=", sum(A_all == a), ":
",
      paste(capture.output(summary(log_q_all[A_all == a])), collapse="
"), "
",
      sep="", file=log_file, append=TRUE)
}
cat("
summary(1/(1+q_all)):
",
    paste(capture.output(summary(plogis(-log_q_all))), collapse="
"), "
",
    file=log_file, append=TRUE)

predict_caret_regression_final <- function(fit, newdata) {
  pred <- predict(fit$finalModel, data = newdata)
  if (!is.null(pred$predictions)) {
    as.numeric(pred$predictions)
  } else {
    as.numeric(predict(fit, newdata = newdata))
  }
}

#----------Computing estimators-------
#num_cores <- detectCores()
#print(paste("Available cores:", num_cores))
if (!is.null(getOption("cl"))) { stopCluster(getOption("cl")) }
cl <- makeCluster(min(limit,tmax))  # 使用核心数减 1，保留一个核心给操作系统
registerDoParallel(cl)
cat("Simulation begun\n")

num_mu_est=10;num_pi_est=10;num_q_est=10;num_p_est=6;num_totest=num_pi_est
grid=data.frame(est_mu=1:num_mu_est,est_pi=1:num_pi_est)
t=1;k=1

replication_ids <- 1:100
results_par=foreach (t=replication_ids, .combine=rbind, .packages=c("ranger", "MASS", "matrixStats","glmnet","RCAL", "mvtnorm", "caret", "densratio","randomForestSRC","butcher","logger","filelock","corpcor","FNN")) %dopar% {
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
  
  X=X_all[(t-1)*N+1:N,];M=M_all[(t-1)*N+1:N,];colnames(X)=paste0("X", 1:ncol(X));colnames(M)=paste0("M", 1:ncol(M))
  true_q_M <- sort(unique(which(delta1_all != 0)))
  # true_q_X <- q_density_x_vars(X, true_q_M)
  true_q_X <- multi_y_var_selection(X, M[, true_q_M])
  true_mu_M <- sort(unique(which(beta_mu_M != 0)))
  A=A_all[(t-1)*N+1:N];Y=Y_all[(t-1)*N+1:N]; A_fac=factor(A, levels=c(0,1), labels=c("Control", "Treatment"))
  pi=pi_all[(t-1)*N+1:N];
  p=p_all[(t-1)*N+1:N];
  q=q_all[(t-1)*N+1:N];
  f0=f_0_all[(t-1)*N+1:N];#可以删掉
  f1=f_1_all[(t-1)*N+1:N];
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
  set.seed(seed)  
  
  log_appender(appender_file(log_file))
  log_formatter(formatter_paste)
  log_threshold(INFO)  
  
  #estimation of pi, p, q, mu, tau1, tau
  
  # Retain the original pre-method RNG consumption and method call order.
  archive_qr1_entry_rng <- get(".Random.seed", envir = .GlobalEnv)
  # ------------MQR1-----------
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
  
    # max(abs(beta_q[1]+cbind(X[index,],M[index,])%*%beta_q[-1]-beta_q_i-X[index,]%*%beta_q_X-M[index,]%*%beta_q_M))
    # grad_q=colSums(as.numeric(A[index_train2]/func_link(alpha_a[1]+X[index_train2,]%*%alpha_a[-1])*exp(beta_q[1]+cbind(X[index_train2,],M[index_train2,])%*%beta_q[-1])-(1-A[index_train2])/(1-func_link(alpha_b[1]+X[index_train2,]%*%alpha_b[-1])))*cbind(X[index_train1,],M[index_train1,]))/length(index_train2)
    # grad_q_phi=colSums(as.numeric(A[index_train2]/func_link(alpha_a[1]+X[index_train2,]%*%alpha_a[-1]))*X[index_train2,]%*%Phi-as.numeric((1-A[index_train2])/(1-func_link(alpha_b[1]+X[index_train2,]%*%alpha_b[-1])))*M[index_train2,])/length(index_train2)
  
    #mu
    # lambda_mu=c()
    fit_lasso_mu=cv.glmnet(x=S_poly[index_train3,],y=Y[index_train3],weights=A[index_train3]*exp(beta_q[1]+S_poly[index_train3,]%*%beta_q[-1])*(1+exp(-alpha_a[1]-X_poly[index_train3,]%*%alpha_a[-1])),family="gaussian",alpha=1)
    pred_mu[1,index]=predict(fit_lasso_mu,newx=S_poly[index,],s="lambda.min")
    lambda_mu=c(lambda_mu,fit_lasso_mu$lambda.min/sqrt(log(d1+d2)/length(index_train3)))
    # beta_mu=as.vector(coef(fit_lasso_mu$glmnet.fit, s=fit_lasso_mu$lambda.min))
    # summary(pred_mu[1,index]-mu[index])
  
    # summary(A[index_A]*q[index_A]/(pi[index_A])*(Y[index_A]-mu[index_A]))
    # summary(A[index_A]*q[index_A]/(pi[index_A])*(Y[index_A]-mu[index_A])-A[index_A]*pred_q[1,index_A]/(pred_pi_a[1,index_A])*(Y[index_A]-pred_mu[1,index_A]))
    # summary(A[index_A]*pred_q[1,index_A]/(pred_pi_a[1,index_A])*(Y[index_A]-pred_mu[1,index_A]))
  
    #tau1
    index_train4_A0=index_train4[which(A[index_train4]==0)]
    find_mu=which(fit_lasso_mu$lambda==fit_lasso_mu$lambda.min)
    a0_mu=fit_lasso_mu$glmnet.fit$a0[find_mu]
    alpha_mu_X=fit_lasso_mu$glmnet.fit$beta[1:d1,find_mu]
    alpha_mu_M=fit_lasso_mu$glmnet.fit$beta[(d1+1):(d1+d2),find_mu]
    pred_mu_new1=S_poly[index_train4_A0,(d1+1):(d1+d2)]%*%alpha_mu_M+a0_mu
    #pred_mu_new1=predict(fit_lasso_mu,newx=S_poly[index_train4,],s="lambda.min")
    if (sd(pred_mu_new1)==0 || length(pred_mu_new1)<8){pred_tau1[1,index]=X_poly[index,]%*%alpha_mu_X+
      rep(mean(pred_mu_new1),length(index)); lambda_tau1=c(lambda_tau1,-1)} else{
        fit_lasso_tau1=cv.glmnet(x=X_poly[index_train4_A0,],y=pred_mu_new1, weights=exp(alpha_b[1]+X_poly[index_train4_A0,]%*%alpha_b[-1]),family="gaussian",alpha=1)
        pred_tau1[1,index]=X_poly[index,]%*%alpha_mu_X+
          predict(fit_lasso_tau1,newx=X_poly[index,],s="lambda.min");
        lambda_tau1=c(lambda_tau1,fit_lasso_tau1$lambda.min/sqrt(log(d1)/length(index_train4)))}
    # beta_tau1=as.vector(coefficients=coef(fit_lasso_tau1$glmnet.fit, s=fit_lasso_tau1$lambda.min))
    # pred_tau1_test=pred_tau1[1,index]
  
    #tau1-gamma
    # lambda_mu=c()
    # fit_lasso_tau1=cv.glmnet(x=X[index_train4,],y=M[index_train4,]%*%beta_mu[(1+d1+1):(1+d1+d2)], weights=(1-A[index_train4])*exp(alpha_b[1]+X[index_train4,]%*%alpha_b[-1]),family="gaussian",alpha=1)
    # pred_tau1[1,index]=beta_mu[1]+X[index,]%*%beta_mu[2:(1+d1)]+predict(fit_lasso_tau1,newx=X[index,],s="lambda.min")
    # lambda_tau1=c(lambda_tau1,fit_lasso_tau1$lambda.min/sqrt(log(d1)/length(index_train4)))
  
    # check: estimation of tau1 on A==0
    # index_A0=index[A[index]==0]
    # summary(pred_mu[1,index_A0]-pred_tau1[1,index_A0])
    # summary((1-A[index_A0])/(1-pred_pi_b[1,index_A0])*(pred_mu[1,index_A0]-pred_tau1[1,index_A0]))
    # summary(A[index]*pred_q[1,index]/(pred_pi_a[1,index])*(Y[index]-pred_mu[1,index])+(1-A[index])/(1-pred_pi_b[1,index])*(pred_mu[1,index]-pred_tau1[1,index]))
    # summary(A[index]*pred_q[1,index]/(pred_pi_a[1,index])*(Y[index]-pred_mu[1,index])+(1-A[index])/(1-pred_pi_b[1,index])*(pred_mu[1,index]-pred_tau1_test))
  
    #tau
    pred_mu_new2=S_poly[index_train5,(d1+1):(d1+d2)]%*%alpha_mu_M+a0_mu
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
    #test: estimation of tau
    # summary(abs(pred_tau[1,index]-tau[index]))
  
    #test the summation of terms that should be offset to zero
    # summary(A[index_A]*pred_q[1,index_A]/(pred_pi_a[1,index_A])*(Y[index_A]-pred_mu[1,index_A])+A[index_A]/(pred_pi_a[1,index_A])*(pred_tau1[1,index_A]-pred_tau[1,index_A]))
  
    # index_max_0=index[A[index]*pred_q[1,index]/(pred_pi_a[1,index])*(Y[index]-pred_mu[1,index])+A[index]/(pred_pi_a[1,index])*(pred_tau1[1,index]-pred_tau[1,index])==max(A[index]*pred_q[1,index]/(pred_pi_a[1,index])*(Y[index]-pred_mu[1,index])+A[index]/(pred_pi_a[1,index])*(pred_tau1[1,index]-pred_tau[1,index]))];cbind(pred_pi_a[1,index_max_0],pred_pi_b[1,index_max_0],pred_q[1,index_max_0],pred_mu[1,index_max_0],pred_tau1[1,index_max_0],pred_tau[1,index_max_0]);cbind(pi[index_max_0],q[index_max_0],mu[index_max_0],tau[index_max_0])
    #test theta
    # pred_theta[1,index]=A[index]*pred_q[1,index]/(pred_pi_a[1,index])*(Y[index]-pred_mu[1,index])+(1-A[index])/(1-pred_pi_b[1,index])*(pred_mu[1,index]-pred_tau1[1,index])+A[index]/(pred_pi_a[1,index])*(pred_tau1[1,index]-pred_tau[1,index])+pred_tau[1,index];summary(pred_theta[1,index])
    # pred_theta[7,index]=A[index]*pred_q[6,index]/(pred_pi_a[7,index])*(Y[index]-pred_mu[7,index])+(1-A[index])/(1-pred_pi_a[7,index])*(pred_mu[7,index]-pred_tau[7,index])+pred_tau[7,index];
  }
  lock_obj=lock(log_file, timeout=5)
  unlock(lock_obj)
  log_info(paste("Time=", t, "MQR1 accomplished"))
  set.seed(seed)

  #------------MQR2------------------
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
  
    #For test of imputation error of q
    # grad_pi_a=colSums(as.numeric(1-A[index_nk]/plogis(alpha_a[1]+X[index_nk,]%*%alpha_a[-1]))*X[index_nk,])/length(index_nk)
    # grad_pi_b=colSums(as.numeric(1-(1-A[index_nk])/(1-plogis(alpha_b[1]+X[index_nk,]%*%alpha_b[-1])))*X[index_train1,])/length(index_nk)
    # iw=as.numeric(A[index_nk]/plogis(alpha_a[1]+X[index_nk,]%*%alpha_a[-1])-(1-A[index_nk])/(1-plogis(alpha_b[1]+X[index_nk,]%*%alpha_b[-1])))
    # mean(iw)
    # (colSums(iw*X[index_nk,])/length(index_nk))
    # iw_oracle=as.numeric(A[index_nk]/pi[index_nk]-(1-A[index_nk])/(1-pi[index_nk]))
    # mean(iw_oracle)
    # mean(colSums(iw_oracle*X[index_nk,])/length(index_nk))
  
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
  
    # grad_q=colSums(as.numeric(A[index_nk]/func_link(alpha_a[1]+X[index_nk,]%*%alpha_a[-1])*exp(beta_q[1]+cbind(X[index_nk,],M[index_nk,])%*%beta_q[-1])-(1-A[index_nk])/(1-func_link(alpha_b[1]+X[index_nk,]%*%alpha_b[-1])))*cbind(X[index_train1,],M[index_train1,]))/length(index_nk)
    # grad_q_phi=colSums(as.numeric(A[index_nk]/func_link(alpha_a[1]+X[index_nk,]%*%alpha_a[-1]))*X[index_nk,]%*%Phi-as.numeric((1-A[index_nk])/(1-func_link(alpha_b[1]+X[index_nk,]%*%alpha_b[-1])))*M[index_nk,])/length(index_nk)
  
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
    pred_mu_new1=a0_mu+S_poly[index_nk_A0,(d1+1):(d1+d2)]%*%alpha_mu_M
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
  
    # summary(1-A[index])/(1-pi[index])*(mu[index]-tau1[index])
    # check: estimation of tau1 on A==0
    # index_A0=index[A[index]==0]
    # summary(pred_mu[2,index_A0]-pred_tau1[2,index_A0])
    # summary((1-A[index_A0])/(1-pred_pi_b[2,index_A0])*(pred_mu[2,index_A0]-pred_tau1[2,index_A0]))
    # summary(A[index]*pred_q[2,index]/(pred_pi_a[2,index])*(Y[index]-pred_mu[2,index])+(1-A[index])/(1-pred_pi_b[2,index])*(pred_mu[2,index]-pred_tau1[2,index]))
    # summary(A[index]*pred_q[2,index]/(pred_pi_a[2,index])*(Y[index]-pred_mu[2,index])+(1-A[index])/(1-pred_pi_b[2,index])*(pred_mu[2,index]-pred_tau1_test))
  
    #tau
    pred_mu_new1=a0_mu+S_poly[index_nk_A,(d1+1):(d1+d2)]%*%alpha_mu_M
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
    #test: estimation of tau
    # summary(abs(pred_tau[2,index]-tau[index]))
  
    #test the summation of terms that should be offset to zero
    # summary(A[index_A]*pred_q[2,index_A]/(pred_pi_a[2,index_A])*(Y[index_A]-pred_mu[2,index_A])+A[index_A]/(pred_pi_a[2,index_A])*(pred_tau1[2,index_A]-pred_tau[2,index_A]))
  
    # index_max_0=index[A[index]*pred_q[2,index]/(pred_pi_a[2,index])*(Y[index]-pred_mu[2,index])+A[index]/(pred_pi_a[2,index])*(pred_tau1[2,index]-pred_tau[2,index])==max(A[index]*pred_q[2,index]/(pred_pi_a[2,index])*(Y[index]-pred_mu[2,index])+A[index]/(pred_pi_a[2,index])*(pred_tau1[2,index]-pred_tau[2,index]))];cbind(pred_pi_a[2,index_max_0],pred_pi_b[2,index_max_0],pred_q[2,index_max_0],pred_mu[2,index_max_0],pred_tau1[2,index_max_0],pred_tau[2,index_max_0]);cbind(pi[index_max_0],q[index_max_0],mu[index_max_0],tau[index_max_0])
    #test theta
    # pred_theta[2,index]=A[index]*pred_q[2,index]/(pred_pi_a[2,index])*(Y[index]-pred_mu[2,index])+(1-A[index])/(1-pred_pi_b[2,index])*(pred_mu[2,index]-pred_tau1[2,index])+A[index]/(pred_pi_a[2,index])*(pred_tau1[2,index]-pred_tau[2,index])+pred_tau[2,index];summary(pred_theta[2,index])
  }
  lock_obj=lock(log_file, timeout=5)
  unlock(lock_obj)
  log_info(paste("Time=", t, "MQR2 accomplished"))
  set.seed(seed)

  # Preserve the QR1 entry state of its published QR-only batch.
  assign(".Random.seed", archive_qr1_entry_rng, envir = .GlobalEnv)
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
    lasso_pi <- lasso_var_selection(X[index_nk12,],A[index_nk12],task=1)
    selected_vars_pi <- lasso_pi$selected_vars
    pi_rf=caret::train(x=X[index_nk12,selected_vars_pi,drop=FALSE], y=A_fac[index_nk12], method= "ranger", trControl=ctrl_pi, tuneGrid=generate_tune_grid(length(index_nk12),length(selected_vars_pi),1), metric="ROC", importance="none", num.trees=generate_num_tree(length(index_nk12)))#,class.weights = c("Control" = 0.5, "Treatment" = 0.5)
    prediction_pi=predict(pi_rf$finalModel, data=X[index,selected_vars_pi,drop=FALSE])
    pred_pi_a[3,index]=pmax(pmin(prediction_pi$predictions[,2], 0.99), 0.01)
    pred_pi_b[3,index]=pred_pi_a[3,index]
    mtry_pi=c(mtry_pi,pi_rf$bestTune$mtry)
    minsize_pi=c(minsize_pi,pi_rf$bestTune$min.node.size)
    
    # f semiparametrically
    selected_M <- select_different_M_components_ks(
      X[index_nk12, ],
      M[index_nk12, ],
      A[index_nk12],
      max_selected = 8,
      min_selected = 1
    )
    
    x_vars <- multi_y_var_selection(
      X[index_nk12, ],
      M[index_nk12, selected_M, drop = FALSE]
    )$union_vars
    basis_sel <- select_el_drm_basis(
      X = X[index_nk12, x_vars, drop = FALSE],
      M = M[index_nk12, selected_M, drop = FALSE],
      A = A[index_nk12],
      stability_reps = 10,
      subsample_frac = 0.8,
      cv_folds = 5,
      alpha = 0.5,
      stability_cutoff = 0.8,
      seed = seed + 20000 * t + 200 * k
    )
    
    selected_basis <- basis_sel$selected_basis
    
    log_info(paste(
      "QR1 basis fold", k,
      "; selected_M =", paste(selected_M, collapse = ","),
      "; selected_basis =", paste(selected_basis, collapse = ","),
      "; successful_reps =", basis_sel$successful_reps,
      "; quadratic_freq =",
      paste(
        names(basis_sel$quadratic_frequency),
        round(basis_sel$quadratic_frequency, 2),
        sep = ":",
        collapse = ","
      )
    ))
    
    f_el_drm_sdr1 <- el_drm_fit_density_ratio(
      X[index_nk12, x_vars, drop = FALSE],
      M[index_nk12, selected_M, drop = FALSE],
      A[index_nk12],
      basis_names = selected_basis,
      ridge_grid = sdr_ridge_grid,
      ridge_cv_folds = 5,
      ridge_cv_seed = seed + 10000 * t + 100 * k + 2
    )
    log_ratio_index <- el_drm_predict_log_ratio(
      f_el_drm_sdr1,
      X[index, x_vars, drop = FALSE],
      M[index, selected_M, drop = FALSE]
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
      "theta_dim =", length(f_el_drm_sdr1$theta)
    ))
    
    #mu
    lasso_mu<-lasso_var_selection(cbind(X[index_nk1_A,], M[index_nk1_A,]), Y[index_nk1_A], task = 0)
    selected_vars_mu<-lasso_mu$selected_vars
    mu_rf=caret::train(x=cbind(X[index_nk1_A,],M[index_nk1_A,])[,selected_vars_mu,drop=FALSE], y=Y[index_nk1_A], method= "ranger", tuneGrid=generate_tune_grid(length(index_nk1_A),length(selected_vars_mu),0), trControl=ctrl_mu, importance="none", num.trees=generate_num_tree(length(index_nk1_A)) )
    pred_mu[3,index]=predict(mu_rf$finalModel, data=cbind(X[index, ],M[index,])[,selected_vars_mu,drop=FALSE])$predictions
    mtry_mu=c(mtry_mu,mu_rf$bestTune$mtry)
    minsize_mu=c(minsize_mu,mu_rf$bestTune$min.node.size)
    # summary(pred_mu[3,index]-mu[index])
    
    #tau1
    prediction_mu_1=predict(mu_rf$finalModel, data=cbind(X[index_nk2_A0,],M[index_nk2_A0,])[,selected_vars_mu,drop=FALSE])
    lasso_tau1<-lasso_var_selection(X[index_nk2_A0,], prediction_mu_1$predictions, task = 0)
    selected_vars_tau1<-lasso_tau1$selected_vars
    tau1_rf=caret::train(x=X[index_nk2_A0,selected_vars_tau1,drop=FALSE], y=prediction_mu_1$predictions, method= "ranger", tuneGrid=generate_tune_grid(length(index_nk2_A0),length(selected_vars_tau1),0), trControl=ctrl_mu, importance="none", num.trees=generate_num_tree(length(index_nk2_A0)))
    mtry_tau1=c(mtry_tau1,tau1_rf$bestTune$mtry)
    minsize_tau1=c(minsize_tau1,tau1_rf$bestTune$min.node.size)
    pred_tau1[3,index]=predict(tau1_rf$finalModel, data=X[index,selected_vars_tau1,drop=FALSE])$predictions
    # summary(pred_tau1[3,index])
    
    #tau
    # prediction_pi_1=pmax(pmin(predict(pi_rf$finalModel, data=X[index_nk3_A,selected_vars_pi,drop=FALSE])$predictions[, 2], 0.99), 0.01)
    # prediction_p_1=pmax(pmin(predict(p_rf$finalModel, data=cbind(X[index_nk3_A,],M[index_nk3_A, ])[,selected_vars_p,drop=FALSE])$predictions[, 2], 0.99), 0.01)
    # prediction_q_1=(1-prediction_p_1)/prediction_p_1*prediction_pi_1/(1-prediction_pi_1)
    
    log_ratio_index_nk3_A <- el_drm_predict_log_ratio(f_el_drm_sdr1,
                                                      X[index_nk3_A, x_vars, drop = FALSE],
                                                      M[index_nk3_A, selected_M, drop = FALSE])
    prediction_q_1 <- normalized_truncate_q(
      log_q = log_ratio_index_nk3_A,
      ref_log_q = ref_log_q_sdr1
    )
    # prediction_q_1=exp(pmax(pmin(log_ratio_index_nk3_A, log(100)), log(0.01)))
    # prediction_q_1=exp(log_ratio_index_nk3_A)
    
    prediction_mu_2=predict(mu_rf$finalModel, data=cbind(X[index_nk3_A,],M[index_nk3_A,])[,selected_vars_mu,drop=FALSE])
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
    
    # q using Random forest
    selected_M <- select_different_M_components_ks(
      X[index_nk, ],
      M[index_nk, ],
      A[index_nk],
      max_selected = 8,
      min_selected = 1
    )
    x_vars <- multi_y_var_selection(X[index_nk, ], M[index_nk, selected_M, drop = FALSE])$union_vars
    basis_sel <- select_el_drm_basis(
      X = X[index_nk, x_vars, drop = FALSE],
      M = M[index_nk, selected_M, drop = FALSE],
      A = A[index_nk],
      stability_reps = 10,
      subsample_frac = 0.8,
      cv_folds = 5,
      alpha = 0.5,
      stability_cutoff = 0.8,
      seed = seed + 20000 * t + 200 * k
    )
    
    selected_basis <- basis_sel$selected_basis
    
    log_info(paste(
      "QR2 basis fold", k,
      "; selected_M =", paste(selected_M, collapse = ","),
      "; selected_basis =", paste(selected_basis, collapse = ","),
      "; successful_reps =", basis_sel$successful_reps,
      "; quadratic_freq =",
      paste(
        names(basis_sel$quadratic_frequency),
        round(basis_sel$quadratic_frequency, 2),
        sep = ":",
        collapse = ","
      )
    ))
    
    f_el_drm <- el_drm_fit_density_ratio(
      X[index_nk, x_vars, drop = FALSE],
      M[index_nk, selected_M, drop = FALSE],
      A[index_nk],
      basis_names = selected_basis,
      ridge_grid = sdr_ridge_grid,
      ridge_cv_folds = 5,
      ridge_cv_seed = seed + 10000 * t + 100 * k + 2
    )
    log_ratio_index <- el_drm_predict_log_ratio(f_el_drm,
                                                X[index, x_vars, drop = FALSE],
                                                M[index, selected_M, drop = FALSE])
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
      "theta_dim =", length(f_el_drm$theta)
    ))
    
    #mu
    lasso_mu<-lasso_var_selection(cbind(X[index_nk_A,], M[index_nk_A,]), Y[index_nk_A], task = 0)
    selected_vars_mu<-lasso_mu$selected_vars
    mu_rf=caret::train(x=cbind(X[index_nk_A,],M[index_nk_A,])[,selected_vars_mu,drop=FALSE], y=Y[index_nk_A], method= "ranger", tuneGrid=generate_tune_grid(length(index_nk_A),length(selected_vars_mu),0), trControl=ctrl_mu, importance="none", num.trees=generate_num_tree(length(index_nk_A)) )
    pred_mu[4,index]=predict(mu_rf$finalModel, data=cbind(X[index,],M[index,])[,selected_vars_mu,drop=FALSE])$predictions
    mtry_mu=c(mtry_mu,mu_rf$bestTune$mtry)
    minsize_mu=c(minsize_mu,mu_rf$bestTune$min.node.size)
    # summary(pred_mu[3,index]-mu[index])
    
    #tau1
    prediction_mu_1=predict(mu_rf$finalModel, data=cbind(X[index_nk_A0,],M[index_nk_A0,])[,selected_vars_mu,drop=FALSE])
    lasso_tau1<-lasso_var_selection(X[index_nk_A0,], prediction_mu_1$predictions, task = 0)
    selected_vars_tau1<-lasso_tau1$selected_vars
    tau1_rf=caret::train(x=X[index_nk_A0,selected_vars_tau1,drop=FALSE], y=prediction_mu_1$predictions, method= "ranger", tuneGrid=generate_tune_grid(length(index_nk_A0),length(selected_vars_tau1),0), trControl=ctrl_mu, importance="none", num.trees=generate_num_tree(length(index_nk_A0)) )
    pred_tau1[4,index]=predict(tau1_rf$finalModel, data=X[index,selected_vars_tau1,drop=FALSE])$predictions
    mtry_tau1=c(mtry_tau1,tau1_rf$bestTune$mtry)
    minsize_tau1=c(minsize_tau1,tau1_rf$bestTune$min.node.size)
    # summary(pred_tau1[3,index])
    
    #tau: use full outer-training-fold predictions for the pseudo-outcome.
    prediction_q_1 <- normalized_truncate_q(
      log_q = ref_log_q_sdr2_full,
      ref_log_q = ref_log_q_sdr2_full
    )
    
    prediction_mu_2_preds <- predict_caret_regression_final(
      mu_rf,
      newdata = cbind(X[index_nk_A,], M[index_nk_A,])[, selected_vars_mu, drop = FALSE]
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
    # beta_mu=as.vector(coefficients=coef(fit_lasso_mu$glmnet.fit, s=fit_lasso_mu$lambda.min))
    # summary(pred_mu[4,index]-mu[index])
  
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
      
    pred_p[4,index]=pmax(pmin(predict(p_rf_TR$finalModel, data=cbind(X[index,],M[index,])[,selected_vars_p,drop=FALSE])$predictions[, 2],0.99), 0.01)
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
    
    #f semiparametrically
    selected_M <- 1:ncol(M)
    ds=length(selected_M)
    x_vars <- multi_y_var_selection(X[index_nk, ], M[index_nk, selected_M, drop = FALSE])$union_vars
    selected_A1=x_vars#multi_y_var_selection(X[index_nk_A,], M[index_nk_A, selected_M])$union_vars
    selected_A0=x_vars#multi_y_var_selection(X[index_nk_A0,], M[index_nk_A0, selected_M])$union_vars
    
    # A=1: inner-CV tuning of mtry and nodesize
    rf_mean_A1 <- fit_mv_rf_cf(
      X_train = X[index_nk_A, selected_A1, drop = FALSE],
      M_train = M[index_nk_A, selected_M, drop = FALSE],
      K = 3,
      ntree_tune = 100,
      ntree_final = 50,
      seed = seed + 10000 * t + 100 * k + 1
    )
    
    mean_rf_A <- rf_mean_A1$fit
    mean_train_A1 <- rf_mean_A1$oob_pred
    residuals_A1 <- M[index_nk_A, selected_M, drop = FALSE] - mean_train_A1
    
    
    # A=0: independently tune mtry and nodesize
    rf_mean_A0 <- fit_mv_rf_cf(
      X_train = X[index_nk_A0, selected_A0, drop = FALSE],
      M_train = M[index_nk_A0, selected_M, drop = FALSE],
      K = 3,
      ntree_tune = 100,
      ntree_final = 50,
      seed = seed + 10000 * t + 100 * k + 2
    )
    
    mean_rf_A0 <- rf_mean_A0$fit
    mean_train_A0 <- rf_mean_A0$oob_pred
    residuals_A0 <- M[index_nk_A0, selected_M, drop = FALSE] - mean_train_A0
    
    log_info(paste(
      "TR mean RF fold", k,
      "A1: mtry =", rf_mean_A1$best_mtry,
      "nodesize =", rf_mean_A1$best_nodesize,
      "; A0: mtry =", rf_mean_A0$best_mtry,
      "nodesize =", rf_mean_A0$best_nodesize
    ))
    
    residuals_A1 <- sweep(residuals_A1, 2, colMeans(residuals_A1), "-")
    residuals_A0 <- sweep(residuals_A0, 2, colMeans(residuals_A0), "-")
    residuals_pool <- rbind(residuals_A1, residuals_A0)
    
    # cov_result <- compute_shrinkage_cov(residuals_pool, lambda = "auto")
    # sig_M_pool <- cov_result$matrix
    # colnames(sig_M_pool) <- colnames(M)[selected_M]
    # rownames(sig_M_pool) <- colnames(M)[selected_M]
    # sig_M_A1 <- sig_M_pool
    # sig_M_A0 <- sig_M_pool
    
    cov_A1 <- compute_shrinkage_cov(residuals_A1, lambda = "auto")
    cov_A0 <- compute_shrinkage_cov(residuals_A0, lambda = "auto")
    
    sig_M_A1 <- cov_A1$matrix
    sig_M_A0 <- cov_A0$matrix
    
    mean_pred_A0_on_A1 <- predict(mean_rf_A0,newdata = data.frame(X[index_nk_A, selected_A0, drop=FALSE]))
    mean_pred_A0_on_A1 <- sapply(1:length(selected_M),function(i) mean_pred_A0_on_A1$regrOutput[[i]][[1]])
    colnames(mean_pred_A0_on_A1) <- colnames(M)[selected_M]
    
    log_dens_A1_train <- dmvnorm(M[index_nk_A, selected_M, drop=FALSE] - mean_train_A1,sigma = sig_M_A1,log = TRUE)
    
    log_dens_A0_train <- dmvnorm(M[index_nk_A, selected_M, drop=FALSE] - mean_pred_A0_on_A1,sigma = sig_M_A0,log = TRUE)
    
    log_q_train <- log_dens_A0_train - log_dens_A1_train
    
    # stable log(mean(exp(log_q_train)))
    if (!all(is.finite(log_q_train))) {
      stop("Non-finite training density-ratio prediction.")
    }
    
    mean_pred <- predict(mean_rf_A,newdata=data.frame(X[index, selected_A1, drop=FALSE]))
    mean_pred_A=sapply(1:ds, function(i) {mean_pred$regrOutput[[i]][[1]]  })
    
    mean_pred <- predict(mean_rf_A0,newdata=data.frame(X[index, selected_A0, drop=FALSE]))
    mean_pred_A0=sapply(1:ds, function(i) {mean_pred$regrOutput[[i]][[1]]  })
    
    colnames(mean_pred_A0) <- colnames(M)[selected_M]
    colnames(mean_pred_A) <- colnames(M)[selected_M]
    
    log_dens_A1<-dmvnorm(M[index, selected_M]-mean_pred_A, sigma=sig_M_A1, log= TRUE)
    # log_dens_A1<-log(pmax(pmin(exp(log_dens_A1), 0.999), 0.001))
    log_dens_A0<-dmvnorm(M[index, selected_M]-mean_pred_A0, sigma=sig_M_A0, log= TRUE)
    # log_dens_A0<-log(pmax(pmin(exp(log_dens_A0), 0.999), 0.001))
    # pred_q[6,index]=pmax(pmin(exp(log_dens_A0 - log_dens_A1),100),0.01)
    
    log_q_test_raw <- log_dens_A0 - log_dens_A1
    
    pred_q[6,index] <- normalized_truncate_q(
      log_q = log_q_test_raw,
      ref_log_q = log_q_train
    )
    
    
    if (!all(is.finite(pred_q[6,index]))) {
      stop("Non-finite density ratio after exponentiation.")
    }
    # pred_q[6,index]=normalized_truncate_q(log_dens_A0, log_dens_A1,
    #                                       treated = A[index_nk] == 1)
    # pred_q[6,index]=exp(log_dens_A0 - log_dens_A1)
    
    # tau_1: reuse the same estimated F0(M|X)
    m_vars <- as.numeric(sub("^M", "", grep("^M", selected_vars_mu, value = TRUE)))
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
        rmvnorm(n_mc, mean=mean_pred_tau1[i, ], sigma=sig_M_tau1)
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
  set.seed(seed)

  #-------------------MSE check-------------
  mse_pi=rowMeans((pred_pi_a-matrix(rep(pi,each=num_pi_est),ncol=N))^2)
  mse_q=rowMeans((pred_q[,which(A==1)]-matrix(rep(q[which(A==1)],each=num_q_est),nrow=num_q_est))^2)
  mse_p=rowMeans((pred_p[,which(A==1)]-matrix(rep(p[which(A==1)],each=num_p_est),nrow=num_p_est))^2)
  mse_mu=rowMeans((pred_mu-matrix(rep(mu,each=num_mu_est),nrow=num_mu_est))^2)
  pred_tau[5:9,]=pred_tau1[5:9,]
  mse_tau=rowMeans((pred_tau-matrix(rep(tau,each=num_mu_est),nrow=num_mu_est))^2)    
  mse_tau1=rowMeans((pred_tau1-matrix(rep(tau,each=num_mu_est),nrow=num_mu_est))^2) 
  res_mu_tau=rowMeans((pred_q[,which(A==1)]*(Y[A==1]-pred_mu[,which(A==1)])+pred_tau1[,which(A==1)]-pred_tau[,which(A==1)]))
  
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
  theta_hat-theta;sigma_hat
  
  result=c(pred_theta=theta_hat,pred_sig2=sigma_hat,mse_tau,mse_mu,mse_q,mse_pi,mse_tau1,res_mu_tau,sel.rho_pi_a,sel.rho_pi_b,sel.rho_q,sel.rho_p,lambda_mu,lambda_tau1,lambda_tau,mtry_pi,mtry_p,mtry_mu,mtry_tau1,mtry_tau,minsize_pi,minsize_p,minsize_mu,minsize_tau1,minsize_tau)
  # result_list=list(result=result)                               
  temp_file <- tempfile(pattern = "temp_", fileext = ".rds")
  saveRDS(result, file = temp_file)
  file.rename(temp_file, file.path(res_dir, paste0("Res1_", current_time, "_", t, ".rds")))
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
order_new=c(1,2,5,8,3,4,6,7,9,10)
mse_pi=mse_pi[order_new];mse_q=mse_q[order_new];mse_mu=mse_mu[order_new];mse_tau=mse_tau[order_new];mse_tau1=mse_tau1[order_new];res_mu_tau=res_mu_tau[order_new]

# save(d1,d2,N,results_par, theta, pred_theta, pred_sig2, AL,AC, ESD,ASD, mse_pi, mse_q, mse_mu, mse_tau, mse_tau1,res_mu_tau, bias, file = file_name)

file_name <- paste0("Result_", current_time, "_", N, "_", d1, d2, "_", tmax, ".RData"); save(Phi,d1,d2,N,beta_pi_h,beta_m0_h,beta_mu_h,beta_mu_h_at_N_ref,N_ref_mu_nonlinear,mu_nonlinear_decay_rate,epsilon_sd,delta1_all,delta2_all,scale_XM,scale_Me,scale_mue,scale_xsq,scale_sub, scale_mix, func_link, Phi, Sigma, Sigma0, Sigma1, s_q, weight0, weight1, dgp_config, results_par, theta, pred_theta, pred_sig2, AL,AC, ESD,ASD, mse_pi, mse_q, mse_mu, mse_tau, mse_tau1,res_mu_tau, bias, file = file_name)
cat("bias=\n",paste(capture.output(bias), collapse = "\n"), file = log_file, append = TRUE)
cat("\n RMSE=\n",paste(capture.output(RMSE), collapse = "\n"), file = log_file, append = TRUE)
cat("\n AC=\n",paste(capture.output(AC), collapse = "\n"), file = log_file, append = TRUE)     
cat("\n AL=\n",paste(capture.output(AL), collapse = "\n"), file = log_file, append = TRUE) 
log_info(paste("\n", "Important data is saved in", file_name,"\n"))                              
cat(theta,"\n",bias,"\n",RMSE,"\n",AC,"\n",sig2,"\n",sum(A_all==1),"\n", sum(A_all==0),"\n")  

save(seed,num_totest,Phi,d1,d2,N,beta_pi_h,beta_m0_h,beta_mu_h,beta_mu_h_at_N_ref,N_ref_mu_nonlinear,mu_nonlinear_decay_rate,epsilon_sd,delta1_all,delta2_all,scale_XM,scale_Me,scale_mue,scale_xsq,scale_sub, scale_mix, func_link, Phi, Sigma, Sigma0, Sigma1, s_q, weight0, weight1, dgp_config, results_par, theta, pred_theta, pred_sig2, bias,RMSE, AL,AC, ESD,ASD, mse_pi, mse_q, mse_mu, mse_tau, mse_tau1,res_mu_tau, file = file_name) 
#  range_lambda_pi_a, range_lambda_pi_b, range_lambda_q, range_lambda_p, range_lambda_mu, range_lambda_tau1, range_lambda_tau,
cat("数据已保存到文件:", file_name, "\n")
cat(proc.time()-ptm)

# est=c("MQR1","MQR2","O-TR'","MTR2 (internal comparator)","QR1","QR2","TR","O-TR","O-TR2","Oracle")
# table_result=function(x){
#   paste(est[x],"&",format(round(bias[x],digits=3),nsmall=3,scientific=FALSE),"&",format(round(RMSE[x],digits=3),nsmall=3,scientific=FALSE),"&",
#         format(round(AL[x],digits=3),nsmall=3,scientific=FALSE),"&",format(round(AC[x],digits=2),nsmall=2,scientific=FALSE),"\\",sep="")
# }
# result=apply(t(1:num_totest),2,table_result)
# result=list(result=c("estimator&Bias&RMSE&AL&AC&ESD&ASD",result),est=est,
#             mse_prop_pi=round(mse_pi/as.numeric(var(pi_all)),digits=3),
#             mse_prop_q=round(mse_q/var(q_all[which(A_all==1)]),digits=3),
#             mse_prop_mu=round(mse_mu/var(mu_all[which(A_all==1)]),digits=3),
#             mse_prop_tau=round(mse_tau/var(tau_all[which(A_all==1)]),digits=3))
# # result
# proc.time()-ptm
