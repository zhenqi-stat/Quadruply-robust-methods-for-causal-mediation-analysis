# 协变量扩充函数
expand_covariates <- function(data, X_conti, X_binar, remove_zero_var = TRUE) {
  
  # 输入验证
  if (!is.data.frame(data)) {
    stop("data must be a data.frame")
  }
  
  # 确保所有变量存在
  missing_conti <- setdiff(X_conti, names(data))
  missing_binar <- setdiff(X_binar, names(data))
  
  if (length(missing_conti) > 0 || length(missing_binar) > 0) {
    stop("Missing variables in data: ",
         ifelse(length(missing_conti) > 0, paste(missing_conti, collapse = ", "), ""),
         ifelse(length(missing_binar) > 0, paste(missing_binar, collapse = ", "), ""))
  }
  
  # 1. 连续变量的线性项
  result <- data[, X_conti, drop = FALSE]
  colnames(result) <- X_conti
  
  # 2. 连续变量的二次项
  conti_sq <- as.matrix(data[, X_conti, drop = FALSE]^2)
  colnames(conti_sq) <- paste0(X_conti, "_sq")
  result <- cbind(result, conti_sq)
  
  # 3. 二元变量的线性项
  binar_linear <- as.matrix(data[, X_binar, drop = FALSE])
  colnames(binar_linear) <- X_binar
  result <- cbind(result, binar_linear)
  
  # 4. 连续变量之间的两两交互项（线性×线性）
  n_conti <- length(X_conti)
  if (n_conti >= 2) {
    conti_inter_names <- character(0)
    conti_inter_matrix <- matrix(0, nrow = nrow(data), ncol = 0)
    
    for (i in 1:(n_conti - 1)) {
      for (j in (i + 1):n_conti) {
        inter <- data[[X_conti[i]]] * data[[X_conti[j]]]
        conti_inter_matrix <- cbind(conti_inter_matrix, inter)
        conti_inter_names <- c(conti_inter_names, 
                               paste0(X_conti[i], "_", X_conti[j]))
      }
    }
    
    colnames(conti_inter_matrix) <- conti_inter_names
    result <- cbind(result, conti_inter_matrix)
  }
  
  # 5. 二元变量之间的两两交互项
  n_binar <- length(X_binar)
  if (n_binar >= 2) {
    binar_inter_names <- character(0)
    binar_inter_matrix <- matrix(0, nrow = nrow(data), ncol = 0)
    
    for (i in 1:(n_binar - 1)) {
      for (j in (i + 1):n_binar) {
        inter <- data[[X_binar[i]]] * data[[X_binar[j]]]
        binar_inter_matrix <- cbind(binar_inter_matrix, inter)
        binar_inter_names <- c(binar_inter_names,
                               paste0(X_binar[i], "_", X_binar[j]))
      }
    }
    
    colnames(binar_inter_matrix) <- binar_inter_names
    result <- cbind(result, binar_inter_matrix)
  }
  
  # 6. 连续变量与二元变量的交互项（连续变量的线性项×二元变量）
  if (n_conti > 0 && n_binar > 0) {
    conti_binar_names <- character(0)
    conti_binar_matrix <- matrix(0, nrow = nrow(data), ncol = 0)
    
    for (c_var in X_conti) {
      for (b_var in X_binar) {
        inter <- data[[c_var]] * data[[b_var]]
        conti_binar_matrix <- cbind(conti_binar_matrix, inter)
        conti_binar_names <- c(conti_binar_names,
                               paste0(c_var, "_", b_var))
      }
    }
    
    colnames(conti_binar_matrix) <- conti_binar_names
    result <- cbind(result, conti_binar_matrix)
  }
  
  # 7. 连续变量的二次项与二元变量的交互项（连续变量的平方×二元变量）
  if (n_conti > 0 && n_binar > 0) {
    conti_sq_binar_names <- character(0)
    conti_sq_binar_matrix <- matrix(0, nrow = nrow(data), ncol = 0)
    
    for (c_var in X_conti) {
      c_sq_var <- paste0(c_var, "_sq")
      for (b_var in X_binar) {
        inter <- (data[[c_var]]^2) * data[[b_var]]
        conti_sq_binar_matrix <- cbind(conti_sq_binar_matrix, inter)
        conti_sq_binar_names <- c(conti_sq_binar_names,
                                  paste0(c_sq_var, "_", b_var))
      }
    }
    
    colnames(conti_sq_binar_matrix) <- conti_sq_binar_names
    result <- cbind(result, conti_sq_binar_matrix)
  }
  
  # 8. 可选：去除方差为0的列
  if (remove_zero_var) {
    variances <- apply(result, 2, var, na.rm = TRUE)
    zero_var_cols <- which(variances == 0)
    
    if (length(zero_var_cols) > 0) {
      cat("Removing", length(zero_var_cols), "columns with zero variance\n")
      result <- result[, -zero_var_cols, drop = FALSE]
    }
  }
  
  # 返回结果
  cat("Generated", ncol(result), "expanded covariates\n")
  return(result)
}


expand_covariates_ck <- function(data, X_conti, X_binar, remove_zero_var = TRUE) {
  
  # 输入验证
  if (!is.data.frame(data)) {
    stop("data must be a data.frame")
  }
  
  # 确保所有变量存在
  missing_conti <- setdiff(X_conti, names(data))
  missing_binar <- setdiff(X_binar, names(data))
  if (length(missing_conti) > 0 || length(missing_binar) > 0) {
    stop("Missing variables in data: ", 
         paste(c(missing_conti, missing_binar), collapse = ", "))
  }
  
  n <- nrow(data)
  
  # 1. 构建连续部分多项式基 (15项)
  # 包括: 常数项1, 一次项, 平方项, 两两交互项 (不包括常数项与变量的交互，因为那即是一次项)
  C_matrix <- matrix(1, nrow = n, ncol = 1)  # 常数项
  colnames(C_matrix) <- "const"
  
  # 一次项
  C_linear <- as.matrix(data[, X_conti, drop = FALSE])
  colnames(C_linear) <- X_conti
  C_matrix <- cbind(C_matrix, C_linear)
  
  # 平方项
  C_sq <- as.matrix(data[, X_conti, drop = FALSE]^2)
  colnames(C_sq) <- paste0(X_conti, "_sq")
  C_matrix <- cbind(C_matrix, C_sq)
  
  # 两两交互项 (连续变量之间)
  n_conti <- length(X_conti)
  if (n_conti >= 2) {
    C_inter_names <- character(0)
    C_inter_matrix <- matrix(0, nrow = n, ncol = 0)
    for (i in 1:(n_conti - 1)) {
      for (j in (i + 1):n_conti) {
        inter <- data[[X_conti[i]]] * data[[X_conti[j]]]
        C_inter_matrix <- cbind(C_inter_matrix, inter)
        C_inter_names <- c(C_inter_names, paste0(X_conti[i], "_", X_conti[j]))
      }
    }
    colnames(C_inter_matrix) <- C_inter_names
    C_matrix <- cbind(C_matrix, C_inter_matrix)
  }
  
  # 2. 构建二元部分多项式基 (29项)
  # 包括: 常数项1, 一次项, 两两交互项 (减去平方项，即二元变量自身的平方等于自身，故不单独列出)
  B_matrix <- matrix(1, nrow = n, ncol = 1)  # 常数项
  colnames(B_matrix) <- "const"
  
  # 一次项
  B_linear <- as.matrix(data[, X_binar, drop = FALSE])
  colnames(B_linear) <- X_binar
  B_matrix <- cbind(B_matrix, B_linear)
  
  # 两两交互项 (二元变量之间)
  n_binar <- length(X_binar)
  if (n_binar >= 2) {
    B_inter_names <- character(0)
    B_inter_matrix <- matrix(0, nrow = n, ncol = 0)
    for (i in 1:(n_binar - 1)) {
      for (j in (i + 1):n_binar) {
        inter <- data[[X_binar[i]]] * data[[X_binar[j]]]
        B_inter_matrix <- cbind(B_inter_matrix, inter)
        B_inter_names <- c(B_inter_names, paste0(X_binar[i], "_", X_binar[j]))
      }
    }
    colnames(B_inter_matrix) <- B_inter_names
    B_matrix <- cbind(B_matrix, B_inter_matrix)
  }
  
  # 3. 笛卡尔积: C的每一列与B的每一列相乘
  nC <- ncol(C_matrix)
  nB <- ncol(B_matrix)
  result_matrix <- matrix(0, nrow = n, ncol = nC * nB)
  result_names <- character(nC * nB)
  
  idx <- 1
  for (i in 1:nC) {
    for (j in 1:nB) {
      result_matrix[, idx] <- C_matrix[, i] * B_matrix[, j]
      result_names[idx] <- paste0(colnames(C_matrix)[i], "_x_", colnames(B_matrix)[j])
      idx <- idx + 1
    }
  }
  colnames(result_matrix) <- result_names
  
  # 4. 去掉常数项 (const_x_const) 因为通常模型会包含截距，且此项全为1
  const_idx <- which(result_names == "const_x_const")
  if (length(const_idx) > 0) {
    result_matrix <- result_matrix[, -const_idx, drop = FALSE]
    result_names <- result_names[-const_idx]
  }
  
  # 5. 可选：去除方差为0的列
  if (remove_zero_var) {
    variances <- apply(result_matrix, 2, var, na.rm = TRUE)
    zero_var_cols <- which(variances == 0)
    if (length(zero_var_cols) > 0) {
      cat("Removing", length(zero_var_cols), "columns with zero variance\n")
      result_matrix <- result_matrix[, -zero_var_cols, drop = FALSE]
    }
  }
  
  # 输出信息
  cat("Generated", ncol(result_matrix), "expanded covariates (from", 
      length(X_conti), "continuous and", length(X_binar), "binary variables)\n")
  
  return(result_matrix)
}

expand_covariates_200 <- function(data, X_conti, X_binar, remove_zero_var = TRUE) {
  
  # 输入验证
  if (!is.data.frame(data)) {
    stop("data must be a data.frame")
  }
  
  # 确保所有变量存在
  missing_conti <- setdiff(X_conti, names(data))
  missing_binar <- setdiff(X_binar, names(data))
  if (length(missing_conti) > 0 || length(missing_binar) > 0) {
    stop("Missing variables in data: ", 
         paste(c(missing_conti, missing_binar), collapse = ", "))
  }
  
  n <- nrow(data)
  
  # 1. 连续变量线性项 (p)
  result <- data[, X_conti, drop = FALSE]
  colnames(result) <- X_conti
  
  # 2. 连续变量二次项 (p)
  conti_sq <- as.matrix(data[, X_conti, drop = FALSE]^2)
  colnames(conti_sq) <- paste0(X_conti, "_sq")
  result <- cbind(result, conti_sq)
  
  # 3. 连续变量线性交互项 (C(p,2))
  p <- length(X_conti)
  if (p >= 2) {
    conti_inter_names <- character(0)
    conti_inter_matrix <- matrix(0, nrow = n, ncol = 0)
    for (i in 1:(p - 1)) {
      for (j in (i + 1):p) {
        inter <- data[[X_conti[i]]] * data[[X_conti[j]]]
        conti_inter_matrix <- cbind(conti_inter_matrix, inter)
        conti_inter_names <- c(conti_inter_names, 
                               paste0(X_conti[i], "_", X_conti[j]))
      }
    }
    colnames(conti_inter_matrix) <- conti_inter_names
    result <- cbind(result, conti_inter_matrix)
  }
  
  # 4. 二元变量线性项 (q)
  q <- length(X_binar)
  binar_linear <- as.matrix(data[, X_binar, drop = FALSE])
  colnames(binar_linear) <- X_binar
  result <- cbind(result, binar_linear)
  
  # 5. 二元变量交互项 (C(q,2))
  if (q >= 2) {
    binar_inter_names <- character(0)
    binar_inter_matrix <- matrix(0, nrow = n, ncol = 0)
    for (i in 1:(q - 1)) {
      for (j in (i + 1):q) {
        inter <- data[[X_binar[i]]] * data[[X_binar[j]]]
        binar_inter_matrix <- cbind(binar_inter_matrix, inter)
        binar_inter_names <- c(binar_inter_names,
                               paste0(X_binar[i], "_", X_binar[j]))
      }
    }
    colnames(binar_inter_matrix) <- binar_inter_names
    result <- cbind(result, binar_inter_matrix)
  }
  
  # 6. 连续变量线性项 × 二元变量 (p × q)
  conti_binar_names <- character(0)
  conti_binar_matrix <- matrix(0, nrow = n, ncol = 0)
  for (c_var in X_conti) {
    for (b_var in X_binar) {
      inter <- data[[c_var]] * data[[b_var]]
      conti_binar_matrix <- cbind(conti_binar_matrix, inter)
      conti_binar_names <- c(conti_binar_names, 
                             paste0(c_var, "_", b_var))
    }
  }
  colnames(conti_binar_matrix) <- conti_binar_names
  result <- cbind(result, conti_binar_matrix)
  
  # 7. 连续变量二次项 × 二元变量 (p × q)
  conti_sq_binar_names <- character(0)
  conti_sq_binar_matrix <- matrix(0, nrow = n, ncol = 0)
  for (c_var in X_conti) {
    c_sq <- data[[c_var]]^2
    for (b_var in X_binar) {
      inter <- c_sq * data[[b_var]]
      conti_sq_binar_matrix <- cbind(conti_sq_binar_matrix, inter)
      conti_sq_binar_names <- c(conti_sq_binar_names, 
                                paste0(c_var, "_sq_", b_var))
    }
  }
  colnames(conti_sq_binar_matrix) <- conti_sq_binar_names
  result <- cbind(result, conti_sq_binar_matrix)
  
  # 8. 连续变量线性项 × 二元变量交互项 (p × C(q,2))
  # 先获取二元变量交互项列名（已在第5步生成）
  if (q >= 2) {
    binar_inter_cols <- colnames(binar_inter_matrix)
    conti_binar_inter_names <- character(0)
    conti_binar_inter_matrix <- matrix(0, nrow = n, ncol = 0)
    for (c_var in X_conti) {
      for (b_inter in binar_inter_cols) {
        # 提取交互项中的两个二元变量
        b_vars <- strsplit(b_inter, "_")[[1]]
        inter_value <- data[[c_var]] * data[[b_vars[1]]] * data[[b_vars[2]]]
        conti_binar_inter_matrix <- cbind(conti_binar_inter_matrix, inter_value)
        conti_binar_inter_names <- c(conti_binar_inter_names, 
                                     paste0(c_var, "_", b_inter))
      }
    }
    colnames(conti_binar_inter_matrix) <- conti_binar_inter_names
    result <- cbind(result, conti_binar_inter_matrix)
  }
  
  # 9. 可选：去除方差为0的列
  if (remove_zero_var) {
    variances <- apply(result, 2, var, na.rm = TRUE)
    zero_var_cols <- which(variances == 0)
    if (length(zero_var_cols) > 0) {
      cat("Removing", length(zero_var_cols), "columns with zero variance\n")
      result <- result[, -zero_var_cols, drop = FALSE]
    }
  }
  
  # 输出信息
  cat("Generated", ncol(result), "expanded covariates (from", 
      length(X_conti), "continuous and", length(X_binar), "binary variables)\n")
  
  return(result)
}

# expand_covariates_with_mediation <- function(data, X_conti, X_binar, M_conti, 
#                                              add_conti_binar_interactions = FALSE,
#                                              add_M_X_interactions = TRUE,
#                                              remove_zero_var = TRUE) {
#   
#   # 输入验证
#   if (!is.data.frame(data)) {
#     stop("data must be a data.frame")
#   }
#   
#   # 确保所有变量存在
#   missing_vars <- setdiff(c(X_conti, X_binar, M_conti), names(data))
#   if (length(missing_vars) > 0) {
#     stop("Missing variables in data: ", paste(missing_vars, collapse = ", "))
#   }
#   
#   n <- nrow(data)
#   
#   # ========== 第一部分: X_expanded (仅使用X_conti和X_binar) ==========
#   
#   # 1.1 连续变量X_conti的线性项和二次项
#   X_linear <- as.matrix(data[, X_conti, drop = FALSE])
#   X_sq <- as.matrix(data[, X_conti, drop = FALSE]^2)
#   
#   # 1.2 连续变量X_conti之间的两两交互项(线性×线性)
#   X_inter_names <- character(0)
#   X_inter_matrix <- matrix(0, nrow = n, ncol = 0)
#   n_X_conti <- length(X_conti)
#   
#   if (n_X_conti >= 2) {
#     for (i in 1:(n_X_conti - 1)) {
#       for (j in (i + 1):n_X_conti) {
#         inter <- X_linear[, i] * X_linear[, j]
#         X_inter_matrix <- cbind(X_inter_matrix, inter)
#         X_inter_names <- c(X_inter_names, paste0(X_conti[i], "_", X_conti[j]))
#       }
#     }
#   }
#   
#   # 1.3 二元变量X_binar的线性项
#   B_linear <- as.matrix(data[, X_binar, drop = FALSE])
#   
#   # 1.4 二元变量X_binar之间的两两交互项
#   B_inter_names <- character(0)
#   B_inter_matrix <- matrix(0, nrow = n, ncol = 0)
#   n_X_binar <- length(X_binar)
#   
#   if (n_X_binar >= 2) {
#     for (i in 1:(n_X_binar - 1)) {
#       for (j in (i + 1):n_X_binar) {
#         inter <- B_linear[, i] * B_linear[, j]
#         B_inter_matrix <- cbind(B_inter_matrix, inter)
#         B_inter_names <- c(B_inter_names, paste0(X_binar[i], "_", X_binar[j]))
#       }
#     }
#   }
#   
#   # 初始化X_expanded
#   X_expanded <- cbind(
#     X_linear, X_sq, X_inter_matrix, B_linear, B_inter_matrix
#   )
#   
#   colnames(X_expanded) <- c(
#     X_conti, 
#     paste0(X_conti, "_sq"),
#     X_inter_names,
#     X_binar,
#     B_inter_names
#   )
#   
#   # 1.5 可选: 连续变量(X_conti)与二元变量的交互项
#   if (add_conti_binar_interactions) {
#     # 连续变量的线性项×二元变量
#     XB_lin_inter_names <- character(0)
#     XB_lin_inter_matrix <- matrix(0, nrow = n, ncol = 0)
#     
#     # 连续变量的二次项×二元变量
#     XB_sq_inter_names <- character(0)
#     XB_sq_inter_matrix <- matrix(0, nrow = n, ncol = 0)
#     
#     for (c_var in X_conti) {
#       for (b_var in X_binar) {
#         # 线性项交互
#         inter_lin <- data[[c_var]] * data[[b_var]]
#         XB_lin_inter_matrix <- cbind(XB_lin_inter_matrix, inter_lin)
#         XB_lin_inter_names <- c(XB_lin_inter_names, paste0(c_var, "_", b_var))
#         
#         # 二次项交互
#         inter_sq <- (data[[c_var]]^2) * data[[b_var]]
#         XB_sq_inter_matrix <- cbind(XB_sq_inter_matrix, inter_sq)
#         XB_sq_inter_names <- c(XB_sq_inter_names, paste0(c_var, "_sq_", b_var))
#       }
#     }
#     
#     # 添加到X_expanded
#     X_expanded <- cbind(
#       X_expanded, 
#       XB_lin_inter_matrix, 
#       XB_sq_inter_matrix
#     )
#     
#     colnames(X_expanded) <- c(
#       colnames(X_expanded)[1:(ncol(X_expanded) - ncol(XB_lin_inter_matrix) - ncol(XB_sq_inter_matrix))],
#       XB_lin_inter_names,
#       XB_sq_inter_names
#     )
#   }
#   
#   # ========== 第二部分: M_expanded (使用M_conti和X_conti) ==========
#   
#   # 2.1 中介变量M_conti的线性项和二次项
#   M_linear <- as.matrix(data[, M_conti, drop = FALSE])
#   M_sq <- as.matrix(data[, M_conti, drop = FALSE]^2)
#   
#   # 2.2 中介变量M_conti之间的两两交互项(线性×线性)
#   M_inter_names <- character(0)
#   M_inter_matrix <- matrix(0, nrow = n, ncol = 0)
#   n_M_conti <- length(M_conti)
#   
#   if (n_M_conti >= 2) {
#     for (i in 1:(n_M_conti - 1)) {
#       for (j in (i + 1):n_M_conti) {
#         inter <- M_linear[, i] * M_linear[, j]
#         M_inter_matrix <- cbind(M_inter_matrix, inter)
#         M_inter_names <- c(M_inter_names, paste0(M_conti[i], "_", M_conti[j]))
#       }
#     }
#   }
#   
#   # 初始化M_expanded
#   M_expanded <- cbind(M_linear, M_sq, M_inter_matrix)
#   colnames(M_expanded) <- c(
#     M_conti,
#     paste0(M_conti, "_sq"),
#     M_inter_names
#   )
#   
#   # 2.3 M的线性项与X的连续变量交互
#   if (add_M_X_interactions) {
#     # M的线性项与X_conti的线性项交互
#     M_X_lin_inter_names <- character(0)
#     M_X_lin_inter_matrix <- matrix(0, nrow = n, ncol = 0)
#     
#     # M的线性项与X_conti的二次项交互
#     M_X_sq_inter_names <- character(0)
#     M_X_sq_inter_matrix <- matrix(0, nrow = n, ncol = 0)
#     
#     # M的线性项与X_binar交互
#     M_B_inter_names <- character(0)
#     M_B_inter_matrix <- matrix(0, nrow = n, ncol = 0)
#     
#     for (m_var in M_conti) {
#       # 与X_conti线性项交互
#       for (c_var in X_conti) {
#         inter_lin <- data[[m_var]] * data[[c_var]]
#         M_X_lin_inter_matrix <- cbind(M_X_lin_inter_matrix, inter_lin)
#         M_X_lin_inter_names <- c(M_X_lin_inter_names, paste0(m_var, "_", c_var))
#       }
#       
#       # 与X_conti二次项交互
#       for (c_var in X_conti) {
#         inter_sq <- data[[m_var]] * (data[[c_var]]^2)
#         M_X_sq_inter_matrix <- cbind(M_X_sq_inter_matrix, inter_sq)
#         M_X_sq_inter_names <- c(M_X_sq_inter_names, paste0(m_var, "_", c_var, "_sq"))
#       }
#       
#       # 与X_binar交互
#       for (b_var in X_binar) {
#         inter_bin <- data[[m_var]] * data[[b_var]]
#         M_B_inter_matrix <- cbind(M_B_inter_matrix, inter_bin)
#         M_B_inter_names <- c(M_B_inter_names, paste0(m_var, "_", b_var))
#       }
#     }
#     
#     # 添加到M_expanded
#     M_expanded <- cbind(
#       M_expanded, 
#       M_X_lin_inter_matrix, 
#       M_X_sq_inter_matrix,
#       M_B_inter_matrix
#     )
#     
#     colnames(M_expanded) <- c(
#       colnames(M_expanded)[1:(ncol(M_expanded) - ncol(M_X_lin_inter_matrix) - ncol(M_X_sq_inter_matrix) - ncol(M_B_inter_matrix))],
#       M_X_lin_inter_names,
#       M_X_sq_inter_names,
#       M_B_inter_names
#     )
#   }
#   
#   # 2.4 可选: 如果需要，M的二次项与X的交互（默认不添加）
#   # 这部分可以根据需要添加，但根据用户要求，目前不添加
#   
#   # ========== 第三部分: 去除方差为0的列 ==========
#   
#   if (remove_zero_var) {
#     # 处理X_expanded
#     X_variances <- apply(X_expanded, 2, var, na.rm = TRUE)
#     X_zero_var_cols <- which(X_variances == 0)
#     if (length(X_zero_var_cols) > 0) {
#       cat("Removing", length(X_zero_var_cols), "columns with zero variance from X_expanded\n")
#       X_expanded <- X_expanded[, -X_zero_var_cols, drop = FALSE]
#     }
#     
#     # 处理M_expanded
#     M_variances <- apply(M_expanded, 2, var, na.rm = TRUE)
#     M_zero_var_cols <- which(M_variances == 0)
#     if (length(M_zero_var_cols) > 0) {
#       cat("Removing", length(M_zero_var_cols), "columns with zero variance from M_expanded\n")
#       M_expanded <- M_expanded[, -M_zero_var_cols, drop = FALSE]
#     }
#   }
#   
#   # ========== 第四部分: 输出 ==========
#   
#   cat("\n=== 变量生成统计 ===\n")
#   cat("X_expanded: ", ncol(X_expanded), "个变量\n")
#   cat("M_expanded: ", ncol(M_expanded), "个变量\n")
#   cat("总计: ", ncol(X_expanded) + ncol(M_expanded), "个变量\n")
#   
#   # 计算变量数量
#   n_X_linear <- length(X_conti)
#   n_X_sq <- length(X_conti)
#   n_X_inter <- ifelse(n_X_conti >= 2, choose(n_X_conti, 2), 0)
#   n_B_linear <- length(X_binar)
#   n_B_inter <- ifelse(n_X_binar >= 2, choose(n_X_binar, 2), 0)
#   
#   n_XB_lin_inter <- ifelse(add_conti_binar_interactions, n_X_conti * n_X_binar, 0)
#   n_XB_sq_inter <- ifelse(add_conti_binar_interactions, n_X_conti * n_X_binar, 0)
#   
#   n_M_linear <- length(M_conti)
#   n_M_sq <- length(M_conti)
#   n_M_inter <- ifelse(n_M_conti >= 2, choose(n_M_conti, 2), 0)
#   
#   n_M_X_lin_inter <- ifelse(add_M_X_interactions, n_M_conti * n_X_conti, 0)
#   n_M_X_sq_inter <- ifelse(add_M_X_interactions, n_M_conti * n_X_conti, 0)
#   n_M_B_inter <- ifelse(add_M_X_interactions, n_M_conti * n_X_binar, 0)
#   
#   cat("\n=== 变量构成分解 ===\n")
#   cat("X_expanded 构成:\n")
#   cat("  - X_conti 线性项: ", n_X_linear, "\n")
#   cat("  - X_conti 二次项: ", n_X_sq, "\n")
#   cat("  - X_conti 两两交互(线性×线性): ", n_X_inter, "\n")
#   cat("  - X_binar 线性项: ", n_B_linear, "\n")
#   cat("  - X_binar 两两交互: ", n_B_inter, "\n")
#   if (add_conti_binar_interactions) {
#     cat("  - X_conti线性项 × X_binar: ", n_XB_lin_inter, "\n")
#     cat("  - X_conti二次项 × X_binar: ", n_XB_sq_inter, "\n")
#   }
#   
#   cat("\nM_expanded 构成:\n")
#   cat("  - M_conti 线性项: ", n_M_linear, "\n")
#   cat("  - M_conti 二次项: ", n_M_sq, "\n")
#   cat("  - M_conti 两两交互(线性×线性): ", n_M_inter, "\n")
#   if (add_M_X_interactions) {
#     cat("  - M_conti线性项 × X_conti线性项: ", n_M_X_lin_inter, "\n")
#     cat("  - M_conti线性项 × X_conti二次项: ", n_M_X_sq_inter, "\n")
#     cat("  - M_conti线性项 × X_binar: ", n_M_B_inter, "\n")
#   }
#   
#   # 返回结果
#   return(list(
#     X_expanded = X_expanded,
#     M_expanded = M_expanded,
#     all_expanded = cbind(X_expanded, M_expanded)
#   ))
# }

expand_covariates_with_mediation <- function(data, X_conti, X_binar, M_conti,
                                             add_conti_binar_interactions = FALSE,
                                             add_M_X_interactions = TRUE,
                                             add_M_X_interactions_extended = TRUE,  # 新增参数控制是否添加扩展的M-X交互
                                             remove_zero_var = TRUE) {

  # 输入验证
  if (!is.data.frame(data)) {
    stop("data must be a data.frame")
  }

  # 确保所有变量存在
  missing_vars <- setdiff(c(X_conti, X_binar, M_conti), names(data))
  if (length(missing_vars) > 0) {
    stop("Missing variables in data: ", paste(missing_vars, collapse = ", "))
  }

  n <- nrow(data)

  # ========== 第一部分: X_expanded (仅使用X_conti和X_binar) ==========

  # 1.1 连续变量X_conti的线性项和二次项
  X_linear <- as.matrix(data[, X_conti, drop = FALSE])
  X_sq <- as.matrix(data[, X_conti, drop = FALSE]^2)

  # 1.2 连续变量X_conti之间的两两交互项(线性×线性)
  X_inter_names <- character(0)
  X_inter_matrix <- matrix(0, nrow = n, ncol = 0)
  n_X_conti <- length(X_conti)

  if (n_X_conti >= 2) {
    for (i in 1:(n_X_conti - 1)) {
      for (j in (i + 1):n_X_conti) {
        inter <- X_linear[, i] * X_linear[, j]
        X_inter_matrix <- cbind(X_inter_matrix, inter)
        X_inter_names <- c(X_inter_names, paste0(X_conti[i], "_", X_conti[j]))
      }
    }
  }

  # 1.3 二元变量X_binar的线性项
  B_linear <- as.matrix(data[, X_binar, drop = FALSE])

  # 1.4 二元变量X_binar之间的两两交互项
  B_inter_names <- character(0)
  B_inter_matrix <- matrix(0, nrow = n, ncol = 0)
  n_X_binar <- length(X_binar)

  if (n_X_binar >= 2) {
    for (i in 1:(n_X_binar - 1)) {
      for (j in (i + 1):n_X_binar) {
        inter <- B_linear[, i] * B_linear[, j]
        B_inter_matrix <- cbind(B_inter_matrix, inter)
        B_inter_names <- c(B_inter_names, paste0(X_binar[i], "_", X_binar[j]))
      }
    }
  }

  # 初始化X_expanded
  X_expanded <- cbind(
    X_linear, X_sq, X_inter_matrix, B_linear, B_inter_matrix
  )

  colnames(X_expanded) <- c(
    X_conti,
    paste0(X_conti, "_sq"),
    X_inter_names,
    X_binar,
    B_inter_names
  )

  # 1.5 可选: 连续变量(X_conti)与二元变量的交互项
  if (add_conti_binar_interactions) {
    # 连续变量的线性项×二元变量
    XB_lin_inter_names <- character(0)
    XB_lin_inter_matrix <- matrix(0, nrow = n, ncol = 0)

    # 连续变量的二次项×二元变量
    XB_sq_inter_names <- character(0)
    XB_sq_inter_matrix <- matrix(0, nrow = n, ncol = 0)

    # 新增: X_conti两两交互项×二元变量
    XB_conti_inter_inter_names <- character(0)
    XB_conti_inter_inter_matrix <- matrix(0, nrow = n, ncol = 0)

    for (c_var in X_conti) {
      for (b_var in X_binar) {
        # 线性项交互
        inter_lin <- data[[c_var]] * data[[b_var]]
        XB_lin_inter_matrix <- cbind(XB_lin_inter_matrix, inter_lin)
        XB_lin_inter_names <- c(XB_lin_inter_names, paste0(c_var, "_", b_var))

        # 二次项交互
        inter_sq <- (data[[c_var]]^2) * data[[b_var]]
        XB_sq_inter_matrix <- cbind(XB_sq_inter_matrix, inter_sq)
        XB_sq_inter_names <- c(XB_sq_inter_names, paste0(c_var, "_sq_", b_var))
      }
    }

    # 添加: X_conti两两交互×X_binar
    if (n_X_conti >= 2) {
      for (i in 1:(n_X_conti - 1)) {
        for (j in (i + 1):n_X_conti) {
          x_inter <- X_linear[, i] * X_linear[, j]
          inter_name <- paste0(X_conti[i], "_", X_conti[j])

          for (b_var in X_binar) {
            inter_inter <- x_inter * data[[b_var]]
            XB_conti_inter_inter_matrix <- cbind(XB_conti_inter_inter_matrix, inter_inter)
            XB_conti_inter_inter_names <- c(XB_conti_inter_inter_names,
                                            paste0(inter_name, "_", b_var))
          }
        }
      }
    }

    # 添加到X_expanded
    X_expanded <- cbind(
      X_expanded,
      XB_lin_inter_matrix,
      XB_sq_inter_matrix,
      XB_conti_inter_inter_matrix
    )

    colnames(X_expanded) <- c(
      colnames(X_expanded)[1:(ncol(X_expanded) - ncol(XB_lin_inter_matrix) -
                                ncol(XB_sq_inter_matrix) -
                                ncol(XB_conti_inter_inter_matrix))],
      XB_lin_inter_names,
      XB_sq_inter_names,
      XB_conti_inter_inter_names
    )
  }

  # ========== 第二部分: M_expanded (使用M_conti和X_conti) ==========

  # 2.1 中介变量M_conti的线性项和二次项
  M_linear <- as.matrix(data[, M_conti, drop = FALSE])
  M_sq <- as.matrix(data[, M_conti, drop = FALSE]^2)

  # 2.2 中介变量M_conti之间的两两交互项(线性×线性)
  M_inter_names <- character(0)
  M_inter_matrix <- matrix(0, nrow = n, ncol = 0)
  n_M_conti <- length(M_conti)

  if (n_M_conti >= 2) {
    for (i in 1:(n_M_conti - 1)) {
      for (j in (i + 1):n_M_conti) {
        inter <- M_linear[, i] * M_linear[, j]
        M_inter_matrix <- cbind(M_inter_matrix, inter)
        M_inter_names <- c(M_inter_names, paste0(M_conti[i], "_", M_conti[j]))
      }
    }
  }

  # 初始化M_expanded
  M_expanded <- cbind(M_linear, M_sq, M_inter_matrix)
  colnames(M_expanded) <- c(
    M_conti,
    paste0(M_conti, "_sq"),
    M_inter_names
  )

  # 2.3 M的线性项与X的交互
  if (add_M_X_interactions) {
    # M的线性项与X_conti的线性项交互
    M_X_lin_inter_names <- character(0)
    M_X_lin_inter_matrix <- matrix(0, nrow = n, ncol = 0)

    # M的线性项与X_conti的二次项交互
    M_X_sq_inter_names <- character(0)
    M_X_sq_inter_matrix <- matrix(0, nrow = n, ncol = 0)

    # M的线性项与X_binar交互
    M_B_inter_names <- character(0)
    M_B_inter_matrix <- matrix(0, nrow = n, ncol = 0)

    # 新增: M的线性项与X_conti两两交互的交互
    M_X_conti_inter_inter_names <- character(0)
    M_X_conti_inter_inter_matrix <- matrix(0, nrow = n, ncol = 0)

    # 新增: M的线性项与X_binar两两交互的交互
    M_X_binar_inter_inter_names <- character(0)
    M_X_binar_inter_inter_matrix <- matrix(0, nrow = n, ncol = 0)

    for (m_var in M_conti) {
      # 与X_conti线性项交互
      for (c_var in X_conti) {
        inter_lin <- data[[m_var]] * data[[c_var]]
        M_X_lin_inter_matrix <- cbind(M_X_lin_inter_matrix, inter_lin)
        M_X_lin_inter_names <- c(M_X_lin_inter_names, paste0(m_var, "_", c_var))
      }

      # 与X_conti二次项交互
      for (c_var in X_conti) {
        inter_sq <- data[[m_var]] * (data[[c_var]]^2)
        M_X_sq_inter_matrix <- cbind(M_X_sq_inter_matrix, inter_sq)
        M_X_sq_inter_names <- c(M_X_sq_inter_names, paste0(m_var, "_", c_var, "_sq"))
      }

      # 与X_binar交互
      for (b_var in X_binar) {
        inter_bin <- data[[m_var]] * data[[b_var]]
        M_B_inter_matrix <- cbind(M_B_inter_matrix, inter_bin)
        M_B_inter_names <- c(M_B_inter_names, paste0(m_var, "_", b_var))
      }
    }

    # 新增: M的线性项与X_conti两两交互的交互
    if (add_M_X_interactions_extended && n_X_conti >= 2) {
      for (m_var in M_conti) {
        for (i in 1:(n_X_conti - 1)) {
          for (j in (i + 1):n_X_conti) {
            inter <- data[[m_var]] * X_linear[, i] * X_linear[, j]
            M_X_conti_inter_inter_matrix <- cbind(M_X_conti_inter_inter_matrix, inter)
            M_X_conti_inter_inter_names <- c(M_X_conti_inter_inter_names,
                                             paste0(m_var, "_", X_conti[i], "_", X_conti[j]))
          }
        }
      }
    }

    # 新增: M的线性项与X_binar两两交互的交互
    if (add_M_X_interactions_extended && n_X_binar >= 2) {
      for (m_var in M_conti) {
        for (i in 1:(n_X_binar - 1)) {
          for (j in (i + 1):n_X_binar) {
            inter <- data[[m_var]] * B_linear[, i] * B_linear[, j]
            M_X_binar_inter_inter_matrix <- cbind(M_X_binar_inter_inter_matrix, inter)
            M_X_binar_inter_inter_names <- c(M_X_binar_inter_inter_names,
                                             paste0(m_var, "_", X_binar[i], "_", X_binar[j]))
          }
        }
      }
    }

    # 添加到M_expanded
    M_expanded <- cbind(
      M_expanded,
      M_X_lin_inter_matrix,
      M_X_sq_inter_matrix,
      M_B_inter_matrix
    )

    # 如果有扩展的交互项，也添加进去
    if (add_M_X_interactions_extended) {
      if (ncol(M_X_conti_inter_inter_matrix) > 0) {
        M_expanded <- cbind(M_expanded, M_X_conti_inter_inter_matrix)
      }
      if (ncol(M_X_binar_inter_inter_matrix) > 0) {
        M_expanded <- cbind(M_expanded, M_X_binar_inter_inter_matrix)
      }
    }

    # 更新列名
    current_colnames <- colnames(M_expanded)[1:(ncol(M_expanded) -
                                                  ncol(M_X_lin_inter_matrix) -
                                                  ncol(M_X_sq_inter_matrix) -
                                                  ncol(M_B_inter_matrix) -
                                                  ifelse(add_M_X_interactions_extended,
                                                         ncol(M_X_conti_inter_inter_matrix) +
                                                           ncol(M_X_binar_inter_inter_matrix), 0))]

    colnames(M_expanded) <- c(
      current_colnames,
      M_X_lin_inter_names,
      M_X_sq_inter_names,
      M_B_inter_names,
      if (add_M_X_interactions_extended) M_X_conti_inter_inter_names else NULL,
      if (add_M_X_interactions_extended) M_X_binar_inter_inter_names else NULL
    )
  }

  # ========== 第三部分: 去除方差为0的列 ==========

  if (remove_zero_var) {
    # 处理X_expanded
    X_variances <- apply(X_expanded, 2, var, na.rm = TRUE)
    X_zero_var_cols <- which(X_variances == 0)
    if (length(X_zero_var_cols) > 0) {
      cat("Removing", length(X_zero_var_cols), "columns with zero variance from X_expanded\n")
      X_expanded <- X_expanded[, -X_zero_var_cols, drop = FALSE]
    }

    # 处理M_expanded
    M_variances <- apply(M_expanded, 2, var, na.rm = TRUE)
    M_zero_var_cols <- which(M_variances == 0)
    if (length(M_zero_var_cols) > 0) {
      cat("Removing", length(M_zero_var_cols), "columns with zero variance from M_expanded\n")
      M_expanded <- M_expanded[, -M_zero_var_cols, drop = FALSE]
    }
  }

  # ========== 第四部分: 输出 ==========

  cat("\n=== 变量生成统计 ===\n")
  cat("X_expanded: ", ncol(X_expanded), "个变量\n")
  cat("M_expanded: ", ncol(M_expanded), "个变量\n")
  cat("总计: ", ncol(X_expanded) + ncol(M_expanded), "个变量\n")

  # 计算变量数量
  n_X_linear <- length(X_conti)
  n_X_sq <- length(X_conti)
  n_X_inter <- ifelse(n_X_conti >= 2, choose(n_X_conti, 2), 0)
  n_B_linear <- length(X_binar)
  n_B_inter <- ifelse(n_X_binar >= 2, choose(n_X_binar, 2), 0)

  n_XB_lin_inter <- ifelse(add_conti_binar_interactions, n_X_conti * n_X_binar, 0)
  n_XB_sq_inter <- ifelse(add_conti_binar_interactions, n_X_conti * n_X_binar, 0)
  n_XB_conti_inter_inter <- ifelse(add_conti_binar_interactions && n_X_conti >= 2,
                                   choose(n_X_conti, 2) * n_X_binar, 0)

  n_M_linear <- length(M_conti)
  n_M_sq <- length(M_conti)
  n_M_inter <- ifelse(n_M_conti >= 2, choose(n_M_conti, 2), 0)

  n_M_X_lin_inter <- ifelse(add_M_X_interactions, n_M_conti * n_X_conti, 0)
  n_M_X_sq_inter <- ifelse(add_M_X_interactions, n_M_conti * n_X_conti, 0)
  n_M_B_inter <- ifelse(add_M_X_interactions, n_M_conti * n_X_binar, 0)

  n_M_X_conti_inter_inter <- ifelse(add_M_X_interactions && add_M_X_interactions_extended && n_X_conti >= 2,
                                    n_M_conti * choose(n_X_conti, 2), 0)
  n_M_X_binar_inter_inter <- ifelse(add_M_X_interactions && add_M_X_interactions_extended && n_X_binar >= 2,
                                    n_M_conti * choose(n_X_binar, 2), 0)

  cat("\n=== 变量构成分解 ===\n")
  cat("X_expanded 构成:\n")
  cat("  - X_conti 线性项: ", n_X_linear, "\n")
  cat("  - X_conti 二次项: ", n_X_sq, "\n")
  cat("  - X_conti 两两交互(线性×线性): ", n_X_inter, "\n")
  cat("  - X_binar 线性项: ", n_B_linear, "\n")
  cat("  - X_binar 两两交互: ", n_B_inter, "\n")
  if (add_conti_binar_interactions) {
    cat("  - X_conti线性项 × X_binar: ", n_XB_lin_inter, "\n")
    cat("  - X_conti二次项 × X_binar: ", n_XB_sq_inter, "\n")
    cat("  - X_conti两两交互 × X_binar: ", n_XB_conti_inter_inter, "\n")
  }

  cat("\nM_expanded 构成:\n")
  cat("  - M_conti 线性项: ", n_M_linear, "\n")
  cat("  - M_conti 二次项: ", n_M_sq, "\n")
  cat("  - M_conti 两两交互(线性×线性): ", n_M_inter, "\n")
  if (add_M_X_interactions) {
    cat("  - M_conti线性项 × X_conti线性项: ", n_M_X_lin_inter, "\n")
    cat("  - M_conti线性项 × X_conti二次项: ", n_M_X_sq_inter, "\n")
    cat("  - M_conti线性项 × X_binar: ", n_M_B_inter, "\n")
    if (add_M_X_interactions_extended) {
      cat("  - M_conti线性项 × X_conti两两交互: ", n_M_X_conti_inter_inter, "\n")
      cat("  - M_conti线性项 × X_binar两两交互: ", n_M_X_binar_inter_inter, "\n")
    }
  }

  # 返回结果
  return(list(
    X_expanded = X_expanded,
    M_expanded = M_expanded,
    all_expanded = cbind(X_expanded, M_expanded)
  ))
}
