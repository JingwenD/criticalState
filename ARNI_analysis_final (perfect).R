#!/usr/bin/env Rscript
# ARNI网络推断分析 - 完整修复版
# 修复：
# 1. 使用原始基因表达数据（而非PCA）
# 2. 选择高变基因
# 3. 正确的ARNI公式：dx_i/dt ≈ Σ A_ij × x_j
# 4. Seurat v5兼容

library(Seurat)
library(pracma)
library(doParallel)
library(foreach)

# 设置并行计算
n_cores <- parallel::detectCores() - 4
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# ==============================================================================
# 数据准备函数
# ==============================================================================
prepare_data_for_arni <- function(seurat_obj, sample_id, n_genes = 100, 
                                   n_time_bins = 50, use_pca = FALSE) {
  
  # 筛选样本
  subset_idx <- seurat_obj$orig.ident == sample_id
  n_cells <- sum(subset_idx)
  
  if (n_cells == 0) {
    stop("没有找到匹配的细胞，样本ID: ", sample_id)
  }
  
  cat("  个体细胞数:", n_cells, "\n")
  
  # 获取标准化表达数据（Seurat v5兼容）
  if ("scale.data" %in% SeuratObject::Layers(seurat_obj@assays$Spatial)) {
    scale_data <- seurat_obj@assays$Spatial$scale.data
    cat("  使用scale.data层 (Seurat v5格式)\n")
  } else if ("scale.data" %in% names(seurat_obj@assays$Spatial@layers)) {
    scale_data <- seurat_obj@assays$Spatial@layers$scale.data
    cat("  使用scale.data层 (Seurat v4格式)\n")
  } else {
    stop("未找到scale.data层")
  }
  
  expr_data <- scale_data[, subset_idx, drop = FALSE]
  cat("  表达矩阵维度:", dim(expr_data), "\n")
  
  # 获取空间坐标
  coords_x <- seurat_obj@meta.data[subset_idx, "x", drop = TRUE]
  coords_y <- seurat_obj@meta.data[subset_idx, "y", drop = TRUE]
  
  # 选择基因
  if ("VariableFeatures" %in% names(seurat_obj@assays$Spatial)) {
    hvg <- seurat_obj@assays$Spatial@VariableFeatures
    cat("  高变基因数量:", length(hvg), "\n")
    n_select <- min(n_genes, length(hvg))
    hvg_top <- hvg[1:n_select]
    
    gene_names <- rownames(expr_data)
    gene_idx <- match(hvg_top, gene_names)
    gene_idx <- gene_idx[!is.na(gene_idx)]
    
    if (length(gene_idx) < 20) {
      gene_idx <- 1:min(100, nrow(expr_data))
    }
  } else {
    # 选择方差最大的基因
    gene_vars <- apply(expr_data, 1, var)
    gene_idx <- order(gene_vars, decreasing = TRUE)[1:n_genes]
  }
  
  expr_selected <- expr_data[gene_idx, ]
  cat("  选择基因数:", nrow(expr_selected), "\n")
  
  # 基于y坐标创建时间箱
  y_quantiles <- seq(0, 1, length.out = n_time_bins + 1)
  y_breaks <- quantile(coords_y, probs = y_quantiles)
  y_breaks[1] <- y_breaks[1] - 0.1
  y_breaks[length(y_breaks)] <- y_breaks[length(y_breaks)] + 0.1
  
  time_bins <- cut(coords_y, breaks = y_breaks, labels = FALSE, include.lowest = TRUE)
  unique_t_bins <- sort(unique(time_bins))
  n_timepoints <- length(unique_t_bins)
  
  cat("  时间箱数:", n_timepoints, "\n")
  
  # 构建时间序列：每个时间箱取平均
  time_series <- matrix(0, nrow = nrow(expr_selected), ncol = n_timepoints)
  
  for (t_idx in 1:n_timepoints) {
    t_bin <- unique_t_bins[t_idx]
    cell_idx <- which(time_bins == t_bin)
    time_series[, t_idx] <- rowMeans(expr_selected[, cell_idx, drop = FALSE], na.rm = TRUE)
  }
  
  cat("  最终数据维度:", dim(time_series), "\n")
  
  return(list(
    expr_data = time_series,
    gene_idx = gene_idx,
    n_timepoints = n_timepoints
  ))
}

# ==============================================================================
# ARNI核心算法
# ==============================================================================
arni_algorithm <- function(X, max_interactions = 20) {
  
  N <- nrow(X)
  M <- ncol(X)
  
  # 计算导数（前向差分）
  DX <- X[, 2:M] - X[, 1:(M-1)]
  X_reduced <- X[, 1:(M-1)]
  
  cat("  X维度:", dim(X), "\n")
  cat("  DX维度:", dim(DX), "\n")
  cat("  X_reduced维度:", dim(X_reduced), "\n")
  
  results <- list()
  
  # 对每个目标变量
  for (target_idx in 1:(N - 1)) {
    cat("  处理目标变量", target_idx, "/", N - 1, "\n")
    
    DX_target <- DX[target_idx, , drop = FALSE]
    remaining <- setdiff(1:N, target_idx)
    
    selected <- c()
    costs <- numeric(max_interactions)
    
    # 贪婪选择
    for (k in 1:max_interactions) {
      best_cost <- Inf
      best_j <- -1
      
      # 并行搜索
      costs_j <- foreach(j = remaining, .combine = c) %dopar% {
        tryCatch({
          current_idx <- c(selected, j)
          
          # 构建回归矩阵
          df <- data.frame(DX = as.vector(DX_target))
          for (idx in current_idx) {
            df[[paste0("X", idx)]] <- as.vector(X_reduced[idx, ])
          }
          
          formula <- as.formula(paste0("DX ~ ", paste0("+ X", current_idx, collapse = "")))
          lm_model <- lm(formula, data = df)
          
          DX_est <- predict(lm_model)
          DIFF <- as.vector(DX_target) - DX_est
          cost <- sqrt(mean(DIFF^2))
          
          return(cost)
        }, error = function(e) {
          return(Inf)
        })
      }
      
      # 找最佳
      best_pos <- which.min(costs_j)
      if (costs_j[best_pos] < Inf) {
        best_j <- remaining[best_pos]
        selected <- c(selected, best_j)
        remaining <- setdiff(remaining, best_j)
        costs[k] <- costs_j[best_pos]
      } else {
        costs[k] <- Inf
        break
      }
      
      # 提前停止
      if (k > 1 && costs[k] >= costs[k - 1]) {
        costs[k] <- Inf
        break
      }
    }
    
    # 计算学习曲线曲率找最优k
    optimal_k <- 1
    if (any(is.finite(costs))) {
      finite_costs <- costs[is.finite(costs)]
      if (length(finite_costs) > 2) {
        curvatures <- diff(diff(finite_costs))
        if (length(curvatures) > 0) {
          optimal_k <- which.min(curvatures) + 1
        }
      }
    }
    
    results[[target_idx]] <- list(
      interactions = selected[1:min(optimal_k, length(selected))],
      cost = costs,
      learning_curve = data.frame(n_interactions = 1:max_interactions, cost = costs),
      optimal_interactions = optimal_k
    )
  }
  
  return(results)
}

# ==============================================================================
# 构建网络矩阵
# ==============================================================================
build_network_matrix <- function(results, n_variables) {
  
  network_matrix <- matrix(0, nrow = n_variables, ncol = n_variables)
  
  for (target_idx in 1:(n_variables - 1)) {
    target_result <- results[[target_idx]]
    
    if (!is.null(target_result$interactions)) {
      for (source_idx in target_result$interactions) {
        network_matrix[target_idx, source_idx] <- 1
      }
    }
  }
  
  return(network_matrix)
}

# ==============================================================================
# 主分析函数
# ==============================================================================
run_arni_analysis <- function(seurat_obj, sample_id, 
                               n_genes = 100,
                               n_time_bins = 50,
                               max_interactions = 20,
                               output_dir = "ARNI_results") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("\n========================================\n")
  cat("ARNI网络推断分析\n")
  cat("样本:", sample_id, "\n")
  cat("参数: n_genes =", n_genes, ", n_time_bins =", n_time_bins, "\n")
  cat("========================================\n\n")
  
  # Step 1: 准备数据
  cat("准备数据...\n")
  ts_data <- prepare_data_for_arni(
    seurat_obj, sample_id,
    n_genes = n_genes,
    n_time_bins = n_time_bins
  )
  
  # Step 2: ARNI分析
  cat("\n执行ARNI算法...\n")
  results <-arni_algorithm(ts_data$expr_data, max_interactions = max_interactions)
  
  # Step 3: 构建网络
  cat("\n构建网络矩阵...\n")
  network_matrix <- build_network_matrix(results, nrow(ts_data$expr_data))
  
  # 保存结果
  saveRDS(results, file = paste0(output_dir, "/arni_results_", sample_id, ".rds"))
  saveRDS(network_matrix, file = paste0(output_dir, "/network_matrix_", sample_id, ".rds"))
  
  # 统计
  n_interactions <- sum(network_matrix > 0)
  cat("\n========================================\n")
  cat("分析完成!\n")
  cat("发现的网络交互数:", n_interactions, "\n")
  cat("========================================\n")
  
  # 输出每个节点的交互
  for (target_idx in 1:(nrow(ts_data$expr_data) - 1)) {
    target_result <- results[[target_idx]]
    if (!is.null(target_result$interactions) && length(target_result$interactions) > 0) {
      cat("  基因", target_idx, ":", target_result$interactions, "\n")
    }
  }
  
  return(list(
    results = results,
    network_matrix = network_matrix,
    sample_id = sample_id,
    ts_data = ts_data
  ))
}

# ==============================================================================
# 主程序
# ==============================================================================
main <- function() {
  
  cat("加载Seurat对象...\n")
  #skin1 <- readRDS("skin1.rds")
  
  # 获取样本ID
  sample_ids <- unique(skin1$orig.ident)
  cat("发现样本:", sample_ids, "\n")
  
  # 分析每个样本
  all_results <- list()
  
  for (sample_id in sample_ids) {
    cat("\n\n========================================\n")
    cat("分析个体:", sample_id, "\n")
    cat("========================================\n")
    
    tryCatch({
      result <- run_arni_analysis(
        skin1,
        sample_id = sample_id,
        n_genes = 100,
        n_time_bins = 50,
        max_interactions = 20,
        output_dir = "ARNI_results"
      )
      
      all_results[[sample_id]] <- result
      
    }, error = function(e) {
      cat("分析", sample_id, "时出错:", e$message, "\n")
    })
  }
  
  # 保存所有结果
  saveRDS(all_results, "ARNI_results/all_results.rds")
  
  cat("\n\n========================================\n")
  cat("所有分析完成!\n")
  cat("结果保存在: ARNI_results/\n")
  cat("========================================\n")
  
  return(all_results)
}

# 运行
results <- main()

# 清理
stopCluster(cl)
