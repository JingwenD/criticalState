#!/usr/bin/env Rscript
# ARNI测试 - 使用原始基因表达数据
# 替代方案：如果PCA效果不好，使用高变基因的表达数据

library(Seurat)
library(pracma)

# 加载数据
#skin1 <- readRDS("skin1.rds")
sample_id <- "lsm"
subset_idx <- skin1$orig.ident == sample_id

cat("=== 使用原始基因表达数据测试ARNI ===\n")
cat("样本:", sample_id, "\n")
cat("筛选前细胞数:", nrow(skin1@meta.data), "\n")
cat("筛选后细胞数:", sum(subset_idx), "\n")

# 获取标准化表达数据（Seurat v5兼容）
if ("scale.data" %in% SeuratObject::Layers(skin1@assays$Spatial)) {
  scale_data <- skin1@assays$Spatial$scale.data
  cat("使用scale.data层\n")
} else if ("scale.data" %in% names(skin1@assays$Spatial@layers)) {
  scale_data <- skin1@assays$Spatial@layers$scale.data
  cat("使用scale.data层 (v4格式)\n")
} else {
  cat("未找到scale.data，使用data层并标准化\n")
  if ("data" %in% SeuratObject::Layers(skin1@assays$Spatial)) {
    data_mat <- skin1@assays$Spatial$data
  } else {
    data_mat <- skin1@assays$Spatial@layers$data
  }
  scale_data <- scale(t(as.matrix(data_mat)))
}

# 筛选样本
expr_data <- scale_data[, subset_idx, drop = FALSE]
cat("表达矩阵维度:", dim(expr_data), "\n")

# 选择高变基因（如果已有）
if ("VariableFeatures" %in% names(skin1@assays$Spatial)) {
  hvg <- skin1@assays$Spatial@VariableFeatures
  cat("高变基因数量:", length(hvg), "\n")
  # 选择前100个高变基因
  n_genes <- min(100, length(hvg))
  hvg_top <- hvg[1:n_genes]
  
  # 找到这些基因在scale.data中的索引
  gene_names <- rownames(expr_data)
  hvg_idx <- match(hvg_top, gene_names)
  hvg_idx <- hvg_idx[!is.na(hvg_idx)]
  
  if (length(hvg_idx) > 20) {
    expr_hvg <- expr_data[hvg_idx, ]
    cat("使用", nrow(expr_hvg), "个高变基因\n")
  } else {
    expr_hvg <- expr_data
    cat("高变基因不足，使用所有基因\n")
  }
} else {
  # 没有高变基因信息，选择方差最大的基因
  gene_vars <- apply(expr_data, 1, var)
  top_idx <- order(gene_vars, decreasing = TRUE)[1:100]
  expr_hvg <- expr_data[top_idx, ]
  cat("选择方差最大的100个基因\n")
}

# 使用前100个细胞测试
M_test <- 100
X_test <- as.matrix(expr_hvg[, 1:M_test])
cat("测试数据维度:", dim(X_test), "\n")

# ==============================================================================
# ARNI分析
# ==============================================================================
cat("\n=== ARNI分析 ===\n")

# 计算导数
DX <- X_test[, 2:M_test] - X_test[, 1:(M_test-1)]
X_reduced <- X_test[, 1:(M_test-1)]

cat("DX维度:", dim(DX), "\n")
cat("X_reduced维度:", dim(X_reduced), "\n")

# 目标变量
target <- 1
DX_target <- DX[target, , drop = FALSE]

# 候选源变量
remaining <- setdiff(1:nrow(X_reduced), target)

cat("目标变量:", target, "\n")
cat("候选源变量数量:", length(remaining), "\n")

# 贪婪选择
selected <- c()
max_interactions <- 5

costs_overall <- numeric(max_interactions)
sources_selected <- list()

for (k in 1:max_interactions) {
  cat("\n搜索第", k, "个交互...\n")
  
  best_cost <- Inf
  best_j <- -1
  
  for (j in remaining) {
    # 当前候选集
    current_idx <- c(selected, j)
    
    # 构建回归
    df_current <- data.frame(DX = as.vector(DX_target))
    for (idx in current_idx) {
      df_current[[paste0("X", idx)]] <- as.vector(X_reduced[idx, ])
    }
    
    formula_current <- as.formula(paste0("DX ~ ", paste0("+ X", current_idx, collapse = "")))
    
    tryCatch({
      lm_current <- lm(formula_current, data = df_current)
      DX_est <- predict(lm_current)
      DIFF <- as.vector(DX_target) - DX_est
      cost <- sqrt(mean(DIFF^2))
      
      if (cost < best_cost) {
        best_cost <- cost
        best_j <- j
      }
    }, error = function(e) {})
  }
  
  if (best_j > 0) {
    selected <- c(selected, best_j)
    remaining <- setdiff(remaining, best_j)
    costs_overall[k] <- best_cost
    sources_selected[[k]] <- best_j
    cat("  选中基因", best_j, ", 代价:", best_cost, "\n")
  } else {
    costs_overall[k] <- Inf
    break
  }
}

# ==============================================================================
# 结果
# ==============================================================================
cat("\n=== 结果 ===\n")
cat("选择的基因索引:", selected, "\n")
cat("代价曲线:", costs_overall[1:length(selected)], "\n")

# 计算R²改善
if (length(selected) > 0) {
  # 零模型（只用均值）
  DX_mean <- mean(DX_target)
  SS_tot <- sum((DX_target - DX_mean)^2)
  
  # 最终模型
  df_final <- data.frame(DX = as.vector(DX_target))
  for (idx in selected) {
    df_final[[paste0("X", idx)]] <- as.vector(X_reduced[idx, ])
  }
  formula_final <- as.formula(paste0("DX ~ ", paste0("+ X", selected, collapse = "")))
  lm_final <- lm(formula_final, data = df_final)
  DX_est_final <- predict(lm_final)
  SS_res <- sum((DX_target - DX_est_final)^2)
  R2_final <- 1 - SS_res / SS_tot
  
  cat("最终模型R²:", R2_final, "\n")
  
  if (R2_final > 0.1) {
    cat("\n✓ 成功找到有效的网络交互!\n")
  } else {
    cat("\n✗ R²仍然较低，可能需要更多优化\n")
  }
}

# ==============================================================================
# 比较：基因 vs PCA
# ==============================================================================
cat("\n=== 基因 vs PCA 比较 ===\n")

# 用PCA数据测试
pca_data <- skin1@reductions$pca@cell.embeddings[subset_idx, 1:20, drop = FALSE]
X_pca <- t(pca_data)[, 1:M_test]

DX_pca <- X_pca[, 2:M_test] - X_pca[, 1:(M_test-1)]
X_pca_reduced <- X_pca[, 1:(M_test-1)]

# 目标变量
DX_target_pca <- DX_pca[target, , drop = FALSE]
remaining_pca <- setdiff(1:20, target)

selected_pca <- c()
for (j in remaining_pca) {
  df_pca <- data.frame(DX = as.vector(DX_target_pca), X = as.vector(X_pca_reduced[j, ]))
  lm_pca <- lm(DX ~ X, data = df_pca)
  DX_est_pca <- predict(lm_pca)
  cost_pca <- sqrt(mean((as.vector(DX_target_pca) - DX_est_pca)^2))
  
  if (j == remaining_pca[1]) {
    cat("PCA - 第一个源变量代价:", cost_pca, "\n")
  }
}

# 基因数据
DX_target_gene <- DX[target, , drop = FALSE]
remaining_gene <- setdiff(1:nrow(X_reduced), target)

selected_gene <- c()
for (j in 1:min(20, length(remaining_gene))) {
  idx <- remaining_gene[j]
  df_gene <- data.frame(DX = as.vector(DX_target_gene), X = as.vector(X_reduced[idx, ]))
  lm_gene <- lm(DX ~ X, data = df_gene)
  DX_est_gene <- predict(lm_gene)
  cost_gene <- sqrt(mean((as.vector(DX_target_gene) - DX_est_gene)^2))
  
  if (j == 1) {
    cat("Gene - 第一个源变量代价:", cost_gene, "\n")
  }
}

cat("\n结论：比较原始基因和PCA数据的效果\n")
