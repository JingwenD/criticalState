#!/usr/bin/env Rscript
# ARNI核心逻辑测试 - 真正正确的版本
# 公式: dx_i/dt ≈ A_ij * x_j

library(Seurat)
library(pracma)

# 加载数据
#skin1 <- readRDS("skin1.rds")
sample_id <- "lsm"
subset_idx <- skin1$orig.ident == sample_id
pca_data <- skin1@reductions$pca@cell.embeddings[subset_idx, 1:20, drop = FALSE]
X <- t(pca_data)  # 20 x 26924

# 使用前100个样本
M_test <- 100
X_test <- X[, 1:M_test]

# ==============================================================================
# Step 1: 计算导数
# ==============================================================================
DX <- X_test[, 2:M_test] - X_test[, 1:(M_test-1)]  # 20 x 99
X_reduced <- X_test[, 1:(M_test-1)]  # 20 x 99

cat("数据维度:\n")
cat("  X (原始值):", dim(X_test), "\n")
cat("  X_reduced (前M-1个):", dim(X_reduced), "\n")
cat("  DX (导数):", dim(DX), "\n")

# ==============================================================================
# Step 2: 正确理解ARNI公式
# ==============================================================================
cat("\n=== ARNI的正确理解 ===\n")
cat("公式: dx_i/dt = Σ_j A_ij × x_j(t)\n")
cat("\n含义:\n")
cat("  dx_i: 目标变量i的导数 (1 x 99)\n")
cat("  x_j: 源变量j的值 (1 x 99)\n")
cat("  A_ij: 交互系数 (标量)\n")
cat("\n")

# ==============================================================================
# Step 3: 测试单个源变量
# ==============================================================================
cat("=== 测试单个源变量 ===\n")

target <- 2  # 目标变量索引（0-based）
source <- 1  # 源变量索引（0-based）

DX_target <- DX[target + 1, , drop = FALSE]  # 1 x 99
X_source <- X_reduced[source + 1, , drop = FALSE]  # 1 x 99

cat("目标变量", target + 1, "的导数维度:", dim(DX_target), "\n")
cat("源变量", source + 1, "的值维度:", dim(X_source), "\n")

# 简单线性回归: DX_target ≈ A * X_source
# 这相当于用源变量的值拟合目标变量的导数
cor_val <- cor(as.vector(DX_target), as.vector(X_source))
cat("DX_target和X_source的相关系数:", cor_val, "\n")

# 线性回归
df <- data.frame(
  DX = as.vector(DX_target),
  X_source = as.vector(X_source)
)
lm_simple <- lm(DX ~ X_source, data = df)
cat("简单线性回归R²:", summary(lm_simple)$r.squared, "\n")
cat("系数A:", coef(lm_simple), "\n")

# ==============================================================================
# Step 4: 测试多个源变量的组合
# ==============================================================================
cat("\n=== 测试多个源变量组合 ===\n")

target <- 2
DX_target <- DX[target + 1, , drop = FALSE]

# 使用多个源变量进行多元回归
# DX_target ≈ A1*X1 + A2*X2 + ...
candidates <- 1:20
candidates <- candidates[candidates != target]

cat("目标变量:", target + 1, "\n")
cat("候选源变量数量:", length(candidates), "\n")

# 多元回归
df_multi <- data.frame(DX = as.vector(DX_target))
for (j in 1:min(10, length(candidates))) {
  idx <- candidates[j]
  df_multi[[paste0("X", j)]] <- as.vector(X_reduced[idx + 1, ])
}

formula_str <- paste0("DX ~ ", paste0("X", 1:min(10, length(candidates)), collapse = "+"))
cat("回归公式:", formula_str, "\n")

lm_multi <- lm(as.formula(formula_str), data = df_multi)
cat("多元回归R²:", summary(lm_multi)$r.squared, "\n")
cat("调整后R²:", summary(lm_multi)$adj.r.squared, "\n")
cat("系数:\n")
print(coef(lm_multi))

# ==============================================================================
# Step 5: 真正的ARNI贪婪选择
# ==============================================================================
cat("\n=== ARNI贪婪选择算法 ===\n")
cat("策略：逐个添加源变量，选择能最大程度降低残差的变量\n")

target <- 2
DX_target <- DX[target + 1, , drop = FALSE]
remaining <- setdiff(1:20, target + 1)
selected <- c()
max_interactions <- 5

costs_overall <- numeric(max_interactions)

for (k in 1:max_interactions) {
  cat("\n搜索第", k, "个交互...\n")
  
  best_cost <- Inf
  best_j <- -1
  
  for (j in remaining) {
    # 当前候选集
    current_idx <- c(selected, j)
    
    # 构建回归矩阵
    df_current <- data.frame(DX = as.vector(DX_target))
    for (idx in current_idx) {
      # idx是0-based，所以idx+1是1-based索引
      mat_idx <- idx + 1
      if (mat_idx <= nrow(X_reduced)) {
        df_current[[paste0("X", idx)]] <- as.vector(X_reduced[mat_idx, ])
      }
    }
    
    formula_current <- as.formula(paste0("DX ~ ", paste0("+ X", current_idx, collapse = "")))
    
    tryCatch({
      lm_current <- lm(formula_current, data = df_current)
      
      # 计算残差
      DX_est <- predict(lm_current)
      DIFF <- as.vector(DX_target) - DX_est
      cost <- sqrt(mean(DIFF^2))
      
      if (cost < best_cost) {
        best_cost <- cost
        best_j <- j
      }
    }, error = function(e) {
      # 如果回归失败，跳过
    })
  }
  
  if (best_j > 0) {
    selected <- c(selected, best_j)
    remaining <- setdiff(remaining, best_j)
    costs_overall[k] <- best_cost
    cat("  选中变量", best_j, ", 代价:", best_cost, "\n")
  } else {
    costs_overall[k] <- Inf
    break
  }
}

# 结果
cat("\n=== ARNI结果 ===\n")
cat("选择的源变量:", selected, "\n")
cat("对应代价:", costs_overall[1:length(selected)], "\n")

# 判断是否找到有效交互
if (length(selected) > 0 && costs_overall[1] < 1.0) {
  cat("\n✓ 成功找到网络交互!\n")
} else {
  cat("\n✗ ARNI未能找到有效的网络交互\n")
  cat("可能原因:\n")
  cat("  1. PCA数据中确实没有因果关系\n")
  cat("  2. 空间数据不适合用时间序列方法分析\n")
  cat("  3. 需要使用原始基因表达数据而非PCA\n")
}

# ==============================================================================
# Step 6: 比较不同目标变量
# ==============================================================================
cat("\n=== 比较不同目标变量的ARNI结果 ===\n")

results_summary <- data.frame(
  target = integer(),
  best_source = integer(),
  min_cost = numeric()
)

for (target_idx in 1:min(10, 20)) {
  DX_t <- DX[target_idx, , drop = FALSE]
  remaining_t <- setdiff(1:20, target_idx)
  
  min_cost_t <- Inf
  best_source_t <- -1
  
  for (j in remaining_t) {
    df_t <- data.frame(
      DX = as.vector(DX_t),
      X = as.vector(X_reduced[j, ])
    )
    
    lm_t <- lm(DX ~ X, data = df_t)
    DX_est_t <- predict(lm_t)
    cost_t <- sqrt(mean((as.vector(DX_t) - DX_est_t)^2))
    
    if (cost_t < min_cost_t) {
      min_cost_t <- cost_t
      best_source_t <- j
    }
  }
  
  results_summary <- rbind(results_summary, data.frame(
    target = target_idx,
    best_source = best_source_t,
    min_cost = min_cost_t
  ))
}

print(results_summary)

cat("\n平均最小代价:", mean(results_summary$min_cost), "\n")
