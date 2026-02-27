#!/usr/bin/env Rscript
# ARNI Analysis Results Visualization - Complete Version (with Spatial Transcriptomics Integration)
# Fixed version with detailed annotations for each plot
#
# ==============================================================================
# OVERVIEW - 图表演释总览
# ==============================================================================
# 本脚本生成7类可视化图表，帮助理解ARNI推断的基因调控网络：
#
# 1. 学习曲线 (Learning Curves) - 评估每个目标需要多少个调控源
# 2. 网络热图 (Network Heatmap) - 展示基因间交互强度矩阵
# 3. 网络图 (Network Graph) - 图论视角的网络拓扑结构
# 4. 空间网络可视化 (Spatial Network) - 结合空间位置分析网络
# 5. Seurat原生图 (Seurat Plots) - 生成高质量空间可视化代码
# 6. 网络统计 (Network Statistics) - 网络拓扑定量特征
# 7. 表达热力图 (Expression Heatmap) - 基因时序表达模式
#
# ==============================================================================

library(Seurat)
library(ggplot2)
library(igraph)
library(pheatmap)
library(viridis)
library(gridExtra)

# ==============================================================================
# Helper Functions
# ==============================================================================

get_sample_colors <- function(n_samples) {
  colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
              "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5")
  if (n_samples <= length(colors)) {
    return(colors[1:n_samples])
  } else {
    return(colorRampPalette(colors)(n_samples))
  }
}

# ==============================================================================
# 1. Learning Curves Visualization (学习曲线可视化)
# ==============================================================================
#
# 【图表目的】
# 学习曲线展示了ARNI算法在选择不同数量交互源时，模型拟合代价的变化趋势。
# 用于评估：
#   1) 每个目标基因需要多少个源基因交互才能充分解释其变化
#   2) 增加更多交互源是否能显著降低模型代价
#   3) 识别最优的交互数量（最小代价点）
#
# 【如何解读】
#   - X轴：选择的交互源数量（1, 2, 3, ...）
#   - Y轴：模型拟合的残差代价（Cost），越低表示拟合越好
#   - 红色点：标记每个目标的最优交互数量（代价最低点）
#   - 曲线下降趋势：说明加入更多交互源能更好地解释目标变化
#   - 曲线趋于平缓：说明再增加交互源的边际效益递减
#
# 【关键指标含义】
#   - optimal_interactions：每个目标基因建议使用的交互数量
#   - 代价下降幅度：反映网络关系的强弱，下降越明显说明交互越重要
#   - 曲线末端代价：如果仍很高，说明可能需要更多交互或数据噪声大
#
# 【生物学解释】
#   - 陡峭下降：目标基因的变化可由少数几个源基因充分解释
#   - 平缓下降：目标基因的变化需要多个源基因协同解释，可能涉及复杂调控
#   - 多个局部最低：可能存在不同尺度的调控模式
#
plot_learning_curves <- function(results, sample_id, output_dir = "ARNI_visualization") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  pdf(file = paste0(output_dir, "/learning_curves_", sample_id, ".pdf"), width = 14, height = 10)
  
  tryCatch({
    par(mar = c(5, 4, 4, 2) + 0.1)
  }, error = function(e) {})
  
  n_display <- min(20, length(results))
  n_cols <- 5
  n_rows <- ceiling(n_display / n_cols)
  
  par(oma = c(1, 1, 1, 1), mfrow = c(n_rows, n_cols), mar = c(3, 3, 2, 1))
  
  for (i in 1:min(length(results), 30)) {
    result <- results[[i]]
    
    if (!is.null(result$learning_curve) && is.data.frame(result$learning_curve)) {
      lc <- result$learning_curve
      
      valid_idx <- is.finite(lc$cost)
      if (sum(valid_idx) > 0) {
        lc_valid <- lc[valid_idx, ]
        
        plot(lc_valid$n_interactions, lc_valid$cost, 
             type = "b", pch = 19, col = "steelblue",
             xlab = "Interactions", ylab = "Cost",
             main = paste("Target", i),
             cex.main = 0.8)
        
        if (!is.null(result$optimal_interactions)) {
          opt_k <- result$optimal_interactions
          if (opt_k <= nrow(lc_valid)) {
            points(lc_valid$n_interactions[opt_k], lc_valid$cost[opt_k], 
                   pch = 19, col = "red", cex = 2)
          }
        }
      }
    }
  }
  
  dev.off()
  cat("Learning curves saved:", paste0(output_dir, "/learning_curves_", sample_id, ".pdf"), "\n")
}

# ==============================================================================
# 2. Network Heatmap (网络热图)
# ==============================================================================
#
# 【图表目的】
# 以热图形式展示ARNI推断的基因间相互作用矩阵。
# 直观显示哪些基因对之间存在调控关系及其强度。
#
# 【如何解读】
#   - 行（Row）：目标基因（被调控的基因）
#   - 列（Column）：源基因（调控其他基因的基因）
#   - 颜色深浅：交互系数大小，深色表示强交互
#   - 聚类分析：行和列会根据相似性自动聚类，揭示功能相关的基因模块
#
# 【关键指标含义】
#   - 矩阵元素值 > 0：源基因表达升高时，目标基因表达趋于升高（激活关系）
#   - 矩阵元素值 < 0：源基因表达升高时，目标基因表达趋于降低（抑制关系）
#   - 绝对值越大：交互作用越强
#   - 聚类结果：揭示协同变化的基因群，可能属于同一生物学通路
#
# 【注意事项】
#   - 如果矩阵太大，默认只显示前50个基因
#   - 聚类可以帮助发现基因模块，但结果受限于显示的基因子集
#   - 对角线值可能较高（自调控），需要结合生物学知识判断
#
# 【生物学解释】
#   - 正交互模块：一起上调的基因群，可能受共同上游因子调控
#   - 负交互模块：此消彼长的基因对，可能存在竞争关系
#   - 聚类边界清晰的模块：功能独立的基因集合
#
plot_network_heatmap <- function(network_matrix, sample_id, output_dir = "ARNI_visualization") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  pdf(file = paste0(output_dir, "/network_heatmap_", sample_id, ".pdf"), width = 10, height = 8)
  
  display_matrix <- network_matrix
  if (nrow(network_matrix) > 50) {
    display_matrix <- network_matrix[1:50, 1:50]
  }
  
  pheatmap(display_matrix,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           color = colorRampPalette(c("white", "steelblue", "darkblue"))(100),
           show_rownames = FALSE,
           show_colnames = FALSE,
           main = paste("ARNI Network Matrix -", sample_id),
           fontsize = 10)
  
  dev.off()
  cat("Network heatmap saved:", paste0(output_dir, "/network_heatmap_", sample_id, ".pdf"), "\n")
}

# ==============================================================================
# 3. Network Graph (网络图 - 使用igraph)
# ==============================================================================
#
# 【图表目的】
# 使用图论可视化方法展示基因调控网络。
# 将基因表示为节点，交互关系表示为有向边。
#
# 【如何解读】
#   - 节点（圆点）：代表基因
#   - 边（箭头）：表示调控方向，从源基因指向目标基因
#   - 节点大小：基于总度（入度+出度），越大表示越可能是Hub基因
#   - 节点颜色：基于入度（被调控次数），颜色越深表示被越多基因调控
#
# 【关键指标含义】
#   - Hub基因（大的节点）：网络中高度连接的基因，可能具有重要生物学功能
#   - 高度数节点：被多个基因调控的靶点，可能处于信号通路的下游
#   - 出度高的节点：影响多个下游基因的调控因子
#   - 边的方向：从A指向B表示A的变化导致B的变化
#
# 【网络特征解释】
#   - 密集区域：可能代表功能相关的基因模块
#   - 孤立节点：与其他基因交互较少的基因
#   - 层级结构：某些基因可能处于调控层级的中枢位置
#   - 双向边：可能存在相互调控的反馈回路
#
# 【生物学意义】
#   - 中央Hub基因：潜在的 master regulator，可作为 biomarker
#   - 边缘基因：特异性调控特定生物学过程的执行者
#   - 密集子网络：功能相关的基因群，可能属于同一 pathway
#
plot_network_graph <- function(network_matrix, sample_id, seurat_obj = NULL, output_dir = "ARNI_visualization") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  pdf(file = paste0(output_dir, "/network_graph_", sample_id, ".pdf"), width = 14, height = 14)
  
  # 获取基因名称
  if (!is.null(seurat_obj)) {
    gene_names <- rownames(seurat_obj@assays$Spatial$scale.data)
    gene_names <- gene_names[1:nrow(network_matrix)]
  } else {
    gene_names <- paste0("Gene_", 1:nrow(network_matrix))
  }
  
  # 如果基因名称数量超过矩阵行数，进行调整
  if (length(gene_names) < nrow(network_matrix)) {
    additional_names <- paste0("Gene_", (length(gene_names) + 1):nrow(network_matrix))
    gene_names <- c(gene_names, additional_names)
  }
  
  g <- graph_from_adjacency_matrix(network_matrix, mode = "directed", weighted = NULL)
  
  # 设置节点名称为基因名称
  V(g)$name <- gene_names
  
  in_degree <- colSums(network_matrix)
  out_degree <- rowSums(network_matrix)
  
  node_size <- (in_degree + out_degree) + 3
  node_color <- in_degree
  
  layout <- layout_with_fr(g)
  
  # 根据基因数量动态调整字体大小
  n_genes <- nrow(network_matrix)
  if (n_genes <= 20) {
    label_cex <- 0.8
  } else if (n_genes <= 50) {
    label_cex <- 0.6
  } else if (n_genes <= 100) {
    label_cex <- 0.4
  } else {
    label_cex <- 0.25
  }
  
  plot(g, layout = layout,
       vertex.size = node_size,
       vertex.color = node_color,
       vertex.label = gene_names,
       vertex.label.cex = label_cex,
       vertex.label.color = "black",
       vertex.label.dist = 0.5,
       vertex.label.degree = -pi/4,
       edge.arrow.size = 0.3,
       edge.width = 0.5,
       edge.color = "gray50",
       main = paste("ARNI Network Graph -", sample_id),
       xlab = "", ylab = "")
  
  dev.off()
  cat("Network graph saved:", paste0(output_dir, "/network_graph_", sample_id, ".pdf"), "\n")
}

# ==============================================================================
# 4. Spatial Transcriptomics Visualization (空间转录组可视化 - 核心功能)
# ==============================================================================
#
# 【图表目的】
# 将ARNI推断的网络交互结果与空间转录组的位置信息结合，
# 从空间维度分析基因调控模式。
#
# 【9个面板详细解读】
#
# Panel 1: 空间分布（按伪时间着色）
#   - 将Y坐标离散化为时间箱，构建伪时间轨迹
#   - 每个点的颜色代表其在伪时间轴上的位置
#   - 用于观察基因表达是否沿空间位置呈现时间动态
#   - 颜色梯度反映空间上的时间进展
#
# Panel 2: 样本空间分布
#   - 显示样本的空间结构，验证细胞在空间上的分布
#   - 确认数据质量：细胞应均匀分布，无明显聚集
#
# Panel 3: 伪时间分布
#   - 显示各时间箱中的细胞数量
#   - 检验时间箱划分是否均匀
#   - 不均匀分布可能反映空间结构特征
#
# Panel 4: 网络矩阵热图（前20节点）
#   - 网络交互矩阵的简化视图
#   - 快速识别强交互的基因对
#   - 对角线：自调控关系
#
# Panel 5: 入度分布
#   - 每个基因被多少其他基因调控的分布
#   - 峰值位置反映网络的大多数节点的入度特征
#   - 高入度基因：被多因素调控的重要靶点
#
# Panel 6: 出度分布
#   - 每个基因调控多少其他基因的分布
#   - 高出度节点可能是关键的调控因子
#   - 少数高出度节点主导网络信息流
#
# Panel 7: Top 15 Hub节点
#   - 按总度排序的前15个基因
#   - 这些基因在网络中高度连接
#   - 可能作为 biomarker 或治疗靶点
#
# Panel 8: 学习曲线示例（目标1）
#   - 展示第一个目标基因的ARNI学习过程
#   - 观察交互选择如何影响模型拟合
#   - 识别最优交互数量
#
# Panel 9: 网络统计摘要
#   - 节点数：网络中的基因数量
#   - 边数：推断的交互关系数量
#   - 密度：实际边数与可能的最大边数之比
#   - 平均入度/出度：网络连接性的平均特征
#   - 最大入度/出度：网络中最具影响力的节点
#
plot_spatial_with_network <- function(seurat_obj, results, network_matrix, sample_id,
                                       output_dir = "ARNI_visualization") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  pdf(file = paste0(output_dir, "/spatial_network_", sample_id, ".pdf"), width = 16, height = 14)
  
  tryCatch({
    par(mar = c(5, 4, 4, 2) + 0.1)
  }, error = function(e) {})
  
  par(oma = c(2, 2, 2, 2), mar = c(4, 4, 3, 1))
  
  subset_idx <- seurat_obj$orig.ident == sample_id
  n_cells <- sum(subset_idx)
  
  coords <- data.frame(
    cell_id = rownames(seurat_obj@meta.data[subset_idx, ]),
    x = seurat_obj@meta.data[subset_idx, "x"],
    y = seurat_obj@meta.data[subset_idx, "y"]
  )
  
  n_time_bins <- 50
  y_quantiles <- seq(0, 1, length.out = n_time_bins + 1)
  y_breaks <- quantile(coords$y, probs = y_quantiles)
  y_breaks[1] <- y_breaks[1] - 0.1
  y_breaks[length(y_breaks)] <- y_breaks[length(y_breaks)] + 0.1
  
  time_bins <- cut(coords$y, breaks = y_breaks, labels = FALSE, include.lowest = TRUE)
  coords$time_bin <- time_bins
  
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(2, 2, 2, 2))
  
  # Panel 1
  plot(coords$x, coords$y, type = "n", xlab = "X", ylab = "Y", 
       main = "Spatial Distribution (Pseudo-time colored)")
  points(coords$x, coords$y, col = coords$time_bin, pch = 19, cex = 0.3)
  
  # Panel 2
  plot(coords$x, coords$y, col = "steelblue", pch = 19, cex = 0.3,
       xlab = "X", ylab = "Y", main = paste("Sample", sample_id, "Spatial Distribution"))
  
  # Panel 3
  hist(time_bins, breaks = n_time_bins, col = "steelblue",
       xlab = "Time Bin", ylab = "Cell Count",
       main = "Pseudo-time Distribution")
  
  # Panel 4
  display_matrix <- network_matrix
  if (nrow(network_matrix) > 30) {
    display_matrix <- network_matrix[1:30, 1:30]
  }
  pheatmap(display_matrix[1:min(20, nrow(display_matrix)), 1:min(20, ncol(display_matrix))],
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           color = colorRampPalette(c("white", "steelblue", "darkblue"))(100),
           show_rownames = TRUE,
           show_colnames = TRUE,
           main = "Network Matrix (Top 20 nodes)",
           fontsize = 8)
  
  # Panel 5
  in_degree <- colSums(network_matrix)
  barplot(table(in_degree), col = viridis::viridis(length(table(in_degree))),
          xlab = "In-degree", ylab = "Node Count",
          main = "In-degree Distribution")
  
  # Panel 6
  out_degree <- rowSums(network_matrix)
  barplot(table(out_degree), col = viridis::viridis(length(table(out_degree))),
          xlab = "Out-degree", ylab = "Node Count",
          main = "Out-degree Distribution")
  
  # Panel 7
  total_degree <- in_degree + out_degree
  top_hubs <- head(sort(total_degree, decreasing = TRUE), 15)
  barplot(top_hubs, col = viridis::viridis(15),
          xlab = "Node", ylab = "Degree",
          main = "Top 15 Hub Nodes")
  
  # Panel 8
  if (length(results) >= 1 && !is.null(results[[1]]$learning_curve)) {
    lc <- results[[1]]$learning_curve
    valid_idx <- is.finite(lc$cost)
    if (sum(valid_idx) > 0) {
      lc_valid <- lc[valid_idx, ]
      plot(lc_valid$n_interactions, lc_valid$cost, 
           type = "b", pch = 19, col = "steelblue",
           xlab = "Interactions", ylab = "Cost",
           main = "Learning Curve Example (Target 1)")
    }
  }
  
  # Panel 9
  g <- graph_from_adjacency_matrix(network_matrix, mode = "directed")
  n_edges <- sum(network_matrix > 0)
  n_nodes <- nrow(network_matrix)
  density <- n_edges / (n_nodes * (n_nodes - 1))
  
  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
  plot.new()
  text(0.5, 0.9, paste("ARNI Network Statistics -", sample_id), cex = 1.5, font = 2)
  text(0.5, 0.75, paste("Nodes:", n_nodes), cex = 1.2)
  text(0.5, 0.65, paste("Edges:", n_edges), cex = 1.2)
  text(0.5, 0.55, paste("Density:", round(density, 4)), cex = 1.2)
  text(0.5, 0.45, paste("Avg In-degree:", round(mean(in_degree), 2)), cex = 1.2)
  text(0.5, 0.35, paste("Avg Out-degree:", round(mean(out_degree), 2)), cex = 1.2)
  
  if (sum(is.finite(in_degree[in_degree > 0])) > 0) {
    text(0.5, 0.25, paste("Max In-degree:", max(in_degree[in_degree > 0])), cex = 1.2)
    text(0.5, 0.15, paste("Max Out-degree:", max(out_degree[out_degree > 0])), cex = 1.2)
  }
  
  par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))
  
  dev.off()
  cat("Spatial network visualization saved:", paste0(output_dir, "/spatial_network_", sample_id, ".pdf"), "\n")
}

# ==============================================================================
# 5. Seurat Native Visualization (使用Seurat原生空间可视化 - 推荐)
# ==============================================================================
#
# 【图表目的】
# 生成高质量的Seurat空间可视化代码，
# 便于用户在RStudio中进行更精细的交互式可视化。
#
# 【生成的内容】
#   1) 时间箱空间分布图：展示伪时间在空间上的分布
#   2) Hub基因空间表达图：Top基因在空间上的表达模式
#   3) 组合图：时间箱与平均表达的叠加
#
# 【如何解读】
#   - 生成的R脚本需要复制到RStudio中运行
#   - 可以进一步调整颜色、大小、布局等参数
#   - 适合生成publication级别的 figures
#
# 【关键基因选择】
#   - 基于网络拓扑选择Top 5 Hub基因
#   - 这些基因在网络中高度连接，可能具有重要生物学意义
#   - 观察这些基因的空间表达模式是否与其网络位置相关
#
# 【生物学应用】
#   - 空间共表达分析：Hub基因是否在相同区域高表达
#   - 空间邻域分析：邻居细胞的基因表达是否相似
#   - 空间异质性：不同区域的调控网络是否不同
#
prepare_seurat_plots <- function(seurat_obj, results, network_matrix, sample_id,
                                  output_dir = "ARNI_visualization") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  subset_obj <- subset(seurat_obj, orig.ident == sample_id)
  
  coords <- data.frame(
    x = subset_obj@meta.data$x,
    y = subset_obj@meta.data$y
  )
  
  n_bins <- 50
  y_quantiles <- seq(0, 1, length.out = n_bins + 1)
  y_breaks <- quantile(coords$y, probs = y_quantiles)
  y_breaks[1] <- y_breaks[1] - 0.1
  y_breaks[length(y_breaks)] <- y_breaks[length(y_breaks)] + 0.1
  
  time_bins <- cut(coords$y, breaks = y_breaks, labels = FALSE, include.lowest = TRUE)
  subset_obj$time_bin <- time_bins
  
  if ("VariableFeatures" %in% names(subset_obj@assays$Spatial)) {
    hvg <- subset_obj@assays$Spatial@VariableFeatures
    features <- hvg[1:min(5, length(hvg))]
  } else {
    scale_data <- subset_obj@assays$Spatial$scale.data
    gene_names <- rownames(scale_data)
    top_genes_idx <- head(order(colSums(network_matrix), decreasing = TRUE), 5)
    top_genes_idx <- top_genes_idx[top_genes_idx <= length(gene_names)]
    features <- gene_names[top_genes_idx]
  }
  
  features_str <- paste0('c("', paste(features, collapse = '", "'), '")')
  
  plot_code_lines <- c(
    '# ============================================',
    '# Seurat Spatial Visualization Code',
    '# Recommended to run in RStudio',
    '# ============================================',
    '',
    paste0('# Load filtered Seurat object'),
    paste0('subset_obj <- readRDS("', output_dir, '/subset_', sample_id, '.rds")'),
    '',
    paste0('# Define genes to visualize'),
    paste0('features <- ', features_str),
    '',
    paste0('# Save temporary data'),
    paste0('saveRDS(subset_obj, file = "', output_dir, '/subset_', sample_id, '.rds")'),
    '',
    '# 1. Time bin spatial visualization',
    'SpatialFeaturePlot(subset_obj, features = "time_bin", ', 
    '                   pt.size = 0.5, crop = TRUE) +',
    '  theme(legend.position = "right")',
    '',
    paste0('ggsave("', output_dir, '/spatial_timebin_', sample_id, '.pdf", width = 10, height = 8)'),
    '',
    '# 2. High variable gene expression',
    'for (gene in features) {',
    '  if (gene %in% rownames(subset_obj)) {',
    '    SpatialFeaturePlot(subset_obj, features = gene,',
    '                       pt.size = 0.5, crop = TRUE) +',
    '      ggtitle(paste("Gene Expression -", gene)) +',
    '      theme(legend.position = "right")',
    '    ',
    paste0('    ggsave("', output_dir, '/spatial_gene_", gene, "_', sample_id, '.pdf", '),
    '           width = 10, height = 8)',
    '  }',
    '}',
    '',
    '# 3. Combined plot',
    'SpatialPlot(list(',
    '  time_bin = subset_obj$time_bin,',
    '  expr = colMeans(subset_obj@assays$Spatial$scale.data)',
    '), pt.size = 0.5) +',
    '  theme(legend.position = "right")'
  )
  plot_code <- paste(plot_code_lines, collapse = '\n')
  
  writeLines(plot_code, con = paste0(output_dir, "/seurat_spatial_plots_", sample_id, ".R"))
  
  cat("Seurat spatial visualization code saved:", paste0(output_dir, "/seurat_spatial_plots_", sample_id, ".R"), "\n")
  cat("  Recommended to run this script in RStudio for high-quality images\n")
}

# ==============================================================================
# 6. Network Statistics Summary (网络统计汇总)
# ==============================================================================
#
# 【图表目的】
# 汇总统计ARNI推断的网络的关键特征，
# 定量描述网络的拓扑属性。
#
# 【6个面板详细解读】
#
# Panel 1: 最优交互数分布
#   - 每个目标基因最优需要的源基因数量分布
#   - 峰值位置：大多数目标基因的典型交互数量
#   - 宽分布：不同目标基因的调控复杂度差异大
#   - 窄分布：网络结构相对均匀
#   - 生物学意义：反映调控网络的局部复杂度
#
# Panel 2: 第一个交互的代价分布
#   - 仅使用最优源基因时的模型代价分布
#   - 低代价：存在主导性的调控关系
#   - 高代价：需要多个源基因协同解释目标变化
#   - 可用于识别"简单"vs"复杂"的调控目标
#
# Panel 3: 入度分布（非零）
#   - 只考虑有入边的节点的入度分布
#   - 揭示被调控基因的分布特征
#   - 幂律分布：少数基因被大量调控，多数基因被少量调控
#   - 峰值在低值：大多数基因只有少数几个调控者
#
# Panel 4: 出度分布（非零）
#   - 只考虑有出边的节点的出度分布
#   - 揭示调控基因的影响力分布
#   - 少数高出度基因可能是关键调控因子
#   - 出度分布反映信息流的集中程度
#
# Panel 5: Top 15 Hub节点
#   - 按总度排序的前15个基因
#   - 这些基因是网络的中心节点
#   - 可能作为 biomarker 或治疗靶点
#   - Hub基因往往参与多个生物学过程
#
# Panel 6: 网络规模
#   - 点图展示节点数vs边数
#   - 标注网络密度
#   - 密集网络：边数接近节点数的平方
#   - 稀疏网络：大多数可能的交互不存在
#   - 密度反映网络的连通性
#
plot_network_statistics <- function(results, network_matrix, sample_id,
                                     output_dir = "ARNI_visualization") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  pdf(file = paste0(output_dir, "/network_statistics_", sample_id, ".pdf"), width = 14, height = 10)
  
  tryCatch({
    par(mar = c(5, 4, 4, 2) + 0.1)
  }, error = function(e) {})
  
  par(oma = c(2, 2, 2, 2), mfrow = c(2, 3), mar = c(4, 4, 3, 1))
  
  optimal_interactions <- sapply(results, function(r) r$optimal_interactions)
  barplot(table(optimal_interactions), 
          col = viridis::viridis(length(table(optimal_interactions))),
          xlab = "Optimal Interactions", ylab = "Target Count",
          main = "Optimal Interactions Distribution")
  
  first_costs <- sapply(results, function(r) {
    if (!is.null(r$cost[1]) && is.finite(r$cost[1])) {
      return(r$cost[1])
    }
    return(NA)
  })
  first_costs <- first_costs[!is.na(first_costs)]
  if (length(first_costs) > 0) {
    hist(first_costs, breaks = 20, col = "coral",
         xlab = "Cost (k=1)", ylab = "Target Count",
         main = "First Interaction Cost Distribution")
  }
  
  in_degree <- colSums(network_matrix)
  barplot(table(in_degree[in_degree > 0]), 
          col = "steelblue",
          xlab = "In-degree", ylab = "Node Count",
          main = "In-degree Distribution (Non-zero)")
  
  out_degree <- rowSums(network_matrix)
  barplot(table(out_degree[out_degree > 0]), 
          col = "coral",
          xlab = "Out-degree", ylab = "Node Count",
          main = "Out-degree Distribution (Non-zero)")
  
  total_degree <- in_degree + out_degree
  top_hubs <- head(sort(total_degree, decreasing = TRUE), 15)
  barplot(top_hubs, 
          col = viridis::viridis(15),
          xlab = "Node Rank", ylab = "Degree",
          main = "Top 15 Hub Nodes")
  
  n_edges <- sum(network_matrix > 0)
  n_nodes <- nrow(network_matrix)
  density <- n_edges / (n_nodes * (n_nodes - 1))
  
  plot(n_nodes, n_edges, pch = 19, cex = 2, col = "steelblue",
       xlab = "Node Count", ylab = "Edge Count",
       main = "Network Scale")
  text(n_nodes * 0.7, n_edges * 0.8, 
       paste("Density:", round(density, 4)),
       cex = 1.2)
  
  dev.off()
  cat("Network statistics plot saved:", paste0(output_dir, "/network_statistics_", sample_id, ".pdf"), "\n")
}

# ==============================================================================
# 7. Spatial Expression Heatmap (空间热力图 - 基因表达随时间变化)
# ==============================================================================
#
# 【图表目的】
# 展示高变基因沿伪时间轴的表达变化模式，
# 揭示基因表达的时序动态。
#
# 【如何解读】
#   - 行：基因（按聚类排序）
#   - 列：时间箱（T1-T20）
#   - 颜色：平均表达水平
#     - 蓝色：低表达
#     - 白色：中等表达
#     - 红色：高表达
#   - 聚类：表达模式相似的基因聚在一起
#
# 【关键信息】
#   - 共表达模块：一起上调或下调的基因群
#   - 时序特异性：在特定时间高表达的基因
#   - 表达趋势：基因表达随时间单调增加或减少
#   - 批次效应：某些基因在所有时间点都高表达
#
# 【与ARNI网络的关联】
#   - 聚类结果可能与网络模块对应
#   - 时序变化反映潜在的因果关系
#   - Hub基因的表达模式通常较为显著
#
# 【分析方法】
#   - 识别早期响应基因：在T1-T5高表达的基因
#   - 识别晚期响应基因：在T15-T20高表达的基因
#   - 识别持续表达基因：在所有时间点都高表达的基因
#   - 识别瞬时表达基因：在特定时间窗口高表达的基因
#
plot_expression_heatmap <- function(seurat_obj, results, network_matrix, sample_id,
                                     output_dir = "ARNI_visualization") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  pdf(file = paste0(output_dir, "/expression_heatmap_", sample_id, ".pdf"), width = 12, height = 10)
  
  subset_idx <- seurat_obj$orig.ident == sample_id
  
  if ("VariableFeatures" %in% names(seurat_obj@assays$Spatial)) {
    hvg <- seurat_obj@assays$Spatial@VariableFeatures[1:20]
  } else {
    gene_names_all <- rownames(seurat_obj@assays$Spatial$scale.data)
    top_genes_idx <- head(order(colSums(network_matrix), decreasing = TRUE), 20)
    top_genes_idx <- top_genes_idx[top_genes_idx <= length(gene_names_all)]
    hvg <- gene_names_all[top_genes_idx]
  }
  
  if (length(hvg) == 0) {
    cat("  ⚠ No valid genes found, skipping heatmap\n")
    dev.off()
    return()
  }
  
  if ("scale.data" %in% SeuratObject::Layers(seurat_obj@assays$Spatial)) {
    scale_data <- seurat_obj@assays$Spatial$scale.data
  } else {
    scale_data <- seurat_obj@assays$Spatial@layers$scale.data
  }
  
  coords_y <- seurat_obj@meta.data[subset_idx, "y"]
  n_bins <- 20
  y_quantiles <- seq(0, 1, length.out = n_bins + 1)
  y_breaks <- quantile(coords_y, probs = y_quantiles)
  y_breaks[1] <- y_breaks[1] - 0.1
  y_breaks[length(y_breaks)] <- y_breaks[length(y_breaks)] + 0.1
  
  time_bins <- cut(coords_y, breaks = y_breaks, labels = FALSE, include.lowest = TRUE)
  
  gene_names <- rownames(scale_data)
  gene_idx <- match(hvg, gene_names)
  gene_idx <- gene_idx[!is.na(gene_idx)]
  
  if (length(gene_idx) == 0) {
    cat("  ⚠ No matching genes found in scale.data, skipping heatmap\n")
    dev.off()
    return()
  }
  
  expr_by_time <- matrix(0, nrow = length(gene_idx), ncol = n_bins)
  rownames(expr_by_time) <- gene_names[gene_idx]
  colnames(expr_by_time) <- paste0("T", 1:n_bins)
  
  for (t in 1:n_bins) {
    cell_idx <- which(time_bins == t)
    if (length(cell_idx) > 0) {
      subset_data <- scale_data[gene_idx, subset_idx, drop = FALSE]
      expr_by_time[, t] <- rowMeans(subset_data[, cell_idx, drop = FALSE])
    }
  }
  
  if (nrow(expr_by_time) > 0 && ncol(expr_by_time) > 0) {
    pheatmap(expr_by_time,
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             color = colorRampPalette(c("blue", "white", "red"))(100),
             xlab = "Pseudo-time",
             ylab = "Gene",
             main = paste("Gene Expression Over Time -", sample_id),
             fontsize = 8)
  } else {
    cat("  ⚠ Expression matrix is empty, skipping heatmap\n")
  }
  
  dev.off()
  cat("Expression heatmap saved:", paste0(output_dir, "/expression_heatmap_", sample_id, ".pdf"), "\n")
}

# ==============================================================================
# Main Visualization Function (主可视化函数)
# ==============================================================================
#
# 【功能说明】
# 整合所有图表演示，批量处理所有样本的可视化。
# 自动调用上述所有绘图函数，生成完整的分析报告。
#
generate_visualization_report <- function(seurat_obj, all_results, 
                                          output_dir = "ARNI_visualization") {
  
  cat("\n========================================\n")
  cat("ARNI Results Visualization\n")
  cat("========================================\n")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  for (sample_id in names(all_results)) {
    cat("\n--- Visualizing sample:", sample_id, "---\n")
    
    result <- all_results[[sample_id]]
    
    network_file <- paste0("ARNI_results/network_matrix_", sample_id, ".rds")
    if (file.exists(network_file)) {
      network_matrix <- readRDS(network_file)
      
      cat("  Generating learning curves...\n")
      tryCatch({
        plot_learning_curves(result$results, sample_id, output_dir)
      }, error = function(e) {
        cat("  ⚠ Learning curves generation failed:", e$message, "\n")
      })
      
      cat("  Generating network heatmap...\n")
      tryCatch({
        plot_network_heatmap(network_matrix, sample_id, output_dir)
      }, error = function(e) {
        cat("  ⚠ Network heatmap generation failed:", e$message, "\n")
      })
      
      cat("  Generating network graph...\n")
      tryCatch({
        plot_network_graph(network_matrix, sample_id, seurat_obj, output_dir)
      }, error = function(e) {
        cat("  ⚠ Network graph generation failed:", e$message, "\n")
      })
      
      cat("  Generating spatial network visualization...\n")
      tryCatch({
        plot_spatial_with_network(seurat_obj, result$results, network_matrix, sample_id, output_dir)
      }, error = function(e) {
        cat("  ⚠ Spatial network visualization generation failed:", e$message, "\n")
        cat("     Trying simplified version...\n")
        tryCatch({
          plot_network_statistics(result$results, network_matrix, sample_id, output_dir)
        }, error = function(e2) {
          cat("  ⚠ Simplified version also failed:", e2$message, "\n")
        })
      })
      
      cat("  Generating network statistics...\n")
      tryCatch({
        plot_network_statistics(result$results, network_matrix, sample_id, output_dir)
      }, error = function(e) {
        cat("  ⚠ Network statistics generation failed:", e$message, "\n")
      })
      
      cat("  Generating expression heatmap...\n")
      tryCatch({
        plot_expression_heatmap(seurat_obj, result$results, network_matrix, sample_id, output_dir)
      }, error = function(e) {
        cat("  ⚠ Expression heatmap generation failed:", e$message, "\n")
      })
      
      cat("  Preparing Seurat visualization code...\n")
      tryCatch({
        prepare_seurat_plots(seurat_obj, result$results, network_matrix, sample_id, output_dir)
      }, error = function(e) {
        cat("  ⚠ Seurat visualization code generation failed:", e$message, "\n")
      })
      
    } else {
      cat("  Warning: Network matrix file not found", network_file, "\n")
    }
  }
  
  cat("\n========================================\n")
  cat("Visualization Complete!\n")
  cat("========================================\n")
  cat("All plots saved to:", output_dir, "/\n")
  cat("Recommended to run the following scripts for high-quality Seurat spatial plots:\n")
  for (sample_id in names(all_results)) {
    cat("  source('", output_dir, "/seurat_spatial_plots_", sample_id, ".R')\n")
  }
}

# ==============================================================================
# Main Program (主程序)
# ==============================================================================
main <- function() {
  
  cat("Loading Seurat object...\n")
  skin1 <- readRDS("skin1.rds")
  
  cat("Loading ARNI analysis results...\n")
  all_results <- readRDS("ARNI_results/all_results.rds")
  
  cat("Found samples:", names(all_results), "\n")
  
  generate_visualization_report(skin1, all_results, output_dir = "ARNI_visualization")
}

main()