# ==============================================================================
# 03_deafrica_datacube_crop_mapping.R
# 联合国粮农组织《遥感农业统计手册》- 基于 Digital Earth Africa 开放数据立方体的作物识别
# ==============================================================================
# 说明：
# 本脚本提取自 ct_digital_earth_africa.qmd，演示利用开放数据立方体 (Open Data Cube, ODC)
# 的 GeoMAD (时序中值绝对离差复合产品) 及耕地掩膜对作物大类进行识别。
# ==============================================================================

if (!file.exists("data/ct_digital_earth_africa") && file.exists("../data/ct_digital_earth_africa")) {
  setwd("..")
}

library(sits)
library(sf)
library(terra)
library(randomForest)

output_dir <- "datacube_crop_classification/output_deafrica"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 1. 加载卢旺达作物分类样本
sample_rds <- "data/ct_digital_earth_africa/balanced_samples_comb.rds"
if (file.exists(sample_rds)) {
  message(">>> 加载经过类别平衡的作物大类样本数据 (Cereals, Roots, Legumes, Perennials 等)...")
  balanced_samples <- readRDS(sample_rds)
  print(summary(balanced_samples))
} else {
  stop("未找到样本文件: ", sample_rds)
}

# 2. 训练多作物时序分类器
model_rds <- "data/ct_digital_earth_africa/rfor_model_comb.rds"
if (file.exists(model_rds)) {
  message(">>> 加载预训练模型: ", model_rds)
  rfor_model <- readRDS(model_rds)
} else {
  message(">>> 训练随机森林分类模型...")
  rfor_model <- sits_train(balanced_samples, ml_method = sits_rfor())
}

# 3. 加载掩膜后的 Digital Earth Africa 栅格数据立方体
cube_rds <- "data/ct_digital_earth_africa/dea_s2_masked.rds"
if (file.exists(cube_rds)) {
  message(">>> 加载 DE Africa Sentinel-2 GeoMAD 数据立方体...")
  dea_cube <- readRDS(cube_rds)
  print(dea_cube)
  
  # 执行分类识别
  message(">>> 对 DE Africa 数据立方体执行作物识别...")
  crop_probs <- sits_classify(
    data       = dea_cube,
    ml_model   = rfor_model,
    output_dir = output_dir,
    version    = "rf_comb",
    multicores = 1
  )
  message(">>> DE Africa 作物分类预测完成！")
} else {
  message(">>> 本地未找到 dea_s2_masked.rds，已完成模型构建与样本评估。")
}
