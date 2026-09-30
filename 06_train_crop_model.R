# ==============================================================================
# 06_train_crop_model.R
# 联合国粮农组织《遥感农业统计手册》- 多时相/多光谱作物识别模型高精度训练与制图系统 (优化升级版)
# ==============================================================================
# 核心提升：
# 1. 遥感特征工程：融合 6 原始波段与 NDVI、EVI、MNDWI、LSWI、NDBI、NDTI 6 大核心指数；
# 2. 高性能多线程：采用 C++ 级 ranger 算法，训练提速 10~20 倍，内存占用减半；
# 3. 联合国标准精度台账：输出 OA、PA(查全率)、UA(查准率)、F1-Score 与 Kappa 系数；
# 4. 空间上下文后处理：应用 3x3 空间众数滤波平滑，消除椒盐斑点，提升田块连通性与交并比；
# 5. 四宫格质检可视化：自动生成【假彩色/真值/预测/空间误差差分】高分辨率诊断成果图。
# ==============================================================================

# 1. 环境与依赖库加载
r_ver <- sprintf("%s.%s", R.version$major, substr(R.version$minor, 1, 1))
user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R", "win-library", r_ver)
if (dir.exists(user_lib)) {
  .libPaths(c(user_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(tibble)
  library(dplyr)
})

# 动态加载 ranger (优先) 或 randomForest (兜底)
use_ranger <- requireNamespace("ranger", quietly = TRUE)
if (use_ranger) {
  suppressPackageStartupMessages(library(ranger))
  message(">>> [引擎] 已启用高性能多线程 ranger 随机森林引擎。")
} else {
  suppressPackageStartupMessages(library(randomForest))
  message(">>> [引擎] 使用原生 randomForest 引擎。")
}

message("=================================================================")
message(">>> 遥感农作物识别 - 机器学习模型专用训练与制图系统 (优化升级版)")
message("=================================================================")

output_dir <- if (dir.exists("datacube_crop_classification")) {
  "datacube_crop_classification/models"
} else {
  "models"
}
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 2. 检查并定位样本文件
sample_candidates <- c(
  "datacube_crop_classification/data/my_crop_samples.csv",
  "data/my_crop_samples.csv",
  "datacube_crop_classification/data/my_crop_samples.shp",
  "data/my_crop_samples.shp"
)

user_sample_file <- NULL
for (f in sample_candidates) {
  if (file.exists(f)) {
    user_sample_file <- f
    break
  }
}

# 3. 寻找影像数据立方体
candidate_dirs <- c(
  "datacube_crop_classification/data/sdc30_cubes",
  "data/sdc30_cubes",
  "data"
)
tif_files <- character(0)
for (d in candidate_dirs) {
  if (dir.exists(d)) {
    files <- list.files(d, pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE)
    files <- files[!grepl("_(NDVI|EVI|Cropland_Mask|CropTypeMap|Diagnostic)\\.tif$", files)]
    if (length(files) > 0) {
      tif_files <- files
      break
    }
  }
}

if (length(tif_files) == 0) {
  stop("未在 data/sdc30_cubes 下找到数据立方体 TIF 影像！")
}

target_tif <- tif_files[1]
message(">>> 使用的数据立方体: ", basename(target_tif))
cube_raster <- terra::rast(target_tif)
if (nlyr(cube_raster) >= 6) {
  names(cube_raster)[1:6] <- c("Blue", "Green", "Red", "NIR", "SWIR1", "SWIR2")
}

feature_cols <- c("Blue", "Green", "Red", "NIR", "SWIR1", "SWIR2",
                  "NDVI", "EVI", "MNDWI", "LSWI", "NDBI", "NDTI")

# 4. 获取样点数据并提取特征
if (!is.null(user_sample_file)) {
  message(">>> 成功加载地面真实训练样点: ", user_sample_file)
  if (grepl("\\.csv$", user_sample_file, ignore.case = TRUE)) {
    df_raw <- read.csv(user_sample_file)
    # 若样本超 8 万条，分层均衡抽样以控制训练耗时并保证高表征力
    if (nrow(df_raw) > 80000) {
      set.seed(42)
      df_samples <- df_raw %>%
        group_by(label) %>%
        slice_sample(n = 40000) %>%
        ungroup()
      message(">>> 样点总数较丰富 (", nrow(df_raw), " 点)，执行最优分层抽样: ", nrow(df_samples), " 点。")
    } else {
      df_samples <- df_raw
    }

    # 检查 CSV 是否已自带直出的 12 维特征
    has_features <- all(feature_cols %in% names(df_samples))
    if (has_features) {
      message(">>> [极速直出通道] 检测到样本文件已包含完整 12 维遥感多光谱波段与植被指数！")
      message(">>> 跳过耗时的空间栅格提取，直通极速训练...")
      sample_features <- df_samples
      sample_features$label <- as.factor(sample_features$label)
      sample_features <- na.omit(sample_features)
    } else {
      samples_vect <- terra::vect(df_samples, geom = c("lon", "lat"), crs = "EPSG:4326")
    }
  } else {
    samples_sf <- sf::st_read(user_sample_file, quiet = TRUE)
    samples_vect <- terra::vect(samples_sf)
  }
} else {
  stop("未检测到真实样点文件，请先运行 prepare_samples.py 提取样点！")
}

# 5. 若未自带特征，则从数据立方体提取像元多光谱特征并计算物候指数
if (!exists("sample_features")) {
  message(">>> 正在从数据立方体提取像元多光谱特征...")
  sample_features <- terra::extract(cube_raster, samples_vect, df = TRUE)
  sample_features$label <- as.factor(samples_vect$label)
  sample_features <- na.omit(sample_features)

  message(">>> 正在构建遥感作物识别核心光谱与物候指数 (NDVI, EVI, MNDWI, LSWI, NDBI, NDTI)...")
  sample_features <- sample_features %>%
    mutate(
      b = Blue / 10000,
      g = Green / 10000,
      r = Red / 10000,
      nir = NIR / 10000,
      s1 = SWIR1 / 10000,
      s2 = SWIR2 / 10000,
      # 植被生长与生物量指数
      NDVI  = (nir - r) / (nir + r + 1e-6),
      EVI   = 2.5 * (nir - r) / (nir + 6 * r - 7.5 * b + 1.0 + 1e-6),
      # 水分与水体指数 (区分沟渠、水田与洼地)
      MNDWI = (g - s1) / (g + s1 + 1e-6),
      LSWI  = (nir - s1) / (nir + s1 + 1e-6),
      # 不透水面与建筑指数 (区分城镇、道路与村庄)
      NDBI  = (s1 - nir) / (s1 + nir + 1e-6),
      # 耕作/麦茬/秸秆残留指数 (区分收割麦田与常年裸地)
      NDTI  = (s1 - s2) / (s1 + s2 + 1e-6)
    )
}

message(">>> 样本各类别数量分布:")
print(table(sample_features$label))

# 6. 划分训练集 (80%) 与独立测试集 (20%)
set.seed(123)
train_indices <- sample(1:nrow(sample_features), size = 0.8 * nrow(sample_features))
train_data <- sample_features[train_indices, ]
test_data  <- sample_features[-train_indices, ]

# 7. 训练高维特征随机森林分类器
message("\n>>> 开始训练多维特征随机森林作物分类模型...")
features_formula <- as.formula(paste("label ~", paste(feature_cols, collapse = " + ")))

t_start <- Sys.time()
if (use_ranger) {
  rf_model <- ranger::ranger(
    formula = features_formula,
    data = train_data,
    num.trees = 150,
    importance = "impurity",
    probability = FALSE,
    num.threads = max(1, parallel::detectCores() - 1)
  )
  train_duration <- round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1)
  message(">>> ranger 模型训练成功！耗时: ", train_duration, " 秒 (OOB 误差: ", round(rf_model$prediction.error * 100, 2), "%)")
  message(">>> 各特征重要性评分 (Importance Top 6):")
  imp_sorted <- sort(ranger::importance(rf_model), decreasing = TRUE)
  print(round(head(imp_sorted, 6), 2))
} else {
  rf_model <- randomForest::randomForest(
    features_formula,
    data = train_data,
    ntree = 150,
    importance = TRUE
  )
  train_duration <- round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1)
  message(">>> randomForest 模型训练成功！耗时: ", train_duration, " 秒")
}

# 8. 独立测试集精度评估 (Confusion Matrix & Accuracy)
message("\n=================================================================")
message(">>> 独立测试集 (Test Set) 精度验证评估报告")
message("=================================================================")
if (use_ranger) {
  test_pred <- predict(rf_model, data = test_data)$predictions
} else {
  test_pred <- predict(rf_model, newdata = test_data)
}

conf_matrix <- table(真实标签 = test_data$label, 预测标签 = test_pred)
print(conf_matrix)

# 联合国精度指标计算
total_n <- sum(conf_matrix)
oa <- sum(diag(conf_matrix)) / total_n

# 针对耕地类 (Cropland)
crop_label_name <- rownames(conf_matrix)[grepl("^Cropland", rownames(conf_matrix))][1]
pa_crop <- conf_matrix[crop_label_name, crop_label_name] / sum(conf_matrix[crop_label_name, ])
ua_crop <- conf_matrix[crop_label_name, crop_label_name] / sum(conf_matrix[, crop_label_name])
f1_crop <- 2 * (pa_crop * ua_crop) / (pa_crop + ua_crop)

# Cohen's Kappa 计算
row_sums <- rowSums(conf_matrix)
col_sums <- colSums(conf_matrix)
pe <- sum(row_sums * col_sums) / (total_n^2)
kappa <- (oa - pe) / (1 - pe)

message("\n-----------------------------------------------------------------")
message(sprintf(">>> 【总体分类精度 (Overall Accuracy)】  : %.2f%%", oa * 100))
message(sprintf(">>> 【耕地查全率 (Producer's Acc / PA)】: %.2f%%", pa_crop * 100))
message(sprintf(">>> 【耕地查准率 (User's Acc / UA)】    : %.2f%%", ua_crop * 100))
message(sprintf(">>> 【耕地综合质量 (F1-Score)】         : %.2f%%", f1_crop * 100))
message(sprintf(">>> 【Kappa 一致性系数】                : %.4f", kappa))
message("-----------------------------------------------------------------")

# 9. 保存训练好的模型
model_save_path <- file.path(output_dir, "crop_rf_model.rds")
saveRDS(rf_model, model_save_path)
message(">>> 训练好的模型已保存至: ", model_save_path)

# 10. 全景栅格推理预测与空间上下文平滑
message("\n>>> 正在准备多维全景栅格特征堆叠 (计算全幅 NDVI, EVI, MNDWI, LSWI, NDBI, NDTI)...")
b_r   <- cube_raster[["Blue"]] / 10000
g_r   <- cube_raster[["Green"]] / 10000
r_r   <- cube_raster[["Red"]] / 10000
nir_r <- cube_raster[["NIR"]] / 10000
s1_r  <- cube_raster[["SWIR1"]] / 10000
s2_r  <- cube_raster[["SWIR2"]] / 10000

ndvi_r  <- (nir_r - r_r) / (nir_r + r_r + 1e-6); names(ndvi_r) <- "NDVI"
evi_r   <- 2.5 * (nir_r - r_r) / (nir_r + 6 * r_r - 7.5 * b_r + 1.0 + 1e-6); names(evi_r) <- "EVI"
mndwi_r <- (g_r - s1_r) / (g_r + s1_r + 1e-6); names(mndwi_r) <- "MNDWI"
lswi_r  <- (nir_r - s1_r) / (nir_r + s1_r + 1e-6); names(lswi_r) <- "LSWI"
ndbi_r  <- (s1_r - nir_r) / (s1_r + nir_r + 1e-6); names(ndbi_r) <- "NDBI"
ndti_r  <- (s1_r - s2_r) / (s1_r + s2_r + 1e-6); names(ndti_r) <- "NDTI"

full_feature_stack <- c(cube_raster, ndvi_r, evi_r, mndwi_r, lswi_r, ndbi_r, ndti_r)

message(">>> 正在应用训练好的模型对整幅数据立方体进行并行预测制图...")
t_pred_start <- Sys.time()
if (use_ranger) {
  pred_wrapper <- function(model, data, ...) {
    preds <- ranger:::predict.ranger(model, data = as.data.frame(data))$predictions
    as.integer(preds == crop_label_name)
  }
  crop_map_raw <- terra::predict(full_feature_stack, rf_model, fun = pred_wrapper)
} else {
  crop_map_pred <- terra::predict(full_feature_stack, rf_model)
  crop_map_raw <- terra::ifel(crop_map_pred == crop_label_name, 1, 0)
}
pred_duration <- round(as.numeric(difftime(Sys.time(), t_pred_start, units = "secs")), 1)
message(">>> 全像素预测完成！耗时: ", pred_duration, " 秒。")

# 11. 空间上下文滤波平滑 (去除椒盐斑点噪声，强化地块完整性)
message(">>> 正在执行 3x3 空间众数滤波平滑 (消除椒盐斑点，提升田块连续性)...")
crop_map_smooth <- terra::focal(crop_map_raw, w = 3, fun = "modal", na.policy = "omit")

# 规范化编码 (1=耕地, 0=非耕地, 255=NoData)
crop_map_final <- terra::ifel(is.na(crop_map_smooth), 255, crop_map_smooth)

res_dir <- if (dir.exists("datacube_crop_classification")) {
  "datacube_crop_classification/output_sdc30"
} else {
  "output_sdc30"
}
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
crop_out_path <- file.path(res_dir, paste0(tools::file_path_sans_ext(basename(target_tif)), "_CropTypeMap.tif"))
terra::writeRaster(crop_map_final, crop_out_path, datatype = "INT1U", NAflag = 255, overwrite = TRUE)

message(">>> 恭喜！空间平滑二值作物分类图已生成:")
message("    -> 本地文件: ", crop_out_path)

# 12. 自动生成四宫格高分辨率质检诊断成果图
message("\n>>> 正在生成高分辨率四宫格质检对比图 (假彩色/真值/分类/空间误差差分)...")
gt_file_candidate <- file.path(dirname(target_tif), "..", paste0(tools::file_path_sans_ext(basename(target_tif)), "_GroundTruth.tif"))
if (!file.exists(gt_file_candidate)) {
  gt_file_candidate <- "datacube_crop_classification/data/SDC30_V003_50SMF_20210618_GroundTruth.tif"
}

diag_png_path <- file.path(res_dir, "SDC30_Crop_Classification_Diagnostic.png")

if (file.exists(gt_file_candidate)) {
  gt_raster <- terra::rast(gt_file_candidate)
  
  # 计算空间误差差分图:
  # 1 = TP (真阳性: 耕地正确识别, 绿)
  # 2 = TN (真阴性: 非耕地正确识别, 浅灰)
  # 3 = FP (假阳性: 虚报为耕地, 红)
  # 4 = FN (假阴性: 漏报漏识耕地, 橙黄)
  diff_map <- terra::ifel(
    gt_raster == 1 & crop_map_final == 1, 1,
    terra::ifel(gt_raster == 0 & crop_map_final == 0, 2,
    terra::ifel(gt_raster == 0 & crop_map_final == 1, 3,
    terra::ifel(gt_raster == 1 & crop_map_final == 0, 4, 255)))
  )
  
  png(diag_png_path, width = 2400, height = 2400, res = 200)
  par(mfrow = c(2, 2), mar = c(3, 3, 3, 1))
  
  # 1. 标准假彩色合成 (NIR-Red-Green: 植被呈亮红色)
  terra::plotRGB(cube_raster, r = 4, g = 3, b = 2, stretch = "lin",
                 main = "① SDC30 标准假彩色合成 (NIR-Red-Green)")
  
  # 2. 地面参考真值底图
  terra::plot(gt_raster, col = c("#F0F0F0", "#228B22"), legend = FALSE,
              main = "② 地面参考真值底图 (绿色=耕地, 灰白=非耕地)")
  
  # 3. 本次模型预测制图 (空间平滑后)
  terra::plot(crop_map_final, col = c("#F0F0F0", "#228B22"), legend = FALSE,
              main = sprintf("③ 作物识别预测成果图 (平滑后, F1=%.1f%%)", f1_crop * 100))
  
  # 4. 空间误差差分图
  terra::plot(diff_map, col = c("#228B22", "#E8E8E8", "#FF3030", "#FFA500"), legend = FALSE,
              main = "④ 空间误差差分图 (绿=TP对, 灰=TN对, 红=FP虚报, 橙=FN漏报)")
  
  dev.off()
  message(">>> 四宫格质检诊断图已成功生成: ", diag_png_path)
}

message("=================================================================")
message(">>> 全部训练、制图与质检流程执行完毕！")
message("=================================================================")
