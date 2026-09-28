# ==============================================================================
# 06_train_crop_model.R
# 多时相/多光谱作物识别模型专用训练与验证脚本
# ==============================================================================

# 1. 加载用户 R 库与依赖包
r_ver <- sprintf("%s.%s", R.version$major, substr(R.version$minor, 1, 1))
user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R", "win-library", r_ver)
if (dir.exists(user_lib)) {
  .libPaths(c(user_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(randomForest)
  library(tibble)
  library(dplyr)
})

message("=================================================================")
message(">>> 遥感农作物识别 - 机器学习模型专用训练系统")
message("=================================================================")

output_dir <- if (dir.exists("datacube_crop_classification")) {
  "datacube_crop_classification/models"
} else {
  "models"
}
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 2. 检查是否有用户样点文件
sample_candidates <- c(
  "datacube_crop_classification/data/my_crop_samples.shp",
  "data/my_crop_samples.shp",
  "datacube_crop_classification/data/my_crop_samples.csv",
  "data/my_crop_samples.csv"
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
    files <- files[!grepl("_(NDVI|EVI|Cropland_Mask|CropTypeMap)\\.tif$", files)]
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

# 4. 获取样点数据
if (!is.null(user_sample_file)) {
  message(">>> 成功加载地面真实训练样点: ", user_sample_file)
  if (grepl("\\.csv$", user_sample_file, ignore.case = TRUE)) {
    df <- read.csv(user_sample_file)
    samples_sf <- sf::st_as_sf(df, coords = c("lon", "lat"), crs = 4326)
  } else {
    samples_sf <- sf::st_read(user_sample_file, quiet = TRUE)
  }
} else {
  message("-----------------------------------------------------------------")
  message("【提示】: 当前未在 data/ 下检测到真实样点文件 (如 my_crop_samples.shp/csv)。")
  message(">>> 正在基于影像多光谱与物候特征，自动构建高纯度地表覆盖训练集...")
  message("-----------------------------------------------------------------")
  
  red <- cube_raster[["Red"]]
  nir <- cube_raster[["NIR"]]
  swir1 <- cube_raster[["SWIR1"]]
  ndvi <- (nir - red) / (nir + red)
  ndwi <- (cube_raster[["Green"]] - nir) / (cube_raster[["Green"]] + nir)
  
  set.seed(42)
  samples_list <- list()
  
  # 1) 水体 (NDWI 较高或 NDVI < 0)
  mask_water <- (ndwi > 0 | ndvi < 0)
  pts_water <- terra::spatSample(mask_water, size = 300, as.points = TRUE, na.rm = TRUE)
  pts_water <- pts_water[pts_water[[1]] == 1, ]
  if (length(pts_water) > 10) {
    pts_water$label <- "Water"
    samples_list[[length(samples_list) + 1]] <- pts_water
  }
  
  # 2) 农田作物高覆盖区 (0.45 <= NDVI <= 0.75)
  mask_crop <- (ndvi >= 0.45 & ndvi <= 0.75)
  pts_crop <- terra::spatSample(mask_crop, size = 400, as.points = TRUE, na.rm = TRUE)
  pts_crop <- pts_crop[pts_crop[[1]] == 1, ]
  if (length(pts_crop) > 10) {
    pts_crop$label <- "Cropland_Crops"
    samples_list[[length(samples_list) + 1]] <- pts_crop
  }
  
  # 3) 密林/常绿树木 (NDVI > 0.78)
  mask_forest <- (ndvi > 0.78)
  pts_forest <- terra::spatSample(mask_forest, size = 300, as.points = TRUE, na.rm = TRUE)
  pts_forest <- pts_forest[pts_forest[[1]] == 1, ]
  if (length(pts_forest) > 10) {
    pts_forest$label <- "Forest"
    samples_list[[length(samples_list) + 1]] <- pts_forest
  }
  
  # 4) 裸土/建设用地 (0 <= NDVI < 0.25 且短波红外较高)
  mask_bare <- (ndvi >= 0 & ndvi < 0.25 & swir1 > 1500)
  pts_bare <- terra::spatSample(mask_bare, size = 300, as.points = TRUE, na.rm = TRUE)
  pts_bare <- pts_bare[pts_bare[[1]] == 1, ]
  if (length(pts_bare) > 10) {
    pts_bare$label <- "Builtup_BareSoil"
    samples_list[[length(samples_list) + 1]] <- pts_bare
  }
  
  # 合并所有样本
  samples_vect <- do.call(rbind, samples_list)
  samples_sf <- sf::st_as_sf(samples_vect)
}

# 5. 提取多光谱反射率特征
message(">>> 正在从数据立方体提取像元光谱多维特征...")
sample_features <- terra::extract(cube_raster, terra::vect(samples_sf), df = TRUE)
sample_features$label <- as.factor(samples_sf$label)
sample_features <- na.omit(sample_features)

message(">>> 样本各类别数量分布:")
print(table(sample_features$label))

# 6. 划分训练集 (80%) 与独立测试集 (20%)
set.seed(123)
train_indices <- sample(1:nrow(sample_features), size = 0.8 * nrow(sample_features))
train_data <- sample_features[train_indices, ]
test_data  <- sample_features[-train_indices, ]

# 7. 训练随机森林多光谱作物分类器
message("\n>>> 开始训练随机森林 (Random Forest) 作物分类模型...")
features_formula <- as.formula("label ~ Blue + Green + Red + NIR + SWIR1 + SWIR2")
rf_model <- randomForest::randomForest(
  features_formula,
  data = train_data,
  ntree = 150,
  importance = TRUE
)

message(">>> 模型训练成功！各波段重要性评分 (Importance):")
print(round(randomForest::importance(rf_model), 2))

# 8. 独立测试集精度评估 (Confusion Matrix & Accuracy)
message("\n=================================================================")
message(">>> 独立测试集 (Test Set) 精度验证评估报告")
message("=================================================================")
test_pred <- predict(rf_model, newdata = test_data)
conf_matrix <- table(真实标签 = test_data$label, 预测标签 = test_pred)
print(conf_matrix)

overall_acc <- sum(diag(conf_matrix)) / sum(conf_matrix)
message("\n>>> 【总体分类精度 (Overall Accuracy)】: ", sprintf("%.2f%%", overall_acc * 100))

# 9. 保存训练好的模型
model_save_path <- file.path(output_dir, "crop_rf_model.rds")
saveRDS(rf_model, model_save_path)
message("\n>>> 训练好的模型已保存至: ", model_save_path)

# 10. 联动执行：对该影像执行分类并输出真正的分类专题图
message("\n>>> 正在应用训练好的模型对整幅数据立方体进行分类制图...")
res_dir <- if (dir.exists("datacube_crop_classification")) {
  "datacube_crop_classification/output_sdc30"
} else {
  "output_sdc30"
}
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

crop_map <- terra::predict(cube_raster, rf_model)

# 显式重映射为标准二值编码：
# 1 = 耕地作物 (Cropland_Crops)
# 0 = 非耕地 (Builtup_BareSoil / 裸土 / 建筑等)
# 255 = 空值像元 (NoData)
message(">>> 正在将模型分类结果规范化为标准二值编码 (1=耕地作物, 0=非耕地)...")
crop_map_binary <- terra::ifel(
  is.na(crop_map),
  255,
  terra::ifel(crop_map == "Cropland_Crops", 1, 0)
)

crop_out_path <- file.path(res_dir, paste0(tools::file_path_sans_ext(basename(target_tif)), "_CropTypeMap.tif"))
dir.create(dirname(crop_out_path), recursive = TRUE, showWarnings = FALSE)
terra::writeRaster(crop_map_binary, crop_out_path, datatype = "INT1U", NAflag = 255, overwrite = TRUE)

message(">>> 恭喜！标准二值作物分类图已生成:")
message("    -> 本地文件: ", crop_out_path)

# 自动同步更新到 E:\agriculture\gaced30_validation_pipeline (如存在)
val_source_dir <- "E:/agriculture/gaced30_validation_pipeline/data/source_data"
if (dir.exists(val_source_dir)) {
  val_target_file <- file.path(val_source_dir, basename(crop_out_path))
  file.copy(crop_out_path, val_target_file, overwrite = TRUE)
  message(">>> 🚀 [自动同步] 已将最新修正编码的分类图同步至验证管线: ", val_target_file)
}
message("=================================================================")
