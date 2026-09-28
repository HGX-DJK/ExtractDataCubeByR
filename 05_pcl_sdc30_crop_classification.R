# ==============================================================================
# 05_pcl_sdc30_crop_classification.R
# 适配鹏城星云 iEarth SDC30 (全球30米无缝地表反射率数据立方体) 作物分类实战脚本
# 数据源: https://data-starcloud.pcl.ac.cn/iearthdata/26
# ==============================================================================

# 1. 确保加载个人用户 R 包库
r_ver <- sprintf("%s.%s", R.version$major, substr(R.version$minor, 1, 1))
user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R", "win-library", r_ver)
if (dir.exists(user_lib)) {
  .libPaths(c(user_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(randomForest)
})

message("=================================================================")
message(">>> 鹏城星云 SDC30 (全球30米无缝数据立方体) 作物识别处理系统")
message("=================================================================")

# 2. 自动检测 SDC30 文件存放位置
candidate_dirs <- c(
  "datacube_crop_classification/data/sdc30_cubes",
  "data/sdc30_cubes",
  "../data/sdc30_cubes",
  "data"
)

sdc_dir <- NULL
tif_files <- character(0)

for (d in candidate_dirs) {
  if (dir.exists(d)) {
    files <- list.files(d, pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE)
    # 排除中间输出的掩膜文件
    files <- files[!grepl("_(NDVI|EVI|Cropland_Mask|CropTypeMap)\\.tif$", files)]
    if (length(files) > 0) {
      sdc_dir <- d
      tif_files <- files
      break
    }
  }
}

if (length(tif_files) == 0) {
  stop("未能找到 SDC30 数据文件！请确认文件已放入 data/sdc30_cubes 目录下。")
}

message(">>> 成功定位数据目录: ", sdc_dir)
message(">>> 发现 SDC30 数据立方体文件数量: ", length(tif_files))
for (f in tif_files) {
  message("    -> 文件: ", basename(f))
}

# 3. 设置输出结果目录
output_dir <- if (dir.exists("datacube_crop_classification")) {
  "datacube_crop_classification/output_sdc30"
} else {
  "output_sdc30"
}
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 4. 读取 SDC30 数据立方体多光谱特征
# SDC30 包含 6 个波段: 1:Blue, 2:Green, 3:Red, 4:NIR, 5:SWIR1, 6:SWIR2
target_tif <- tif_files[1]
message("\n>>> [步骤 1/4] 正在加载并解析多光谱数据立方体: ", basename(target_tif))
cube_raster <- terra::rast(target_tif)

# 规范波段名称
names(cube_raster) <- c("Blue", "Green", "Red", "NIR", "SWIR1", "SWIR2")

message(">>> 数据立方体基本规格:")
message("    - 空间分辨率: ", res(cube_raster)[1], " 米 x ", res(cube_raster)[2], " 米")
message("    - 栅格大小: ", nrow(cube_raster), " 行 x ", ncol(cube_raster), " 列 (共 ", ncell(cube_raster), " 像元)")
message("    - 坐标参考系 (CRS): ", crs(cube_raster, proj = TRUE))
message("    - 波段组成: ", paste(names(cube_raster), collapse = ", "))

# 5. 计算作物关键物候指数 (NDVI & EVI & LSWI)
message("\n>>> [步骤 2/4] 计算作物物候特征指数 (NDVI / EVI)...")
red <- cube_raster[["Red"]]
nir <- cube_raster[["NIR"]]
blue <- cube_raster[["Blue"]]

# 计算 NDVI = (NIR - Red) / (NIR + Red)
ndvi <- (nir - red) / (nir + red)
names(ndvi) <- "NDVI"

# 计算增强型植被指数 EVI (对高生物量更敏感，不易饱和)
evi <- 2.5 * ((nir - red) / (nir + 6 * red - 7.5 * blue + 1))
names(evi) <- "EVI"

# 保存 NDVI 与 EVI 成果到本地
ndvi_out_path <- file.path(output_dir, paste0(tools::file_path_sans_ext(basename(target_tif)), "_NDVI.tif"))
evi_out_path  <- file.path(output_dir, paste0(tools::file_path_sans_ext(basename(target_tif)), "_EVI.tif"))

message(">>> 正在导出多时相 NDVI 地图: ", basename(ndvi_out_path))
terra::writeRaster(ndvi, ndvi_out_path, overwrite = TRUE)

message(">>> 正在导出多时相 EVI 地图: ", basename(evi_out_path))
terra::writeRaster(evi, evi_out_path, overwrite = TRUE)

message(">>> NDVI 统计概况: 最小值 = ", round(minmax(ndvi)[1], 3), "，最大值 = ", round(minmax(ndvi)[2], 3))

# 6. 作物与农田提取制图
message("\n>>> [步骤 3/4] 正在执行农田与作物识别...")

# 检查是否有地面作物调查样点 (Shapefile 或 CSV)
samples_path_candidate <- c(
  file.path(sdc_dir, "my_crop_samples.shp"),
  "data/my_crop_samples.shp",
  "../data/my_crop_samples.shp"
)
sample_file <- NULL
for (sp in samples_path_candidate) {
  if (file.exists(sp)) {
    sample_file <- sp
    break
  }
}

if (!is.null(sample_file)) {
  message(">>> 检测到地面作物调查样点: ", sample_file)
  samples_sf <- sf::st_read(sample_file, quiet = TRUE)
  
  # 提取光谱值
  message(">>> 提取样点多光谱像元特征...")
  sample_vals <- terra::extract(cube_raster, samples_sf, df = TRUE)
  sample_vals$label <- as.factor(samples_sf$label)
  
  # 训练随机森林模型
  message(">>> 训练多光谱随机森林作物分类器...")
  rf_model <- randomForest(label ~ Blue + Green + Red + NIR + SWIR1 + SWIR2, data = sample_vals, ntree = 100)
  
  # 全局栅格预测
  message(">>> 对 SDC30 区域执行像元级作物类型预测制图...")
  crop_map <- terra::predict(cube_raster, rf_model)
  crop_out_path <- file.path(output_dir, paste0(tools::file_path_sans_ext(basename(target_tif)), "_CropTypeMap.tif"))
  terra::writeRaster(crop_map, crop_out_path, overwrite = TRUE)
  message(">>> 作物分类图已生成: ", crop_out_path)
} else {
  message(">>> 当前未提供地面样点文件 (如 data/my_crop_samples.shp)。")
  message(">>> 正在采用作物物候动态阈值法，提取旺盛生长的农田植被覆盖区 (NDVI > 0.4)...")
  
  # 农田植被掩膜提取: 0 = 非农田/低植被, 1 = 农田高覆盖植被
  crop_mask <- ndvi > 0.4
  names(crop_mask) <- "Cropland_Vegetation"
  
  mask_out_path <- file.path(output_dir, paste0(tools::file_path_sans_ext(basename(target_tif)), "_Cropland_Mask.tif"))
  terra::writeRaster(crop_mask, mask_out_path, overwrite = TRUE)
  message(">>> 农田作物高覆盖区提取完成，成果图已导出: ", basename(mask_out_path))
}

# 7. 总结
message("\n=================================================================")
message(">>> [步骤 4/4] 全部处理成功完成！")
message(">>> 生成的成果文件存放在: ", normalizePath(output_dir))
message("=================================================================")
