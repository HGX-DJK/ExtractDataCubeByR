# ==============================================================================
# 04_multicube_data_extraction.R
# 联合国粮农组织《遥感农业统计手册》- 多数据立方体时序数据提取专有工具脚本
# ==============================================================================
# 说明：
# 本脚本系统演示如何使用 sits 框架从多种“数据立方体”中提取农作物时序数据：
# 1. 单立方体点矢量采样提取 (Point sampling)
# 2. 多边形农田地块采样提取 (Polygon sampling, 支持每地块提取 N 个像元)
# 3. 多源数据立方体融合提取 (Multi-source Cube Fusion: 光学 Sentinel-2 + 高程 DEM + 雷达)
# 4. 跨多瓦片 (Multi-tile) 大范围数据立方体并行提取
# 5. 提取结果导出为标准表格 (CSV / RDS / GeoParquet) 供下游 Python/R 算法使用
# ==============================================================================

# 1. 确保工作目录为项目根目录
if (!file.exists("data/ct_chile") && file.exists("../data/ct_chile")) {
  setwd("..")
}

library(sits)
library(sf)
library(tibble)
library(dplyr)

output_dir <- "datacube_crop_classification/extracted_data"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

message("=================================================================")
message(">>> 场景 1: 从本地多时相数据立方体中根据【点位】提取时序数据")
message("=================================================================")

# 读取地面样点 (包含经纬度和作物标签)
sample_points_path <- "data/ct_chile/Samples/ground_data_19HBA_Chile_eng.shp"
if (file.exists(sample_points_path)) {
  points_sf <- sf::st_read(sample_points_path, quiet = TRUE)
  message(">>> 成功加载点位数据，共计 ", nrow(points_sf), " 个采样点。")
} else {
  # 如果没有本地文件，构建几个虚拟点位做演示
  points_sf <- tibble::tibble(
    longitude = c(-71.25, -71.20, -71.15),
    latitude  = c(-35.40, -35.45, -35.50),
    label     = c("Winter_crops", "Spring_crops", "Deciduous_fruit_tree"),
    id        = 1:3
  ) |> sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
}

# 假设已有规整化立方体（或加载本地已存对象）
cube_file <- "data/ct_chile/sent_19HBA_reg.rds"
if (file.exists(cube_file)) {
  cube_s2 <- readRDS(cube_file)
  message(">>> 成功加载 Sentinel-2 多时相数据立方体: ", cube_file)
  
  # 调用 sits_get_data 从数据立方体中提取时序数据
  message(">>> 正在提取多时相波段时序...")
  extracted_pts <- sits_get_data(
    cube       = cube_s2,
    samples    = points_sf,
    label_attr = "label",
    multicores = 2
  )
  
  # 查看提取结果的结构
  message(">>> 提取完成！数据格式为 sits tibble，前 2 行如下:")
  print(head(extracted_pts, 2))
}

message("\n=================================================================")
message(">>> 场景 2: 多源数据立方体融合提取 (Multi-source Cube Extraction)")
message("=================================================================")
# 例如：同时融合【Sentinel-2 光学时序立方体】与【Copernicus DEM 高程立方体】
dem_cube_file <- "data/ct_chile/dem_19HBA_reg.rds"

if (file.exists(cube_file) && file.exists(dem_cube_file)) {
  cube_s2  <- readRDS(cube_file)
  cube_dem <- readRDS(dem_cube_file)
  
  message(">>> 使用 sits_merge() 将光学时序立方体与 DEM 立方体融合成【多源数据立方体】...")
  # sits_merge 将两者按空间几何对齐，DEM 作为额外波段在各时间步复制
  cube_multi <- sits_merge(cube_s2, cube_dem)
  message(">>> 融合后的立方体波段包含: ", paste(sits_bands(cube_multi), collapse = ", "))
  
  # 从多源数据立方体一次性提取所有光谱波段 + 高程信息
  message(">>> 从多源立方体中提取联合时序特征...")
  extracted_multi <- sits_get_data(
    cube       = cube_multi,
    samples    = points_sf[1:min(50, nrow(points_sf)), ], # 提取前50个点作为演示
    label_attr = "label",
    multicores = 2
  )
  message(">>> 联合时序提取成功！单个样本包含的波段:")
  print(colnames(extracted_multi$time_series[[1]]))
}

message("\n=================================================================")
message(">>> 场景 3: 从【农田多边形地块】(Polygon) 中批量采样提取像元")
message("=================================================================")
# 如果样点是农田矢量多边形 (例如地块)，可以通过 n_sam_pol 指定每个地块采样提取的点数
roi_file <- "data/ct_chile/ROI/ROI_19HBA.shp"
if (file.exists(roi_file) && file.exists(cube_file)) {
  poly_sf <- sf::st_read(roi_file, quiet = TRUE)
  poly_sf$label <- "Study_Area"
  
  message(">>> 从多边形地块内采样提取像元 (例如每地块抽取 20 个像素时序)...")
  extracted_poly <- sits_get_data(
    cube       = cube_s2,
    samples    = poly_sf,
    n_sam_pol  = 20,         # 每个地块抽样 20 个点
    label_attr = "label",
    multicores = 2
  )
  message(">>> 地块内部像元提取成功，总计获得 ", nrow(extracted_poly), " 条像元时间序列。")
}

message("\n=================================================================")
message(">>> 场景 4: 跨多个瓦片 (Multi-tile) 大范围数据立方体提取")
message("=================================================================")
# sits_cube 允许传入多个瓦片，例如 tiles = c("19HBA", "19HBB")
# 此时 sits_get_data 会自动判断每个点落在哪个瓦片，并实现跨瓦片无缝并行提取：
#
# multi_tile_cube <- sits_cube(
#   source = "MPC", collection = "SENTINEL-2-L2A",
#   tiles = c("19HBA", "19HBB"),
#   start_date = "2020-05-01", end_date = "2021-05-30", ...
# )
# all_samples <- sits_get_data(cube = multi_tile_cube, samples = national_points, multicores = 8)

message("\n=================================================================")
message(">>> 场景 5: 导出提取结果为标准文件 (供 Python / 机器学习使用)")
message("=================================================================")

# 将提取到的样本保存为本地 RDS
saveRDS(extracted_pts, file.path(output_dir, "extracted_samples.rds"))

# 扁平化展开为常见的多时相表格 CSV (每一行一个点在某日期的波段值，或按列展开)
flat_ts <- extracted_pts |>
  tidyr::unnest(time_series)

csv_path <- file.path(output_dir, "extracted_timeseries_flat.csv")
readr::write_csv(flat_ts, csv_path)

message(">>> 提取结果已导出: ")
message("   1. 完整 R 对象: ", file.path(output_dir, "extracted_samples.rds"))
message("   2. 扁平化 CSV 表格: ", csv_path)
message(">>> 多数据立方体时序提取演示完成！")
