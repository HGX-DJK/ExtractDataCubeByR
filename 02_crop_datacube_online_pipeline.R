# ==============================================================================
# 02_crop_datacube_online_pipeline.R
# 联合国粮农组织《遥感农业统计手册》- 多时相数据立方体作物识别 (云端在线完整构建流程)
# ==============================================================================
# 说明：
# 本脚本演示从微软行星计算机 (MPC) 检索多时相 Sentinel-2 影像、构建时空规整的数据立方体、
# 叠加 Copernicus DEM、提取多时相作物物候特征、训练模型并对大区域栅格进行作物识别制图的完整流水线。
# 注：首次联网下载与计算较大数据集时需保证网络畅通，并分配充足的内存。
# ==============================================================================

# 1. 确保工作目录为项目根目录
if (!file.exists("data/ct_chile") && file.exists("../data/ct_chile")) {
  setwd("..")
}

# 2. 加载所需 R 包
library(tibble)
library(sits)
library(sf)
library(terra)
library(gdalcubes)
library(randomForest)

# 3. 设定输出目录
dir_work <- "datacube_crop_classification/work_dir"
dir_out  <- file.path(dir_work, "Out")
dir_cube <- file.path(dir_work, "Cube")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_cube, recursive = TRUE, showWarnings = FALSE)

# 4. 加载地面调查样本与感兴趣区 (ROI)
roi_path     <- "data/ct_chile/ROI/ROI_19HBA.shp"
samples_path <- "data/ct_chile/Samples/ground_data_19HBA_Chile_eng.shp"

if (!file.exists(roi_path) || !file.exists(samples_path)) {
  stop("未找到本地 ROI 或 Samples 矢量文件，请检查 data/ct_chile 目录。")
}

roi_test   <- sf::st_read(roi_path)
points_roi <- sf::st_read(samples_path)

# 5. 从云端 (MPC) 构建非规则 Sentinel-2 数据立方体
# 时间范围跨越农业生长季 (2020-05-01 至 2021-05-30)
message(">>> 步骤 1/7: 从云端 MPC 检索并构建 Sentinel-2 多时相数据立方体...")
sent_cube <- sits_cube(
  source     = "MPC",                           
  collection = "SENTINEL-2-L2A",                
  tiles      = "19HBA",                          
  bands      = c("B04", "B08", "B11", "CLOUD"), 
  start_date = "2020-05-01",                    
  end_date   = "2021-05-30"
)

# 6. 数据立方体时空规整化 (Regularization)
# 消除云污染，重采样至统一分辨率，并统一按 16 天间隔插值合成
message(">>> 步骤 2/7: 规整化数据立方体 (P16D 周期合成与去云处理)...")
sent_cube_reg <- sits_regularize(
  cube       = sent_cube,                  
  period     = "P16D",                         
  res        = 10,                             
  roi        = roi_test,                       
  output_dir = dir_cube,
  multicores = 4
)

# 7. 构建并叠加 DEM 高程数据立方体 (辅助区分山区与平原作物)
message(">>> 步骤 3/7: 检索并规整 Copernicus 30m DEM 数据立方体...")
dem_cube <- sits_cube(
  source     = "MPC",             
  collection = "COP-DEM-GLO-30",  
  tiles      = "19HBA"
)

dem_cube_reg <- sits_regularize(
  cube       = dem_cube,
  res        = 10,
  roi        = roi_test,
  output_dir = dir_cube,
  multicores = 4
)

# 合并多时相光学波段与高程波段
cube_19HBA <- sits_merge(sent_cube_reg, dem_cube_reg)

# 8. 计算植被指数 (NDVI: (B08 - B04)/(B08 + B04))
message(">>> 步骤 4/7: 计算时序 NDVI 植被指数...")
cube_19HBA <- sits_apply(
  cube_19HBA,
  NDVI       = ((B08 - B04) / (B08 + B04)),
  output_dir = dir_cube,
  multicores = 4
)

# 9. 采样提取多时相光谱与物候时序数据
message(">>> 步骤 5/7: 从多时相数据立方体中根据样点提取时间序列...")
samples_ts <- sits_get_data(
  cube       = cube_19HBA,            
  samples    = points_roi,              
  multicores = 4
)

# 10. 训练随机森林分类模型 (针对冬季作物、春季作物、果树等)
message(">>> 步骤 6/7: 训练随机森林多时相分类器...")
rfor_model <- sits_train(
  samples   = samples_ts,
  ml_method = sits_rfor(num_trees = 100)
)

# 11. 对数据立方体执行栅格分类预测与后处理
message(">>> 步骤 7/7: 对多时相数据立方体执行全像素识别与制图...")
probs_cube <- sits_classify(
  data       = cube_19HBA,              
  ml_model   = rfor_model,   
  roi        = roi_test,
  output_dir = dir_out,
  version    = "online_v1",
  multicores = 4,
  memsize    = 8
)

# 空间平滑
smooth_cube <- sits_smooth(
  cube           = probs_cube,
  window_size    = 7,
  neigh_fraction = 0.50,
  output_dir     = dir_out
)

# 输出最终分类 GeoTIFF 地图
class_cube <- sits_label_classification(
  cube       = smooth_cube,
  output_dir = dir_out,
  version    = "final_map"
)

message(">>> 在线端到端作物数据立方体识别全流程执行完毕！")
message(">>> 输出分类地图位置: ", dir_out)
