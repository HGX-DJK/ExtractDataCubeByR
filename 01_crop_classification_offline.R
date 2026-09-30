# ==============================================================================
# 01_crop_classification_offline.R
# 联合国粮农组织《遥感农业统计手册》- 多时相数据立方体作物识别 (离线快速运行版)
# ==============================================================================
# 说明：
# 本脚本使用项目自带的本地预处理数据（位于 data/ct_chile/），无需重新下载数十GB的卫星数据。
# 演示从“多时相光谱样本提取 -> 时序模型训练 -> 数据立方体分类 -> 贝叶斯后处理平滑 -> 作物制图”全流程。
# ==============================================================================

# 1. 确保工作目录为项目根目录 (如果从子目录运行则自动切换)
if (!file.exists("data/ct_chile") && file.exists("../data/ct_chile")) {
  setwd("..")
}
message(">>> 当前工作目录: ", getwd())

# 2. 加载所需 R 依赖包
required_packages <- c("tibble", "sits", "sf", "terra", "randomForest")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("缺少必要包: '", pkg, "'。请先运行 install_dependencies.R 进行安装。"))
  }
  library(pkg, character.only = TRUE)
}

# 3. 设置输出目录
output_dir <- "datacube_crop_classification/output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 4. 设置农作物与土地覆盖分类配色表
cl_tbl_eng <- tibble::tibble(name = character(), color = character()) |>
  tibble::add_row(name = "Sand_dunes",           color = "#bababa") |>
  tibble::add_row(name = "Fallow",               color = "#7D4A0C") |>
  tibble::add_row(name = "Deciduous_forest",     color = "#228710") |>
  tibble::add_row(name = "Perennial_forest",     color = "#004529") |>
  tibble::add_row(name = "Winter_crops",         color = "#ffff99") |>
  tibble::add_row(name = "Spring_crops",         color = "#ffd92f") |>
  tibble::add_row(name = "Successive_crops",     color = "#FFA500") |>
  tibble::add_row(name = "Deciduous_fruit_tree", color = "#df65b0") |>
  tibble::add_row(name = "Perennial_fruit_tree", color = "#ce1256") |>
  tibble::add_row(name = "Wetlands",             color = "#93dfe6") |>
  tibble::add_row(name = "Shrublands",           color = "#8B864E") |>
  tibble::add_row(name = "Snow_glaciers",        color = "#00BFFF") |>
  tibble::add_row(name = "Annual_pasture",       color = "#6959CD") |>
  tibble::add_row(name = "Perennial_pasture",    color = "#B452CD") |>
  tibble::add_row(name = "Mature_plantation",    color = "#7CFC00") |>
  tibble::add_row(name = "Clear_cut_plantation", color = "#bf812d") |>
  tibble::add_row(name = "Young_plantation",     color = "#0ecf65") |>
  tibble::add_row(name = "Mountain_grasslands",  color = "#35978f") |>
  tibble::add_row(name = "Grasslands",           color = "#01665e") |>
  tibble::add_row(name = "Bare_soils",           color = "#e0e0e0") |>
  tibble::add_row(name = "Water_bodies",         color = "#0C3BD4") |>
  tibble::add_row(name = "Artificial_surfaces",  color = "#fa0000")

sits_colors_set(cl_tbl_eng)
message(">>> 作物分类配色表加载完毕。")

# 5. 加载多时相卫星时序样本数据 (Sentinel-2 + DEM)
samples_file <- "data/ct_chile/balanced_samples.rds"
if (file.exists(samples_file)) {
  message(">>> 正在加载经过质控与类别平衡的多时相时序样本: ", samples_file)
  balanced_samples <- readRDS(samples_file)
  print(summary(balanced_samples))
} else {
  stop("未找到时序样本文件: ", samples_file)
}

# 6. 训练或加载机器学习时序分类模型 (Random Forest)
model_file <- "data/ct_chile/rfor_model.rds"
if (file.exists(model_file)) {
  message(">>> 检测到预训练好的多时相随机森林模型，直接加载: ", model_file)
  rfor_model <- readRDS(model_file)
} else {
  message(">>> 正在基于多时相样本训练随机森林模型 (sits_rfor)...")
  rfor_model <- sits_train(
    samples   = balanced_samples,
    ml_method = sits_rfor(num_trees = 100)
  )
  saveRDS(rfor_model, file.path(output_dir, "my_rfor_model.rds"))
}

# 7. 对多时相数据立方体进行分类预测 (sits_classify)
# 检查本地是否有保存的概率立方体或原始规整化立方体
probs_file <- "data/ct_chile/probs_19HBA.rds"
if (file.exists(probs_file)) {
  message(">>> 正在加载多时相数据立方体分类概率结果: ", probs_file)
  probs_19HBA <- readRDS(probs_file)
} else if (file.exists("data/ct_chile/sent_19HBA_reg.rds")) {
  message(">>> 正在对规整化多时相数据立方体进行逐像素时序分类预测...")
  sent_reg <- readRDS("data/ct_chile/sent_19HBA_reg.rds")
  probs_19HBA <- sits_classify(
    data       = sent_reg,
    ml_model   = rfor_model,
    output_dir = output_dir,
    version    = "v_offline",
    multicores = 2,
    memsize    = 4
  )
}

# 8. 贝叶斯后处理空间平滑 (去除椒盐噪声，结合空间邻域上下文)
smooth_file <- "data/ct_chile/smooth_19HBA.rds"
if (file.exists(smooth_file)) {
  message(">>> 加载贝叶斯平滑概率立方体: ", smooth_file)
  smooth_19HBA <- readRDS(smooth_file)
} else if (exists("probs_19HBA")) {
  message(">>> 计算贝叶斯空间平滑...")
  smooth_19HBA <- sits_smooth(
    cube           = probs_19HBA,
    window_size    = 7,
    neigh_fraction = 0.50,
    output_dir     = output_dir
  )
}

# 9. 生成最终作物分类专题图 (Labeling)
class_file <- "data/ct_chile/class_19HBA.rds"
if (file.exists(class_file)) {
  message(">>> 加载最终作物分类图对象: ", class_file)
  class_19HBA <- readRDS(class_file)
} else if (exists("smooth_19HBA")) {
  message(">>> 提取最高概率类别，生成最终作物分类图...")
  class_19HBA <- sits_label_classification(
    cube       = smooth_19HBA,
    output_dir = output_dir,
    version    = "final_map"
  )
}

# 10. 输出分类信息与摘要，保存分析成果
if (exists("class_19HBA")) {
  message(">>> 作物分类识别成功完成！")
  labels_vec <- sits_labels(class_19HBA)
  message(">>> 包含的作物与土地覆盖类别 (共 ", length(labels_vec), " 类): ")
  print(labels_vec)
  
  # 导出类别台账与配色表到 output 目录
  cl_export <- cl_tbl_eng |>
    dplyr::filter(name %in% labels_vec)
  write.csv(cl_export, file.path(output_dir, "crop_classes_summary.csv"), row.names = FALSE)
  saveRDS(rfor_model, file.path(output_dir, "crop_rfor_model.rds"))
  message(">>> 成果已保存至: ", normalizePath(output_dir))
  message("    - 类别色系表: ", file.path(output_dir, "crop_classes_summary.csv"))
  message("    - 随机森林模型: ", file.path(output_dir, "crop_rfor_model.rds"))
}

message(">>> 脚本运行结束。")
