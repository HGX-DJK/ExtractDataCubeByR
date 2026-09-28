# ==============================================================================
# install_dependencies.R
# 自动安装多时相数据立方体作物识别所需的 R 依赖包 (自动处理用户权限与国内高速镜像)
# ==============================================================================

# 1. 设置用户个人的独立 R 包目录 (彻底避免 C:/Program Files 写入权限不足报错)
r_ver <- sprintf("%s.%s", R.version$major, substr(R.version$minor, 1, 1))
user_lib <- file.path(Sys.getenv("LOCALAPPDATA"), "R", "win-library", r_ver)

if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
}
.libPaths(c(user_lib, .libPaths()))

message("=================================================================")
message(">>> R 包安装目标目录 (已获写入权限): ", user_lib)
message(">>> 镜像源: 清华大学开源软件镜像站 (CRAN)")
message("=================================================================")

options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

required_packages <- c(
  "tibble",
  "dplyr",
  "sf",
  "terra",
  "randomForest",
  "gdalcubes",
  "rstac",
  "stars",
  "sits"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE, lib.loc = user_lib)) {
    message(">>> [正在安装]: ", pkg, " ... (请稍候)")
    tryCatch({
      install.packages(pkg, lib = user_lib, dependencies = TRUE)
      message(">>> [成功]: ", pkg, " 安装完成！")
    }, error = function(e) {
      message(">>> [错误]: ", pkg, " 安装失败: ", e$message)
    })
  } else {
    message(">>> [OK]: ", pkg, " 已安装在本地。")
  }
}

message("=================================================================")
message(">>> 依赖检查与安装完毕！现在可以重新运行 05_pcl_sdc30_crop_classification.R 了。")
message("=================================================================")
