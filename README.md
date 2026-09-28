# 多时相数据立方体作物识别 R 代码项目

本项目代码提取自联合国粮农组织（FAO）与联合国统计司合著的《遥感农业统计手册》（*UN Handbook on Remote Sensing for Agricultural Statistics*），专门针对**地球观测数据立方体（Earth Observation Data Cubes）**进行农作物识别与制图（Crop Type Mapping）。

---

## 一、 文件清单

| 文件名 | 类型 | 说明 |
| :--- | :--- | :--- |
| `01_crop_classification_offline.R` | R 脚本 | **推荐首选**。离线快速演示脚本，直接读取项目 `data/ct_chile` 内现成的数据立方体与样本，无需耗时下载数 GB 影像。 |
| `02_crop_datacube_online_pipeline.R` | R 脚本 | **端到端在线全流程**。从微软行星计算机 (MPC) 检索 Sentinel-2 影像、16天时间规则化合成、DEM融合、时序采样、训练分类与大尺度制图。 |
| `03_deafrica_datacube_crop_mapping.R` | R 脚本 | 基于 **Digital Earth Africa** 开放数据立方体（GeoMAD 复合产品）的卢旺达作物识别脚本。 |
| `04_multicube_data_extraction.R` | R 脚本 | **多数据立方体数据提取专有工具**。演示多源立方体融合（Sentinel-2 + DEM + 雷达）、点/面地块抽样、跨多瓦片提取及导出 CSV。 |
| `05_pcl_sdc30_crop_classification.R` | R 脚本 | **鹏城星云 iEarth SDC30 数据集专用适配脚本**。演示如何直接接入该平台的全球30米无缝数据立方体进行作物识别。 |
| `06_train_crop_model.R` | R 脚本 | **作物分类模型专用训练与验证脚本**。加载样点、训练 Random Forest、输出混淆矩阵与波段重要性、保存模型并全景制图。 |
| `install_dependencies.R` | R 脚本 | 一键安装所有依赖 R 包（`sits`, `terra`, `sf`, `gdalcubes` 等）。 |

---

## 二、 核心算法与识别作物

1. **核心框架**：基于 R 语言的 `sits` (Satellite Image Time Series) 框架，采用 “Time-first, space-later”（先时序、后空间）范式：
   - **时间维**：利用 Sentinel-2 / Landsat 多时相光谱曲线捕获作物物候生长期特征；
   - **空间维**：利用贝叶斯空间平滑（`sits_smooth`）去除斑点噪声，结合地块与空间邻域。
2. **主要识别作物与覆盖类型**：
   - **冬季作物 (Winter crops)**
   - **春季作物 (Spring crops)**
   - **接茬/连作作物 (Successive crops)**
   - **落叶果树 (Deciduous fruit tree)** 与 **常绿果树 (Perennial fruit tree)**
   - **一年生草场 (Annual pasture)** 与 **多年生草场 (Perennial pasture)**
   - 以及休耕地（Fallow）、森林、水体等背景类。

---

## 三、 运行指南

### 1. 安装 R 语言环境
如果你的电脑尚未安装 R 语言，请前往官网下载并安装：
- [CRAN R 官方下载](https://cloud.r-project.org/)（推荐 R 4.3 或更高版本）
- 可选安装 [RStudio Desktop](https://posit.co/download/rstudio-desktop/)（最推荐的交互式运行环境）

### 2. 安装依赖包
在 R 控制台或命令行中运行：
```bash
Rscript install_dependencies.R
```
或者在 RStudio 中打开 `install_dependencies.R` 并运行。

### 3. 运行作物识别脚本

#### 方式 A：在 RStudio 中运行（最方便）
1. 双击打开根目录的 `UN-Handbook.Rproj`；
2. 打开 `datacube_crop_classification/01_crop_classification_offline.R`；
3. 点击右上角 **Source** 或逐行点击 **Run** 即可查看数据立方体时序曲线和分类制图结果。

#### 方式 B：在命令行运行
在项目根目录下执行：
```bash
Rscript datacube_crop_classification/01_crop_classification_offline.R
```
分类结果将保存在 `datacube_crop_classification/output/` 文件夹中。
