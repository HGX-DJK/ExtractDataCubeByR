# 多时相数据立方体作物识别与制图系统 (优化升级版)

本项目代码提取自联合国粮农组织（FAO）与联合国统计司合著的《遥感农业统计手册》（*UN Handbook on Remote Sensing for Agricultural Statistics*），专门针对**地球观测数据立方体（Earth Observation Data Cubes）**与**鹏城星云 iEarth SDC30 全球30米无缝数据立方体**进行农作物识别、高精度制图与空间质量诊断。

---

## 一、 文件清单与功能说明

| 文件名 | 类型 | 说明 |
| :--- | :--- | :--- |
| `06_train_crop_model.R` | R 脚本 | **核心训练与制图系统 (优化升级版)**。结合 NDVI/EVI/MNDWI/LSWI/NDBI/NDTI 等 12 维特征，多线程 `ranger` 训练，3×3 空间众数滤波去噪，输出全景二值制图与四宫格质检对比图。 |
| `05_pcl_sdc30_crop_classification.R` | R 脚本 | **SDC30 快速适配脚本**。已升级支持自动读取 CSV/SHP 样点、多光谱分类与空间平滑制图。 |
| `01_crop_classification_offline.R` | R 脚本 | **智利研究区离线演示脚本**。展示多时相时序模型、贝叶斯平滑与 22 类农作物土地覆盖分类。 |
| `02_crop_datacube_online_pipeline.R` | R 脚本 | **端到端在线全流程**。从微软行星计算机 (MPC) 检索多时相 Sentinel-2 影像、16天规整化合成与制图。 |
| `03_deafrica_datacube_crop_mapping.R` | R 脚本 | 基于 **Digital Earth Africa** 开放数据立方体（GeoMAD 复合产品）的卢旺达作物识别脚本。 |
| `04_multicube_data_extraction.R` | R 脚本 | **多数据立方体数据提取专有工具**。演示多源立方体融合（Sentinel-2 + DEM）、点/面地块抽样与表格导出。 |
| `prepare_samples.py` | Python 脚本 | **统一真值切片制作与高纯度样点提取工具**。智能自适应：本地有切片则秒级提取纯净样点，无切片则自动从大图裁剪对齐；支持步长调参（`--step`）。 |
| `install_dependencies.R` | R 脚本 | 一键安装所有依赖 R 包（`sits`, `terra`, `sf`, `ranger`, `tidyr`, `readr` 等）。 |

---

## 二、 核心算法与升级特性

1. **多维特征工程 (Feature Engineering)**：
   - 原始反射率：`Blue`, `Green`, `Red`, `NIR`, `SWIR1`, `SWIR2`
   - 植被生长与生物量：`NDVI`, `EVI`
   - 地表水分与沟渠水体：`MNDWI`, `LSWI`
   - 城乡建筑与不透水面：`NDBI`
   - 耕作/麦茬/秸秆残留：`NDTI`（归一化耕作指数，夏收夏种关键特征）
2. **训练加速与空间纯化**：
   - **样本纯化**：形态学侵蚀剥离地块边缘 30 米混合像元，降低边界噪声；
   - **算力加速**：集成 C++ 级 `ranger` 随机森林，训练耗时从数分钟降至 **15 秒**；
3. **空间上下文滤波 (Spatial Smoothing)**：
   - 预测后采用 3×3 空间众数滤波平滑（Modal Filter），消除椒盐噪声，农田地块连续性与交并比显著提高；
4. **可视化质检诊断**：
   - 每次训练自动导出 `output_sdc30/SDC30_Crop_Classification_Diagnostic.png` 高清四宫格质检对比图：
     - ① SDC30 标准假彩色合成 (NIR-Red-Green)
     - ② 地面参考真值底图
     - ③ 作物识别预测成果图 (空间平滑后)
     - ④ 空间误差差分图 (绿=TP正确耕地, 灰=TN正确非耕地, 红=FP虚报, 橙=FN漏报)

---

## 三、 标准运行流程

所有脚本均采用纯相对路径与自适应环境设计，无任何绝对路径绑定，可以在任意路径、任意终端或 IDE（RStudio / VS Code）中执行：

### 步骤 1：准备训练样点（支持【按切片裁剪】与【全图不裁剪】双模式）
```bash
# 模式 A：【针对大图进行裁剪】(局部切片模式)
# 以当前 SDC30 卫星切片几何为边界，从大图中裁剪出 100km×100km 局部切片并提取纯核样点：
python prepare_samples.py --crop --step 10

# 模式 B：【针对大图不进行裁剪】(全图宏观模式)
# 完全不按切片裁剪！直接获取整张大图 (如全国/全域 420 亿像元) 所有地理范围的宏观训练样本：
python prepare_samples.py --no-crop --samples 100000

# 模式 C (默认智能模式)：
# 本地已有切片底图则直接秒级抽样；未检测到切片底图则自动从大图裁剪保存：
python prepare_samples.py --step 10
```
*注：局部模式支持 `--step` 步长与形态学侵蚀；全图不裁剪模式支持 `--samples` 指定目标全域样本规模。*

### 步骤 2：执行模型训练、全幅制图与质检评估
```bash
Rscript 06_train_crop_model.R
```
*也可以直接在 **RStudio** 中打开 `06_train_crop_model.R` 点击 **Source** 或逐步调试运行。*

### 步骤 3：查看产出成果
- **作物分类二值地图**：`output_sdc30/SDC30_V003_50SMF_20210618_CropTypeMap.tif`
- **四宫格质检诊断图**：`output_sdc30/SDC30_Crop_Classification_Diagnostic.png`
- **训练好的模型对象**：`models/crop_rf_model.rds`
