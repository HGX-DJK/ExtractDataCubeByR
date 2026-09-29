# -*- coding: utf-8 -*-
"""
==============================================================================
  从地面真值栅格 GroundTruth.tif 提取所有有效像元坐标，
  转换为 WGS84 经纬度，导出为 my_crop_samples.csv
  供 06_train_crop_model.R 直接读取使用
==============================================================================
  输出格式（与 R 脚本完全兼容）:
    lon   - WGS84 经度
    lat   - WGS84 纬度
    label - "Cropland" 或 "Non_Cropland"
==============================================================================
"""

import os
import sys
import io
import csv
import math
import numpy as np
import rasterio

# 兼容 Windows 控制台中文输出
if sys.platform.startswith("win"):
    try:
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
        sys.stderr.reconfigure(encoding='utf-8', errors='replace')
    except Exception:
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
        sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

# ──────────────────────────── 路径配置 ────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
GT_TIF     = os.path.join(SCRIPT_DIR, "data", "SDC30_V003_50SMF_20210618_GroundTruth.tif")
OUT_CSV    = os.path.join(SCRIPT_DIR, "data", "my_crop_samples.csv")

# ──────────────────────────── 采样控制 ────────────────────────────
# 全量像元约 1360 万个，R randomForest 建议 ≤20 万样本。
# 此处用"规则格网"系统采样（非随机），保证均匀覆盖全幅影像：
#   STEP_SIZE=10  → ~13.4 万点（推荐默认）
#   STEP_SIZE=6   → ~37 万点
#   STEP_SIZE=1   → 全量 1360 万点（CSV ~1.3 GB，需 R 内存 64 GB+）
STEP_SIZE = 10

# ──────────────── UTM zone 50N → WGS84 纯数学转换 ────────────────
# 不依赖 PROJ 数据库（proj.db），直接用 WGS84 椭球参数计算。
# 参考：Snyder 1987 "Map Projections – A Working Manual", p.60-64

def utm_to_wgs84(easting, northing, zone_number=50, northern=True):
    """
    将 UTM 坐标数组批量转换为 WGS84 经纬度（向量化，纯 Python/NumPy）
    easting, northing : numpy 数组，单位米
    返回 (lons, lats) : numpy 数组，单位度
    """
    k0 = 0.9996
    a  = 6378137.0          # WGS84 长半轴
    e  = 0.0818191908426215  # WGS84 第一离心率
    e1sq = e**2 / (1 - e**2)
    lon0 = math.radians((zone_number - 1) * 6 - 180 + 3)

    x = np.asarray(easting,  dtype=np.float64) - 500000.0
    y = np.asarray(northing, dtype=np.float64)
    if not northern:
        y -= 10000000.0

    M  = y / k0
    mu = M / (a * (1 - e**2/4 - 3*e**4/64 - 5*e**6/256))
    e1 = (1 - np.sqrt(1 - e**2)) / (1 + np.sqrt(1 - e**2))

    phi1 = (mu
            + (3*e1/2 - 27*e1**3/32)   * np.sin(2*mu)
            + (21*e1**2/16 - 55*e1**4/32) * np.sin(4*mu)
            + (151*e1**3/96)            * np.sin(6*mu)
            + (1097*e1**4/512)          * np.sin(8*mu))

    N1 = a / np.sqrt(1 - (e * np.sin(phi1))**2)
    T1 = np.tan(phi1)**2
    C1 = e1sq * np.cos(phi1)**2
    R1 = a * (1 - e**2) / (1 - (e * np.sin(phi1))**2)**1.5
    D  = x / (N1 * k0)

    lat_rad = phi1 - (N1 * np.tan(phi1) / R1) * (
          D**2/2
        - (5 + 3*T1 + 10*C1 - 4*C1**2 - 9*e1sq) * D**4/24
        + (61 + 90*T1 + 298*C1 + 45*T1**2 - 252*e1sq - 3*C1**2) * D**6/720
    )
    lon_rad = lon0 + (
          D
        - (1 + 2*T1 + C1) * D**3/6
        + (5 - 2*C1 + 28*T1 - 3*C1**2 + 8*e1sq + 24*T1**2) * D**5/120
    ) / np.cos(phi1)

    return np.degrees(lon_rad), np.degrees(lat_rad)


# ──────────────────────────── 主函数 ──────────────────────────────

def export_samples():
    if not os.path.exists(GT_TIF):
        print(f"[错误] 未找到真值底图：{GT_TIF}")
        print("请先运行 export_ground_truth.py 生成真值底图！")
        return

    print("=================================================================")
    print(">>> 正在从全量地面真值底图中提取训练坐标样点...")
    print(f"    * 输入真值图  : {os.path.basename(GT_TIF)}")
    print(f"    * 输出样点文件: {os.path.basename(OUT_CSV)}")
    print(f"    * 格网步长    : 每 {STEP_SIZE} 个像元取1点（均匀覆盖全图，非随机）")
    print("=================================================================")

    with rasterio.open(GT_TIF) as src:
        arr = src.read(1)   # 读入内存（0.94 MB，秒级）
        tfm = src.transform

    rows, cols = arr.shape

    # 格网行列索引
    row_idx = np.arange(0, rows, STEP_SIZE)
    col_idx = np.arange(0, cols, STEP_SIZE)
    rr, cc  = np.meshgrid(row_idx, col_idx, indexing='ij')
    rr = rr.ravel();  cc = cc.ravel()

    # 像元中心的 UTM 50N 坐标
    xs = tfm.c + (cc + 0.5) * tfm.a   # easting
    ys = tfm.f + (rr + 0.5) * tfm.e   # northing

    labels_raw = arr[rr, cc]
    valid      = labels_raw != 255
    xs, ys, labels_raw = xs[valid], ys[valid], labels_raw[valid]

    n_crop    = int(np.sum(labels_raw == 1))
    n_noncrop = int(np.sum(labels_raw == 0))
    n_total   = n_crop + n_noncrop
    print(f">>> 格网采样完成（步长={STEP_SIZE}，共 {n_total:,} 个有效样点）:")
    print(f"    * 耕地    (Cropland)    : {n_crop:,}  ({n_crop/n_total*100:.1f}%)")
    print(f"    * 非耕地  (Non_Cropland): {n_noncrop:,}  ({n_noncrop/n_total*100:.1f}%)")

    # UTM 50N → WGS84 lon/lat（纯数学，无需 PROJ 数据库）
    print(">>> 正在批量转换坐标（UTM 50N → WGS84 lon/lat）...")
    lons, lats = utm_to_wgs84(xs, ys, zone_number=50, northern=True)

    # 写出 CSV
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    print(f">>> 正在写出 CSV（{n_total:,} 行）...")
    with open(OUT_CSV, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(["lon", "lat", "label"])
        for lon, lat, lbl in zip(lons, lats, labels_raw):
            writer.writerow([f"{lon:.7f}", f"{lat:.7f}",
                             "Cropland" if lbl == 1 else "Non_Cropland"])

    size_mb = os.path.getsize(OUT_CSV) / 1024 / 1024
    print("=================================================================")
    print(">>> 样点文件导出完成！")
    print(f"    * 文件路径 : {OUT_CSV}")
    print(f"    * 文件大小 : {size_mb:.1f} MB")
    print(f"    * 总样点数 : {n_total:,} 个（耕地 {n_crop:,} + 非耕地 {n_noncrop:,}）")
    print(f"    * 经度范围 : {lons.min():.4f}E ~ {lons.max():.4f}E")
    print(f"    * 纬度范围 : {lats.min():.4f}N ~ {lats.max():.4f}N")
    print("=================================================================")
    print(">>> 下一步：运行 R 训练脚本（会自动加载此样点文件）：")
    print("    Rscript 06_train_crop_model.R")
    print("=================================================================")


if __name__ == "__main__":
    export_samples()
