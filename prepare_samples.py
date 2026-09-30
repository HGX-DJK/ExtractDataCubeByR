# -*- coding: utf-8 -*-
"""
==============================================================================
  prepare_samples.py: 地面真值切片制备与作物样点提取统一工具 (全域/局部双模式版)
==============================================================================
  核心功能：
    针对大尺度参考底图 (如 42GB 全国 gengdi.tif 或区域大图)，提供两种采样范式：

    1. 【针对大图进行裁剪】(--crop / 默认切片模式)：
       - 以当前 SDC30 卫星切片 (如 50SMF 瓦片范围) 为几何边界，从大图中裁剪出
         局部切片底图 GroundTruth.tif，提取该切片范围内的局部训练样点；
       - 适用场景：针对当前特定瓦片进行高精度制图与质检对比。

    2. 【针对大图不进行裁剪】(--no-crop / --full-image)：
       - 【完全不按卫星切片做任何裁剪】！
       - 直接在整张大图的全部空间幅员内 (如全国 73°E~135°E, 4°N~53°N) 执行
         全图宏观抽样，直接获取覆盖整张大图所有区域的跨区/全域宏观训练样本！
       - 适用场景：构建大尺度跨区通用训练集、全域作物样本台账。
==============================================================================
  常用命令：
    • 模式 A (局部切片模式，针对大图裁剪)：
        python prepare_samples.py --crop --step 10
    • 模式 B (全域大图模式，不裁剪大图，直接获取大图全域样本)：
        python prepare_samples.py --no-crop --samples 100000
    • 默认智能模式 (已有局部切片直接用，无切片则自动裁剪)：
        python prepare_samples.py --step 10
==============================================================================
"""

import os
import sys
import io
import math
import argparse
import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling, transform_bounds, transform
from rasterio.windows import from_bounds
import pandas as pd
from scipy.ndimage import binary_erosion

# 兼容 Windows 控制台输出编码
if sys.platform.startswith("win"):
    try:
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
        sys.stderr.reconfigure(encoding='utf-8', errors='replace')
    except Exception:
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
        sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

# 自动配置 PROJ 数据库环境变量
def setup_proj_env():
    candidate_proj_dirs = [
        r"D:\software\python\Lib\site-packages\rasterio\proj_data",
        os.path.join(os.path.dirname(sys.executable), "Lib", "site-packages", "rasterio", "proj_data"),
        r"D:\software\EarthVisLabApps\cesiumlab\4.0.15\tools\proj_data",
        r"C:\Program Files\PostgreSQL\share\contrib\postgis\proj"
    ]
    for d in candidate_proj_dirs:
        if os.path.isdir(d) and os.path.exists(os.path.join(d, "proj.db")):
            os.environ['PROJ_DATA'] = d
            os.environ['PROJ_LIB'] = d
            return True
    return False

# 纯数学解析转换兜底
def fallback_utm_to_wgs84(easting, northing, zone_number=50, northern=True):
    k0 = 0.9996
    a  = 6378137.0
    e  = 0.0818191908426215
    e1sq = e**2 / (1 - e**2)
    lon0 = math.radians((zone_number - 1) * 6 - 180 + 3)

    x = np.asarray(easting, dtype=np.float64) - 500000.0
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


def get_default_ref_path():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.environ.get("REF_DATA_PATH"),
        os.path.join(script_dir, "data", "gengdi.tif"),
        os.path.join(script_dir, "..", "data", "gengdi.tif"),
        r"E:\agriculture\gaced30_validation_pipeline\data\reference_data\gengdi.tif"
    ]
    return next((c for c in candidates if c and os.path.exists(c)), None)


# ==============================================================================
# 模式 A: 针对大图进行裁剪 (从大图裁剪切片真值底图 GroundTruth.tif)
# ==============================================================================
def crop_ground_truth_tile(sdc_path, ref_path, out_gt_path):
    if not os.path.exists(sdc_path):
        print(f"[错误] 未找到目标 SDC30 影像文件: {sdc_path}")
        return False

    if not ref_path or not os.path.exists(ref_path):
        ref_path = get_default_ref_path()

    if not ref_path or not os.path.exists(ref_path):
        print("[错误] 未找到大尺度参考真值底图 (如 gengdi.tif)！请通过 --ref 指定。")
        return False

    print("=================================================================")
    print(">>> 【裁剪模式】 正在从大底图中按当前切片范围裁剪对齐切片真值图...")
    print(f"    • 目标切片瓦片 : {os.path.basename(sdc_path)}")
    print(f"    • 来源参考大图 : {ref_path}")
    print(f"    • 输出切片真值 : {out_gt_path}")
    print("=================================================================")

    with rasterio.open(sdc_path) as sdc:
        sdc_crs = sdc.crs
        sdc_bounds = sdc.bounds
        sdc_shape = sdc.shape
        sdc_transform = sdc.transform
        meta = sdc.meta.copy()

    with rasterio.open(ref_path) as ref:
        ref_b = transform_bounds(sdc_crs, ref.crs, *sdc_bounds)
        ref_win = from_bounds(*ref_b, ref.transform)
        col_off = max(0, int(np.floor(ref_win.col_off)) - 10)
        row_off = max(0, int(np.floor(ref_win.row_off)) - 10)
        width = int(np.ceil(ref_win.width)) + 20
        height = int(np.ceil(ref_win.height)) + 20
        pad_win = rasterio.windows.Window(col_off, row_off, width, height)
        src_data = ref.read(1, window=pad_win)
        src_transform = ref.window_transform(pad_win)
        src_crs = ref.crs

    dst_data = np.full(sdc_shape, 255, dtype=np.uint8)
    reproject(
        source=src_data,
        destination=dst_data,
        src_transform=src_transform,
        src_crs=src_crs,
        dst_transform=sdc_transform,
        dst_crs=sdc_crs,
        resampling=Resampling.nearest,
        src_nodata=255,
        dst_nodata=255
    )

    dst_binary = np.where((dst_data == 1) | (dst_data == 10), np.uint8(1), np.where(dst_data == 0, np.uint8(0), np.uint8(255)))

    meta.update({
        'count': 1,
        'dtype': 'uint8',
        'nodata': 255,
        'compress': 'lzw'
    })

    os.makedirs(os.path.dirname(out_gt_path), exist_ok=True)
    with rasterio.open(out_gt_path, 'w', **meta) as dst:
        dst.write(dst_binary, 1)

    print(f">>> 真值底图切片制作成功: {out_gt_path} (大小: {os.path.getsize(out_gt_path) / 1024 / 1024:.2f} MB)")
    return True


# ==============================================================================
# 从局部切片真值底图中提取空间样点 (带侵蚀纯化)
# ==============================================================================
def extract_local_tile_samples(gt_path, out_csv_path, step_size=10, erode_edges=True, sdc_path=None, extract_features=True):
    if not os.path.exists(gt_path):
        print(f"[错误] 未找到切片真值底图: {gt_path}")
        return False

    print("=================================================================")
    print(">>> 正在从【局部切片真值底图】提取高纯度训练样点...")
    print(f"    • 输入切片文件 : {os.path.basename(gt_path)}")
    print(f"    • 输出样点文件 : {os.path.basename(out_csv_path)}")
    print(f"    • 格网采样步长 : 每 {step_size} 个像元取 1 点 (物理间距约 {step_size * 30} 米)")
    print(f"    • 边缘侵蚀纯化 : {'已启用 (剥离地块边界 30米 混合过渡像元)' if erode_edges else '未启用'}")
    if extract_features and sdc_path and os.path.exists(sdc_path):
        print(f"    • 特征直出加速 : 已启用 (秒级直出 6 波段 + 6 大核心遥感指数)")
    print("=================================================================")

    with rasterio.open(gt_path) as src:
        arr = src.read(1)
        tfm = src.transform
        crs = src.crs

    rows, cols = arr.shape
    crop_mask = (arr == 1)
    noncrop_mask = (arr == 0)

    if erode_edges:
        selem = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=bool)
        crop_pure = binary_erosion(crop_mask, structure=selem)
        noncrop_pure = binary_erosion(noncrop_mask, structure=selem)
        print(f">>> 空间纯净化完成: 剥离了 {(np.sum(crop_mask) - np.sum(crop_pure)):,} 个耕地边缘像元与 {(np.sum(noncrop_mask) - np.sum(noncrop_pure)):,} 个非耕地边缘像元。")
    else:
        crop_pure = crop_mask
        noncrop_pure = noncrop_mask

    row_idx = np.arange(0, rows, step_size)
    col_idx = np.arange(0, cols, step_size)
    rr, cc = np.meshgrid(row_idx, col_idx, indexing='ij')
    rr = rr.ravel()
    cc = cc.ravel()

    is_crop = crop_pure[rr, cc]
    is_noncrop = noncrop_pure[rr, cc]
    valid = is_crop | is_noncrop

    rr_v = rr[valid]
    cc_v = cc[valid]
    labels = np.where(is_crop[valid], "Cropland", "Non_Cropland")

    # ──────────────────────────────────────────────────────────────────────────
    # 特征直出与物理质量过滤 (毫秒级 NumPy 矢量化采样)
    # ──────────────────────────────────────────────────────────────────────────
    feature_dict = {}
    if extract_features and sdc_path and os.path.exists(sdc_path):
        print(">>> 正在以 NumPy 毫秒级矩阵切片直出多光谱波段与遥感指数 (NDVI, EVI, MNDWI, LSWI, NDBI, NDTI)...")
        with rasterio.open(sdc_path) as sdc:
            bands_data = sdc.read()  # (6, H, W)
            
        b2 = bands_data[0, rr_v, cc_v]
        b3 = bands_data[1, rr_v, cc_v]
        b4 = bands_data[2, rr_v, cc_v]
        b8 = bands_data[3, rr_v, cc_v]
        b11 = bands_data[4, rr_v, cc_v]
        b12 = bands_data[5, rr_v, cc_v]

        # 遥感物理合理性清洗：剔除无数据黑边像元 (全0) 与极端异常饱和值
        valid_spec = (b2 > 0) & (b4 > 0) & (b8 > 0) & (b2 < 9000)
        n_filtered = int(np.sum(~valid_spec))
        if n_filtered > 0:
            print(f">>> 光谱物理清洗: 剔除了 {n_filtered:,} 个黑边无数据或异常饱和像元。")

        rr_v = rr_v[valid_spec]
        cc_v = cc_v[valid_spec]
        labels = labels[valid_spec]
        b2 = b2[valid_spec]
        b3 = b3[valid_spec]
        b4 = b4[valid_spec]
        b8 = b8[valid_spec]
        b11 = b11[valid_spec]
        b12 = b12[valid_spec]

        # 归一化地表反射率 (0 ~ 1)
        b = b2.astype(np.float32) / 10000.0
        g = b3.astype(np.float32) / 10000.0
        r = b4.astype(np.float32) / 10000.0
        nir = b8.astype(np.float32) / 10000.0
        s1 = b11.astype(np.float32) / 10000.0
        s2 = b12.astype(np.float32) / 10000.0

        eps = 1e-6
        ndvi = (nir - r) / (nir + r + eps)
        evi = 2.5 * (nir - r) / (nir + 6.0 * r - 7.5 * b + 1.0 + eps)
        mndwi = (g - s1) / (g + s1 + eps)
        lswi = (nir - s1) / (nir + s1 + eps)
        ndbi = (s1 - nir) / (s1 + nir + eps)
        ndti = (s1 - s2) / (s1 + s2 + eps)

        feature_dict = {
            'Blue': b2,
            'Green': b3,
            'Red': b4,
            'NIR': b8,
            'SWIR1': b11,
            'SWIR2': b12,
            'NDVI': np.round(ndvi, 6),
            'EVI': np.round(evi, 6),
            'MNDWI': np.round(mndwi, 6),
            'LSWI': np.round(lswi, 6),
            'NDBI': np.round(ndbi, 6),
            'NDTI': np.round(ndti, 6),
        }

    xs = tfm.c + (cc_v + 0.5) * tfm.a
    ys = tfm.f + (rr_v + 0.5) * tfm.e

    setup_proj_env()
    try:
        lons, lats = transform(crs, 'EPSG:4326', xs, ys)
        lons = np.array(lons, dtype=np.float64)
        lats = np.array(lats, dtype=np.float64)
    except Exception:
        lons, lats = fallback_utm_to_wgs84(xs, ys, zone_number=50, northern=True)

    df_data = {'lon': np.round(lons, 7), 'lat': np.round(lats, 7), 'label': labels}
    df_data.update(feature_dict)
    df = pd.DataFrame(df_data)

    n_crop = int(np.sum(labels == "Cropland"))
    n_noncrop = int(np.sum(labels == "Non_Cropland"))
    crop_pct = n_crop / len(df) * 100 if len(df) > 0 else 0
    non_pct = n_noncrop / len(df) * 100 if len(df) > 0 else 0

    os.makedirs(os.path.dirname(out_csv_path), exist_ok=True)
    df.to_csv(out_csv_path, index=False)

    print("=================================================================")
    print(">>> 局部切片样点提取完毕！")
    print(f"    • 成果文件 : {out_csv_path} ({os.path.getsize(out_csv_path) / 1024 / 1024:.2f} MB)")
    print(f"    • 样点总数 : {len(df):,} 个 (耕地: {n_crop:,} [{crop_pct:.1f}%], 非耕地: {n_noncrop:,} [{non_pct:.1f}%])")
    print(f"    • 经纬范围 : [{lons.min():.4f}°E ~ {lons.max():.4f}°E], [{lats.min():.4f}°N ~ {lats.max():.4f}°N]")
    if feature_dict:
        print(f"    • 特征维度 : 包含 12 维完整遥感光谱与植被指数特征 (可直接供 R/Python 训练秒读)！")
    print("=================================================================")
    return True


# ==============================================================================
# 模式 B: 针对大图不进行裁剪 (全图大尺度直接抽样，跨越整个大图所有地理幅员)
# ==============================================================================
def extract_full_big_image_samples(ref_path, out_csv_path, target_sample_count=100000):
    if not ref_path or not os.path.exists(ref_path):
        ref_path = get_default_ref_path()

    if not ref_path or not os.path.exists(ref_path):
        print("[错误] 未找到大尺度参考真值底图 (如 gengdi.tif)！请通过 --ref 指定。")
        return False

    print("=================================================================")
    print(">>> 【全图大图模式】 针对大图【不进行任何切片裁剪】，直接获取大图全域训练样本！")
    print(f"    • 来源大图图层 : {ref_path}")
    print(f"    • 目标样本规模 : 约 {target_sample_count:,} 个样点 (全国/全域大尺度)")
    print(f"    • 输出样点文件 : {out_csv_path}")
    print("=================================================================")

    with rasterio.open(ref_path) as src:
        rows, cols = src.shape
        tfm = src.transform
        crs = src.crs
        bounds = src.bounds
        print(f">>> 大图规格: {rows:,} 行 × {cols:,} 列 (总计 {rows*cols/1e8:.2f} 亿像元)")
        print(f">>> 地理幅员: 经度 [{bounds.left:.2f}° ~ {bounds.right:.2f}°], 纬度 [{bounds.bottom:.2f}° ~ {bounds.top:.2f}°]")

        # 针对 42GB 超大型栅格采用高效流式分行扫描抽样
        num_scan_rows = max(100, int(np.sqrt(target_sample_count * 2)))
        step_row = max(1, rows // num_scan_rows)
        step_col = max(1, cols // num_scan_rows)

        sample_rows = np.arange(step_row // 2, rows, step_row)
        print(f">>> 正在以高并发流式扫描大图 (扫描 {len(sample_rows)} 条横断面，横向步长 {step_col})...")

        samples_list = []
        for idx, r in enumerate(sample_rows):
            win = rasterio.windows.Window(0, int(r), cols, 1)
            row_data = src.read(1, window=win)[0]

            col_indices = np.arange(step_col // 2, cols, step_col)
            vals = row_data[col_indices]
            valid = (vals == 1) | (vals == 0) | (vals == 10)

            c_v = col_indices[valid]
            v_v = vals[valid]

            if len(c_v) > 0:
                xs = tfm.c + (c_v + 0.5) * tfm.a
                ys = tfm.f + (r + 0.5) * tfm.e
                lbls = np.where((v_v == 1) | (v_v == 10), "Cropland", "Non_Cropland")
                samples_list.append(pd.DataFrame({'lon': xs, 'lat': ys, 'label': lbls}))

    if len(samples_list) == 0:
        print("[错误] 未在大图中检索到有效像元！")
        return False

    df_full = pd.concat(samples_list, ignore_index=True)

    # 若坐标不是地理坐标系 (如投影坐标 UTM)，则转换至 WGS84 经纬度；若已经是经纬度则跳过转换，免除 PROJ 警告
    if not crs.is_geographic:
        setup_proj_env()
        print(">>> 正在转换大图坐标至 WGS84 经纬度...")
        try:
            lons, lats = transform(crs, 'EPSG:4326', df_full['lon'].values, df_full['lat'].values)
            df_full['lon'] = np.round(lons, 7)
            df_full['lat'] = np.round(lats, 7)
        except Exception:
            pass

    # 类别平衡优化：全国尺度上非耕地占比高达 90%+，必须做 1:1 分层平衡采样，防止模型严重偏向非耕地
    df_crop = df_full[df_full['label'] == "Cropland"]
    df_non = df_full[df_full['label'] == "Non_Cropland"]
    
    half_target = target_sample_count // 2
    n_sample_crop = min(len(df_crop), half_target) if len(df_crop) > 0 else 0
    # 保持耕地与非耕地 1:1 平衡采样
    n_sample_non = min(len(df_non), n_sample_crop) if n_sample_crop > 0 else min(len(df_non), half_target)

    if n_sample_crop > 0 and n_sample_non > 0:
        df_crop_sampled = df_crop.sample(n=n_sample_crop, random_state=42)
        df_non_sampled = df_non.sample(n=n_sample_non, random_state=42)
        df_full = pd.concat([df_crop_sampled, df_non_sampled]).sample(frac=1.0, random_state=42).reset_index(drop=True)

    n_crop = int(np.sum(df_full['label'] == "Cropland"))
    n_non = int(np.sum(df_full['label'] == "Non_Cropland"))

    os.makedirs(os.path.dirname(out_csv_path), exist_ok=True)
    df_full.to_csv(out_csv_path, index=False)

    print("=================================================================")
    print(">>> 🎉 大图全域样本提取完成 (未进行任何切片裁剪)！")
    print(f"    • 成果文件 : {out_csv_path}")
    print(f"    • 文件大小 : {os.path.getsize(out_csv_path) / 1024 / 1024:.2f} MB")
    print(f"    • 样点总数 : {len(df_full):,} 个 (耕地: {n_crop:,} [50%], 非耕地: {n_non:,} [50%])")
    print(f"    • 全域经纬度跨度: 经度 [{df_full['lon'].min():.2f}°E ~ {df_full['lon'].max():.2f}°E], 纬度 [{df_full['lat'].min():.2f}°N ~ {df_full['lat'].max():.2f}°N]")
    print("-----------------------------------------------------------------")
    print(">>> ⚠️ 【重要空间匹配提示】：")
    print("    当前提取的是【全国/全域大尺度样点】(跨越整个中国)。")
    print("    如果下一步仅使用单个局部切片卫星瓦片 (如 50SMF) 进行训练，绝大部分全国")
    print("    样点会因超出瓦片边界而被提取算法作为 NA 丢弃 (最终仅剩几十个点)！")
    print("    • 如需针对当前单一瓦片 (50SMF) 训练高精度模型，请运行【裁剪模式】：")
    print("      python prepare_samples.py --crop --step 10")
    print("    • 如需训练全国/跨区通用大模型，请确保输入包含全国各瓦片的卫星影像。")
    print("=================================================================")
    return True


# ==============================================================================
# 主入口
# ==============================================================================
def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_sdc = os.path.join(script_dir, "data", "sdc30_cubes", "SDC30_V003_50SMF_20210618.tif")
    default_gt  = os.path.join(script_dir, "data", "SDC30_V003_50SMF_20210618_GroundTruth.tif")
    default_csv = os.path.join(script_dir, "data", "my_crop_samples.csv")

    parser = argparse.ArgumentParser(description="SDC30 作物真值切片制作与高纯度样点提取统一工具")
    
    # 核心开关：针对大图是否进行裁剪
    crop_group = parser.add_mutually_exclusive_group()
    crop_group.add_argument("--crop", action="store_true", help="【针对大图进行裁剪】：以切片范围为边界裁剪出 GroundTruth.tif，提取局部样点")
    crop_group.add_argument("--no-crop", "--full-image", action="store_true", help="【针对大图不进行裁剪】：不按切片裁剪，直接获取整张大图的全域/宏观训练样本")

    parser.add_argument("--step", type=int, default=10, help="局部切片抽样步长 (默认: 10，物理间距300米)")
    parser.add_argument("--samples", type=int, default=100000, help="全图不裁剪模式下的目标样本规模 (默认: 100000)")
    parser.add_argument("--no-erode", action="store_true", help="禁用边缘侵蚀 (保留边缘混合像元)")
    parser.add_argument("--no-features", action="store_true", help="仅提取经纬度坐标与标签，不直出 12 维遥感多光谱波段与指数")
    parser.add_argument("--crop-only", action="store_true", help="仅执行大图裁剪出切片底图，不提取样点")
    parser.add_argument("--sdc", default=default_sdc, help="SDC30 目标切片路径")
    parser.add_argument("--ref", default=None, help="来源大尺度参考真值图层路径 (如 gengdi.tif)")
    parser.add_argument("--gt", default=default_gt, help="切片真值底图保存或读取路径")
    parser.add_argument("--out-csv", default=default_csv, help="输出样点 CSV 路径")
    args = parser.parse_args()

    # ──────────────── 模式 1: 针对大图不进行裁剪 (--no-crop / --full-image) ────────────────
    if args.no_crop:
        extract_full_big_image_samples(
            ref_path=args.ref,
            out_csv_path=args.out_csv,
            target_sample_count=args.samples
        )
        return

    # ──────────────── 模式 2: 仅裁剪底图切片 (--crop-only) ────────────────
    if args.crop_only:
        crop_ground_truth_tile(sdc_path=args.sdc, ref_path=args.ref, out_gt_path=args.gt)
        print(">>> 已指定 --crop-only，切片制作完成。")
        return

    # ──────────────── 模式 3: 针对大图进行裁剪并提取局部样点 (--crop 或 智能模式) ────────────────
    need_crop = args.crop or (not os.path.exists(args.gt))
    if need_crop:
        print(">>> 正在准备切片真值底图 (从大图裁剪对齐)...")
        ok = crop_ground_truth_tile(sdc_path=args.sdc, ref_path=args.ref, out_gt_path=args.gt)
        if not ok and not os.path.exists(args.gt):
            print("[终止] 无法从大图裁剪切片，处理结束。")
            return
    else:
        print(f">>> [快速通道] 检测到本地已存在现成切片底图 ({os.path.basename(args.gt)})，直接秒级提取局部切片样点。")

    # 提取局部切片样点 (默认直出 12 维特征工程波段与指数)
    extract_local_tile_samples(
        gt_path=args.gt,
        out_csv_path=args.out_csv,
        step_size=args.step,
        erode_edges=not args.no_erode,
        sdc_path=args.sdc,
        extract_features=not args.no_features
    )


if __name__ == "__main__":
    main()
