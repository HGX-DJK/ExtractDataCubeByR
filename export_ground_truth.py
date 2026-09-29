# -*- coding: utf-8 -*-
"""
从 42GB 全球/全国参考图层 gengdi.tif 中精确裁剪并重投影出当前 SDC30 切片的全量真值底图
不抽样、100%全像元逐像元几何空间对齐
"""

import os
import sys
import io
import rasterio
from rasterio.warp import reproject, Resampling, transform_bounds
from rasterio.windows import from_bounds
import numpy as np

# 兼容 Windows 控制台输出编码
if sys.platform.startswith("win"):
    try:
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
        sys.stderr.reconfigure(encoding='utf-8', errors='replace')
    except Exception:
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
        sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

def export_full_ground_truth():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    sdc_path = os.path.join(script_dir, "data", "sdc30_cubes", "SDC30_V003_50SMF_20210618.tif")
    ref_path = r"E:\agriculture\gaced30_validation_pipeline\data\reference_data\gengdi.tif"
    out_dir = os.path.join(script_dir, "data")
    out_path = os.path.join(out_dir, "SDC30_V003_50SMF_20210618_GroundTruth.tif")

    if not os.path.exists(sdc_path):
        # 尝试备用路径
        sdc_path = os.path.join(script_dir, "..", "data", "sdc30_cubes", "SDC30_V003_50SMF_20210618.tif")

    if not os.path.exists(ref_path):
        print(f"[错误] 未找到参考真值底图: {ref_path}")
        return

    print("=================================================================")
    print(">>> 正在从 42GB gengdi.tif 中无损提取 100% 全像元严格对齐的真实底图 (不抽样)...")
    print(f"    • 目标数据立方体: {os.path.basename(sdc_path)}")
    print(f"    • 参考来源图层  : {ref_path}")
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

    # 规范化二值编码：1=耕地, 0=非耕地, 255=NoData
    dst_binary = np.where((dst_data == 1) | (dst_data == 10), np.uint8(1), np.where(dst_data == 0, np.uint8(0), np.uint8(255)))

    meta.update({
        'count': 1,
        'dtype': 'uint8',
        'nodata': 255,
        'compress': 'lzw'
    })

    os.makedirs(out_dir, exist_ok=True)
    with rasterio.open(out_path, 'w', **meta) as dst:
        dst.write(dst_binary, 1)

    crop_count = int(np.sum(dst_binary == 1))
    noncrop_count = int(np.sum(dst_binary == 0))
    nodata_count = int(np.sum(dst_binary == 255))
    total_valid = crop_count + noncrop_count

    print(">>> 导出成功！")
    print(f"    • 成果位置: {out_path}")
    print(f"    • 文件大小: {os.path.getsize(out_path) / 1024 / 1024:.2f} MB")
    print(f"    • 影像尺寸: {sdc_shape[0]} × {sdc_shape[1]} (共 {dst_binary.size:,} 像元)")
    print(f"    • 空间分辨率: 30 米 | 坐标系: {sdc_crs}")
    print(f"    • 统计分布 (100%全覆盖普查，零抽样):")
    print(f"        - 真实耕地像元数 (1): {crop_count:,} ({crop_count / total_valid * 100:.2f}%)")
    print(f"        - 真实非耕地像元 (0): {noncrop_count:,} ({noncrop_count / total_valid * 100:.2f}%)")
    if nodata_count > 0:
        print(f"        - 背景无效像元 (255): {nodata_count:,}")
    print("=================================================================")

if __name__ == "__main__":
    export_full_ground_truth()
