@echo off
chcp 65001 >nul
echo =======================================================
echo 正在运行 SDC30 数据立方体作物识别脚本...
echo =======================================================

set "RSCRIPT=C:\Program Files\R\R-4.6.1\bin\Rscript.exe"

if not exist "%RSCRIPT%" (
    echo [错误] 未在 C:\Program Files\R\R-4.6.1\bin\Rscript.exe 找到 Rscript。
    echo 请确认 R 的安装路径。
    pause
    exit /b 1
)

"%RSCRIPT%" 05_pcl_sdc30_crop_classification.R
pause
