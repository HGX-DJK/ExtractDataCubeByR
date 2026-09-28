@echo off
chcp 65001 >nul
echo =======================================================
echo 正在安装 R 遥感/时序数据立方体所需依赖包...
echo =======================================================

set "RSCRIPT=C:\Program Files\R\R-4.6.1\bin\Rscript.exe"

if not exist "%RSCRIPT%" (
    echo [错误] 未在 C:\Program Files\R\R-4.6.1\bin\Rscript.exe 找到 Rscript。
    pause
    exit /b 1
)

"%RSCRIPT%" install_dependencies.R
pause
