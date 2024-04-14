#!/bin/bash
echo "Configuring and building C++ program with CMake..."
mkdir -p ../build && cd ../build  # 创建并进入构建目录
cmake ..  # 生成 Makefile
make      # 构建项目

echo "Running C++ program..."
../bin/data_generator  # 执行生成的可执行文件

echo "Running MATLAB script..."
matlab -batch "run('../matlab/plot_data.m')"

echo "Done."

