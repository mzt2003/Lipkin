#!/bin/bash
echo "Compiling C++ program..."
g++ ../src/main.cpp -o ../bin/data_generator

echo "Running C++ program..."
../bin/data_generator

echo "Running MATLAB script..."
matlab -batch "run('../matlab/plot_data.m')"

echo "Done."
chmod +x scripts/run_process.sh
