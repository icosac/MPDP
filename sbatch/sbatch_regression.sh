#!/bin/sh

#SBATCH -p long
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH -N 1
#SBATCH --output output_regression.txt
#SBATCH --error error_regression.txt

module load CUDA/12.1.1

source /home/enrico.saccon/mpdp/venv/bin/activate
cd /home/enrico.saccon/mpdp/examples/3PMD/prediction/rectangle/NN/Regression
python3 main.py
