#!/bin/sh

#SBATCH -p long
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH --gres gpu:1
#SBATCH --output output_classification.txt
#SBATCH --error error_classification.txt

module load CUDA/12.1.1

source /home/enrico.saccon/mpdp/venv/bin/activate
cd /home/enrico.saccon/mpdp/examples/3PMD/prediction/rectangle/NN/Classification
python3 main.py
