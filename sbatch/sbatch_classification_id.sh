#!/bin/sh

#SBATCH -p long
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --output=/home/enrico.saccon/mpdp_new/sbatch/classification_id/output_classification_id_%j.txt
#SBATCH --error=/home/enrico.saccon/mpdp_new/sbatch/classification_id/error_classification_id_%j.txt
#SBATCH --time=2-00:00:00

module load CUDA/12.1.1

date=$(date +%Y/%m/%d_%H:%M:%S)
echo "${date}"
echo "${date}" 1>&2;

source /home/enrico.saccon/mpdp/venv/bin/activate
cd /home/enrico.saccon/mpdp_new/examples/3PMD/prediction/rectangle/NN/Classification
python3 main.py --yaml-config /home/enrico.saccon/mpdp_new/examples/3PMD/prediction/rectangle/NN/Classification/classification_config.yaml