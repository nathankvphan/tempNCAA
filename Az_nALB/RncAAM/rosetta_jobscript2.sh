#!/bin/bash
#SBATCH --account=PAS2252
#SBATCH --job-name=MEGARun5
#SBATCH --time=90:00:00
#SBATCH --nodes=3
#SBATCH --mem=100
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=ALL
#SBATCH --mail-user=phan.202@buckeyemail.osu.edu      # <-- add this if you want notifications
#SBATCH --output=rosetta_param_output.txt
#SBATCH --error=rosetta_param_error.txt

# Load required modules
module load rosetta/3.12

# Activate conda & necessary environments
source /users/PAS2252/nathankvphan4/miniconda3/etc/profile.d/conda.sh
conda activate ncaa

# Go to working directory
cd /users/PAS2252/nathankvphan4/Projects/ncPPI/RncAAM

echo "Starting at $(date)"

# Run parameterization with srun
srun python example.py

echo "Complete at $(date)"
