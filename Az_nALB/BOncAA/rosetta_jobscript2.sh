#!/bin/bash
#SBATCH --account=PAS2252
#SBATCH --job-name=random
#SBATCH --time=90:00:00
#SBATCH --nodes=1
#SBATCH --mem=150G
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=ALL
#SBATCH --mail-user=phan.202@buckeyemail.osu.edu      # <-- add this if you want notifications
#SBATCH --output=relax5pleasework.txt
#SBATCH --error=relax5errorpleasework.txt

# Load required modules
module load rosetta/3.12

# Activate conda & necessary environments
source /users/PAS2252/nathankvphan4/miniconda3/etc/profile.d/conda.sh
conda activate ncaa

# Go to working directory
cd /users/PAS2252/nathankvphan4/Az_nALB/BOncAA

echo "Starting at $(date)"

# Run parameterization with srun
srun python random_sele.py

echo "Complete at $(date)"
