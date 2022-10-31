#!/usr/bin/env bash
#
#SBATCH -J msprime # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 1-00:00 # 1 days
#SBATCH --mem 50G
#SBATCH -o /project/berglandlab/connor/slim_bottleneck/err/neutralVCF.%A_%a.out # Standard output
#SBATCH -e /project/berglandlab/connor/slim_bottleneck/err/neutralVCF.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

# Run using:
# sbatch --array=0-5 run.msprime.sh

# Modules to load
module load anaconda/2020.11-py3.8
source activate msprime_env

# Start message
echo "Start"
date

# Working directory
wd="/project/berglandlab/connor/slim_bottleneck"

# Population size array
popList=( 1500000 2000000 3000000 5000000 7500000 10000000 )
popListName=( Neutral_1.5m Neutral_2m Neutral_3m Neutral_5m Neutral_7.5m Neutral_10m  )

# Run python script to generate Neutral VCF
python ${wd}/3.run.msprime.py \
${popList[${SLURM_ARRAY_TASK_ID}]} \
${popListName[${SLURM_ARRAY_TASK_ID}]}

# Extract summary statistics
module load vcftools

# vcftools for theta pi
vcftools --vcf ${wd}/${popListName[${SLURM_ARRAY_TASK_ID}]}.vcf \
--window-pi 10000 \
--window-pi-step 5000 \
--out ${wd}/${popListName[${SLURM_ARRAY_TASK_ID}]}

# Finish message
echo "Finish"
date
