#!/bin/bash -l
#SBATCH --job-name=salmon-array
#SBATCH --array=1-60
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=1-05:005:00
#SBATCH --partition brc
#SBATCH --mem-per-cpu=6000MB
#SBATCH --mail-type=END,FAIL
# Autogenerated script from salmon_array_job.R
# date Wed Sep 16 11:47:50 2020
# make sure directory paths exist before running script



module load apps/salmon/1.2.1-singularity



# Parse parameter file to get variables.
number=$SLURM_ARRAY_TASK_ID
paramfile=salmon_param.txt
 
inr1=`sed -n ${number}p $paramfile | awk '{print $1}'`
inr2=`sed -n ${number}p $paramfile | awk '{print $2}'`
outsal=`sed -n ${number}p $paramfile | awk '{print $3}'`

# 9. Run the program.
salmon quant -i /users/k1625253/scratch/old-scratch_rosalind-legacy-import_2020-01-28/Data/MetaData/GenomeIndex/salmon_transcriptome_index/gencode.v35.transcripts_index/ -p 8 -l ISR -1 $inr1 -2 $inr2 -o $outsal


