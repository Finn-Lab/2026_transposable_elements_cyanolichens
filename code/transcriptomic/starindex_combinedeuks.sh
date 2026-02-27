#!/bin/bash

#Submit this script with: sbatch thefilename
#For more details about each parameter, please check SLURM sbatch documentation https://slurm.schedmd.com/sbatch.html

#SBATCH --time=3-0
#SBATCH --ntasks=1   # number of tasks
#SBATCH --cpus-per-task=12   # number of CPUs Per Task i.e if your code is multi-threaded
#SBATCH --nodes=15   # number of nodes
#SBATCH --mem=24G   # memory per node
#SBATCH -J "star_index"   # job name
#SBATCH -o "cyanolichen_rnaseq_analysis/logs/starindex_combined_%A_%a.log"   # job output file
#SBATCH --array=1-2 # job array number

# Specify path to the config file
rnaseq_config=cyanolichen_rnaseq_analysis/rnaseq_config_combofun2.txt

# Extract sample name for array task ID
magfasta=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $rnaseq_config)
magid=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $rnaseq_config)
species=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $rnaseq_config)
tolid=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $6}' $rnaseq_config)

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
magfastapath=cyanolichen_rnaseq_analysis/refs/combined_fasta/${tolid}.fa

metaeukdir=cyanolichen_rnaseq_analysis/refs/combined_gff
metaeukgff=${metaeukgff}/${tolid}.gff

indexdir=cyanolichen_rnaseq_analysis/star_index_combination/${species}/${tolid}

if [ ! -d ${indexdir} ]; then

	mkdir -p ${indexdir}

	STAR \
		--runThreadN 12 --runMode genomeGenerate --genomeDir ${indexdir} \
		--genomeSAindexNbases 11 \
		--genomeFastaFiles ${magfastapath} \
		--sjdbGTFtagExonParentTranscript ${metaeukgff}
fi
