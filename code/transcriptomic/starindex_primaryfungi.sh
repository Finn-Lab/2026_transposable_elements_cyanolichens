#!/bin/bash

#Submit this script with: sbatch thefilename
#For more details about each parameter, please check SLURM sbatch documentation https://slurm.schedmd.com/sbatch.html

#SBATCH --time=12:00:00
#SBATCH --ntasks=1   # number of tasks
#SBATCH --cpus-per-task=12   # number of CPUs Per Task i.e if your code is multi-threaded
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=24G   # memory per node
#SBATCH -J "star_index"   # job name
#SBATCH -o "cyanolichen_rnaseq_analysis/logs/star_primaryfungi_index_v2.log"   # job output file
#SBATCH --array=1-15%1 # job array number

# Specify path to the config file
rnaseq_config=/cyanolichen_rnaseq_analysis/rnaseq_config_11.txt

# Extract sample name for array task ID
magfasta=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $rnaseq_config)
magid=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $rnaseq_config)
species=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $rnaseq_config)

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
magfastapath=${magfasta}

metaeukdir=farm_cyanolichen_analysis/metamdbg_fungi_analysis/metaeuk/
metaeukgff=${metaeukgff}/${magid}/${magid}.gff

indexdir=cyanolichen_rnaseq_analysis/star_index_fungi/${species}/${magid}

if [ ! -d ${indexdir} ]; then

	mkdir -p ${indexdir}

	STAR \
		--runThreadN 12 --runMode genomeGenerate --genomeDir ${indexdir} \
		--genomeSAindexNbases 11 \
		--genomeFastaFiles ${magfastapath} \
		--sjdbGTFtagExonParentTranscript ${metaeukgff}
fi
