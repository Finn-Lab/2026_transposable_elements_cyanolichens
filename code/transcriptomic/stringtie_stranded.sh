#!/bin/bash

#Submit this script with: sbatch thefilename
#For more details about each parameter, please check SLURM sbatch documentation https://slurm.schedmd.com/sbatch.html

#SBATCH --time=3-0
#SBATCH --ntasks=1   # number of tasks
#SBATCH --cpus-per-task=12   # number of CPUs Per Task i.e if your code is multi-threaded
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=64G   # memory per node
#SBATCH -J "string"   # job name
#SBATCH -o "cyanolichen_rnaseq_analysis/logs/stringtiestranded_primaryfungi_v2_%A_%a.log"   # job output file
#SBATCH --array=1-2 # job array number

# Specify path to the config file
rnaseq_config=cyanolichen_rnaseq_analysis/rnaseq_config_combofun2.txt

# Extract sample name for array task ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $rnaseq_config)
magfasta=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $rnaseq_config)
magid=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $7}' $rnaseq_config)
species=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $rnaseq_config)

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
speciesdir=cyanolichen_rnaseq/samples/${species}

outdir=cyanolichen_rnaseq_analysis/star_reads_combination/${species}
indexdir=cyanolichen_rnaseq_analysis/star_index_combination/${species}/${magid}

stringdir=cyanolichen_rnaseq_analysis/stringtie_combo/${species}
mkdir -p ${stringdir}

mergedbam=${stringdir}/${species}.merged.bam

samtools view -b -f 128 -F 16 ${mergedbam} > ${stringdir}/${species}.fwd1.bam
samtools view -b -f 64 -F 32 ${mergedbam} > ${stringdir}/${species}.fwd2.bam
samtools merge -o ${stringdir}/${species}.fwdmerge.bam ${stringdir}/*.fwd*.bam -@ 12

samtools view -b -f 144 ${mergedbam} > ${stringdir}/${species}.rev1.bam
samtools view -b -f 96 ${mergedbam} > ${stringdir}/${species}.rev2.bam
samtools merge -o ${stringdir}/${species}.revmerge.bam ${stringdir}/*.rev*.bam -@ 12

stringtie ${stringdir}/${species}.revmerge.bam -o ${stringdir}/${species}.revstringtie_v2.gtf --rf -l ${species}.STRGR -A ${stringdir}/${species}.revstringtie_gene_abund_v2.txt -p 12 -j 5 -g 5

stringtie ${stringdir}/${species}.fwdmerge.bam -o ${stringdir}/${species}.fwdstringtie_v2.gtf --rf -l ${species}.STRGF -A ${stringdir}/${species}.fwdstringtie_gene_v2abund.txt -p 12 -j 5 -g 5

