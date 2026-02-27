#!/bin/bash

#Submit this script with: sbatch thefilename
#For more details about each parameter, please check SLURM sbatch documentation https://slurm.schedmd.com/sbatch.html

#SBATCH --time=3-0
#SBATCH --ntasks=1   # number of tasks
#SBATCH --cpus-per-task=4   # number of CPUs Per Task i.e if your code is multi-threaded
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=24G   # memory per node
#SBATCH -J "star"   # job name
#SBATCH -o "logs/starstring_combined_array_%A_%a.log"   # job output file
#SBATCH --array=4 # job array number

# Specify path to the config file
rnaseq_config=cyanolichen_rnaseq_analysis/rnaseq_config_combo2.txt

# Extract sample name for array task ID
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $rnaseq_config)
magfasta=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $rnaseq_config)
magid=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $6}' $rnaseq_config)
species=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $rnaseq_config)
tolid=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $7}' $rnaseq_config)

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
speciesdir=cyanolichen_rnaseq/samples/${species}
readdir=${speciesdir}/rna_qc/fastq

read1=${readdir}/${sample}_trimmed_1.fastq
read2=${readdir}/${sample}_trimmed_2.fastq

magfastapath=cyanolichen_rnaseq_analysis/${magfasta}

metaeukdir=cyanolichen_rnaseq_analysis/refs/combined_gff
metaeukgff=${metaeukgff}/${tolid}.gff

outdir=cyanolichen_rnaseq_analysis/star_reads_combination/${species}
indexdir=cyanolichen_rnaseq_analysis/star_index_combination/${species}/${tolid}

if [ ! -d ${outdir}/${sample} ]; then
	mkdir -p ${outdir}/${sample}

#STAR \
#	--runThreadN 12 --runMode genomeGenerate --genomeDir ${indexdir} \
#	--genomeSAindexNbases 11 \
#	--genomeFastaFiles ${magfastapath} \
#	--sjdbGTFtagExonParentTranscript ${metaeukgff}

	STAR \
	    --genomeDir ${indexdir} \
	    --runThreadN 12 \
	    --readFilesIn \
	      ${read1} \
	      ${read2} \
	    --outFileNamePrefix ${outdir}/${sample}/${sample}. \
	    --outSAMstrandField intronMotif \
	    --outMultimapperOrder Random \
	    --outSAMtype BAM Unsorted \
	    --outFilterMultimapNmax 10000 \
	    --limitOutSAMoneReadBytes 10000000

	samtools fixmate -m -u ${outdir}/${sample}/${sample}.Aligned.out.bam - | samtools sort -u - | samtools markdup --write-index - ${outdir}/${sample}/${sample}_fixsort.bam

	STAR \
		--runMode inputAlignmentsFromBAM \
		--inputBAMfile ${outdir}/${sample}/${sample}_fixsort.bam \
		--outWigType bedGraph \
		--outWigNorm RPM  \
		--outFileNamePrefix ${outdir}/${sample}/${sample}.
fi
