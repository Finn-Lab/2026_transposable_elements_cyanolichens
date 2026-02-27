#!/bin/bash

#Submit this script with: sbatch thefilename
#For more details about each parameter, please check SLURM sbatch documentation https://slurm.schedmd.com/sbatch.html

#SBATCH --time=3-0
#SBATCH --ntasks=1   # number of tasks
#SBATCH --cpus-per-task=12   # number of CPUs Per Task i.e if your code is multi-threaded
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=64G   # memory per node
#SBATCH -J "salmon"   # job name
#SBATCH -o "logs/salmon_primaryfungirep_%A_%a.log"   # job output file
#SBATCH --array=1-9 # job array number

# Specify path to the config file
rnaseq_config=cyanolichen_rnaseq_analysis/rnaseq_config_primaryfun.txt

# Extract sample name for array task ID
magfasta=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $rnaseq_config)
magid=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $rnaseq_config)
species=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $rnaseq_config)
stringprefix=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $rnaseq_config)

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
salmon=bin/salmon

speciesdir=cyanolichen_rnaseq/samples/${species}
readsdir=${speciesdir}/rna_qc/fastq

outdir=cyanolichen_rnaseq_analysis/star_reads_fungi/${species}
indexdir=cyanolichen_rnaseq_analysis/star_index_fungi/${species}/${magid}

stringgtf=cyanolichen_rnaseq/stringtie_final_primaryfungi/${stringprefix}.combined.gtf
stringtxfa=cyanolichen_rnaseq/stringtie_final_primaryfungi/${stringprefix}.tx.fa

stringdir=cyanolichen_rnaseq_analysis/stringtie_fungi/${species}
mkdir -p ${stringdir}

salmondir=cyanolichen_rnaseq_analysis/salmon_fungi/${species}
mkdir -p ${salmondir}

mergedbam=${stringdir}/${species}.merged.bam

#generate decoys
bash SalmonTools/scripts/generateDecoyTranscriptome.sh \
	-a ${stringgtf} -g ${magfasta} -t ${stringtxfa} -o ${salmondir} -j 12
	
# index transcripts
salmon index -t ${salmondir}/gentrome.fa -i ${salmondir}/transcripts_index --decoys ${salmondir}/decoys.txt -k 31

for SAMPLE in ${readsdir}/*_1.fastq; do
	sample=`basename $SAMPLE`
	sample=${sample%_trimmed_1.fastq}

	read1=${readsdir}/${sample}_trimmed_1.fastq
	read2=${readsdir}/${sample}_trimmed_2.fastq
	
	salmon quant -i ${salmondir}/transcripts_index -l ISR -1 ${read1} -2 ${read2} \
		--validateMappings -p 12 --gcBias \
		-g ${stringgtf} -o ${salmondir}/transcripts_quant/${sample}
done
