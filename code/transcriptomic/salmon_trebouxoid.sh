#!/bin/bash

#Submit this script with: sbatch thefilename
#For more details about each parameter, please check SLURM sbatch documentation https://slurm.schedmd.com/sbatch.html

#SBATCH --time=3-0
#SBATCH --ntasks=1   # number of tasks
#SBATCH --cpus-per-task=12   # number of CPUs Per Task i.e if your code is multi-threaded
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=64G   # memory per node
#SBATCH -J "salmon"   # job name
#SBATCH -o "cyanolichen_rnaseq_analysis/logs/salmon_trebouxoid_%A_%a.log"   # job output file
#SBATCH --array=1-4%1 # job array number

# Specify path to the config file
rnaseq_config=rnaseq_config_trebouxoid.txt

# Extract sample name for array task ID
#magfasta=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $rnaseq_config)
magid=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $rnaseq_config)
species=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $rnaseq_config)
stringprefix=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $7}' $rnaseq_config)

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
#salmon=bin/salmon

speciesdir=cyanolichen_rnaseq/samples/${species}
readsdir=${speciesdir}/rna_qc/fastq

stringgtf=cyanolichen_rnaseq_analysis/combined_trebouxoid/gffcmp.combined.gtf
stringtxfa=cyanolichen_rnaseq_analysis/combined_trebouxoid/Trebouxoid_symbiont.combined.tx.fa
magfasta=cyanolichen_rnaseq_analysis/${magid}

salmondir=cyanolichen_rnaseq_analysis/salmon_trebouxoid/${species}
mkdir -p ${salmondir}


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
