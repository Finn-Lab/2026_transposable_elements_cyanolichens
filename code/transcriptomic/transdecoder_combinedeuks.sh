#!/bin/bash

#Submit this script with: sbatch thefilename
#For more details about each parameter, please check SLURM sbatch documentation https://slurm.schedmd.com/sbatch.html

#SBATCH --time=3-0
#SBATCH --ntasks=1   # number of tasks
#SBATCH --cpus-per-task=24   # number of CPUs Per Task i.e if your code is multi-threaded
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=64G   # memory per node
#SBATCH -J "transdecoder"   # job name
#SBATCH -o "logs/transdecoderpredict_homology1_combo_%A_%a.log"   # job output file
#SBATCH --array=1-2 # job array number

# Specify path to the config file
rnaseq_config=cyanolichen_rnaseq_analysis/rnaseq_config_combofun2.txt

# Extract sample name for array task ID
magfasta=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $rnaseq_config)
magid=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $6}' $rnaseq_config)
species=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $rnaseq_config)
stringprefix=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $6}' $rnaseq_config)

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
transdecoder=TransDecoder-TransDecoder-v5.7.1
transdecoder_singularity=containers/transdecoder.v5.7.1.simg

stringdir=cyanolichen_rnaseq_analysis/stringtie_combo_final_v2/
stringgtf=${stringdir}/${stringprefix}.combined.gtf
stringtxfa=${stringdir}/${stringprefix}.tx.fa
genemap=cyanolichen_rnaseq_analysis/stringtie_combo_final_v2/${stringprefix}.gene_to_tx.txt

outdir=cyanolichen_rnaseq_analysis/transdecoder_combo/${species}
mkdir -p ${outdir}

uniref90=databases/uniref90.fasta
diamonddb=databases/uniref90_dmnd.dmnd

# construct transcript fasta 
#gtf_to_alignment_gff3.pl ${stringgtf} > ${outdir}/transcripts.gff3

#TransDecoder.LongOrfs -t ${stringtxfa} -S --gene_trans_map ${genemap} --output_dir ${outdir}


mkdir -p ${outdir}/homologysearch
homologydir=${outdir}/homologysearch

#blastp homology search
#blastp -query ${outdir}/*dec*/longest_orfs.pep \
#	-db ${uniref90} -max_target_seqs 1 \
#	-outfmt 6 -evalue 1e-5 -num_threads 24 > ${homologydir}/blastp_uniref90.outfmt6

#diamond blastp -d ${diamonddb} -q ${outdir}/*dec*/longest_orfs.pep \
#	-o ${homologydir}/diamondblastp_uniref90.outfmt6 \
#	-f 6 -k 1 -e 1e-5 -p 24
	
# hmmsearch pfam homology
#hmmsearch --cpu 16 -E 1e-10 --domtblout ${homologydir}/pfam.domtblout \
#	databases/Pfam-A.hmm \
#	${outdir}/*deco*/longest_orfs.pep


#speciesdir=cyanolichen_rnaseq/samples/${species}
#readsdir=${speciesdir}/rna_qc/fastq

TransDecoder.Predict -t ${stringtxfa} --retain_pfam_hits ${homologydir}/pfam.domtblout --retain_blastp_hits ${homologydir}/diamondblastp_uniref90.outfmt6 \
	--output_dir ${outdir}

