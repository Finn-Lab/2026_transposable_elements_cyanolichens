farmcyanodir=farm_cyanolichens
farmcyanowork=farm_cyanolichen_analysis
farmcyanoworksamples=farm_cyanolichen_analysis/samples

for SAMPLE in ${farmcyanodir}/*; do
	sample=`basename $SAMPLE`

	mkdir -p ${farmcyanoworksamples}/${sample}/logs
	logdir=${farmcyanoworksamples}/${sample}/logs

	outdir=${farmcyanoworksamples}/${sample}/qs50_mag_map
	
	if [ ! -d ${outdir} ]; then
	
		mkdir -p ${outdir}
		
		
		readsdir=${SAMPLE}/qc_reads
		ref=${farmcyanowork}/qs50_mags_all.fasta

		sbatch --time=2-0 --mem=16000 -c 1 -J lrmap2ref -o ${logdir}/lrmap2qs50_metamdbg.log --wrap="bash lrmetalichen/scripts/lrmap2ref.sh -t 1 -i ${readsdir}/* -r ${ref} -o ${outdir}/${sample} -c contigs"
	fi
	
done
