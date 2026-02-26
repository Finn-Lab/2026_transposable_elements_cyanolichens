farmcyanodir=farm_cyanolichens
farmcyanowork=farm_cyanolichen_analysis/samples

for SAMPLE in ${farmcyanodir}/*; do
	sample=`basename $SAMPLE`

	mkdir -p ${farmcyanowork}/${sample}/logs
	logdir=${farmcyanowork}/${sample}/logs

	outdir=${farmcyanowork}/${sample}/binning/metabat2_metamdbg

	if [ ! -d ${outdir} ]; then
		
		mkdir -p ${outdir}
		
		bamdir=${farmcyanowork}/${sample}/reads2assembly/metamdbg
		assemblydir=${SAMPLE}/assemblies/metamdbg

		bsub -M 8000 -n 1 -oo ${logdir}/metabat2_metamdbg.log "bash lrmetalichen/scripts/metabat_bin.sh -i  ${assemblydir}/*.fasta.gz -m ${bamdir}/*.bam -o ${outdir}"
	fi
	
done
