farmcyanodir=farm_cyanolichens
farmcyanowork=farm_cyanolichen_analysis/samples

for SAMPLE in ${farmcyanodir}/*; do
	sample=`basename $SAMPLE`

	mkdir -p ${farmcyanowork}/${sample}/logs
	logdir=${farmcyanowork}/${sample}/logs

	outdir=${farmcyanowork}/${sample}/binning/concoct_metamdbg

	if [ ! -d ${outdir} ]; then
		
		mkdir -p ${outdir}
		
		bamdir=${farmcyanowork}/${sample}/reads2assembly/metamdbg
		assemblydir=${SAMPLE}/assemblies/metamdbg

		bsub -M 12000 -n 1 -oo ${logdir}/concoct_metamdbg.log "bash lrmetalichen/scripts/concoct.sh -i  ${assemblydir}/*.fasta.gz -m ${bamdir}/*.bam -o ${outdir}"

	fi
	
done
