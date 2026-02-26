lichens=farm_cyanolichens

for SAMPLE in ${lichens}/*; do
	sample=`basename $SAMPLE`

	if [ ! -d ${SAMPLE}/assemblies/metamdbg ]; then
		mkdir -p ${SAMPLE}/assemblies/metamdbg

		bsub -M 12000 -n 2 -oo ${SAMPLE}/assemblies/metamdbg/${sample}_metamdbg.log "metaMDBG asm ${SAMPLE}/assemblies/metamdbg ${SAMPLE}/qc_reads/*.fasta.gz"
	fi
done
