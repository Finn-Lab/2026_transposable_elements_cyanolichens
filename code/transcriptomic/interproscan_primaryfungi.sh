workdir=cyanolichen_rnaseq_analysis
logdir=${workdir}/logs/interproscan
mkdir -p ${logdir}

interpro=${workdir}/interproscan_primaryfungi

mkdir -p ${interpro}
transdecoder=${workdir}/transdecoder

for SAMPLE in ${transdecoder}/*; do
	sample=`basename $SAMPLE`
	
	queryfile=${SAMPLE}/*pep  

	sed 's/*//g' ${queryfile} > ${SAMPLE}/${sample}.interproquery.faa

	interproquery=${SAMPLE}/${sample}.interproquery.faa
	
	outdir=${interpro}/${sample}
	mkdir -p ${outdir}

	sbatch --time=72:00:00 --mem 48000 -c 10  -o ${logdir}/${sample}.log --wrap="interproscan.sh -i ${interproquery} -appl Pfam,TIGRFAM \
		-cpu 10 -d ${outdir} -iprlookup -pa -goterms"

done
