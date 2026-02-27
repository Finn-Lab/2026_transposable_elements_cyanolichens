workdir=cyanolichen_rnaseq_analysis
logdir=${workdir}/logs/interproscan
mkdir -p ${logdir}

interpro=${workdir}/interproscan_trebouxoid

mkdir -p ${interpro}
transdecoder=${workdir}/transdecoder_trebouxoid
	
queryfile=${transdecoder}/*pep  

sed 's/*//g' ${queryfile} > ${transdecoder}/Trebouxoid_symbiont.combined.interproquery.faa

interproquery=${transdecoder}/Trebouxoid_symbiont.combined.interproquery.faa
	
outdir=${interpro}/
mkdir -p ${outdir}

interproscan.sh -i ${interproquery} -appl Pfam,TIGRFAM -cpu 10 -d ${outdir} -iprlookup -pa -goterms

