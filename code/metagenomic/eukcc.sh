cyanodir=farm_cyanolichen_analysis/samples

for SAMPLE in ${cyanodir}/*
do

  sample=`basename $SAMPLE`
  logdir=${SAMPLE}/logs

  bindir=${SAMPLE}/binning/metabat2_metamdbg/bins
  outdir=${SAMPLE}/eukcc
	
  if [ $(grep -c "TERM_MEM" ${logdir}/eukcc.log) -gt 0 ]; then
  	rm -rf ${outdir}
  	bsub -M 24000 -n 8 -oo ${logdir}/eukcc.log "singularity exec  singularity-cache/quay.io-microbiome-informatics-eukcc-latest.img \
 	 eukcc folder --db databases/eukcc2/v4_iqtree --out ${outdir} --threads 8 ${bindir}"
  fi

done
