bindir=farm_cyanolichen_analysis/samples

for SAMPLE in ${bindir}/*
do
  sample=`basename $SAMPLE`
    if [ ! -d ${SAMPLE}/refined_bins/metamdbg_binref_2bin ]; then   
      mkdir -p ${SAMPLE}/refined_bins/metamdbg_binref_2bin
      outdir=${SAMPLE}/refined_bins/metamdbg_binref_2bin

      logdir=${SAMPLE}/logs
      
      binningdir=${SAMPLE}/binning

      metabat2=${binningdir}/metabat2_metamdbg/bins
      concoct=${binningdir}/concoct_metamdbg/concoct_output/fasta_bins
	  #lrbinner=${binningdir}/lrbinner_metamdbg/binned_contigs
      
      bsub -M 40000 -n 8 -oo ${logdir}/metamdbg_binrefinement_2bin.log "metawrap bin_refinement -o ${outdir} -t 8 -A ${metabat2} -B ${concoct} -c 50 -x 5"

    fi

done
