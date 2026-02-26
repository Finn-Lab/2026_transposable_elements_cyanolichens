bindir=farm_cyanolichen_analysis/dereplicated_genomes/bacteria/metawrap_all
logdir=farm_cyanolichen_analysis/logs

mkdir -p ${logdir}


outdir=farm_cyanolichen_analysis/dereplicated_genomes/bacteria/drep

bsub  -M 25000 -n 8 -oo ${logdir}/drep_feb14.log "singularity exec \
  singularity-cache/quay.io_microbiome-informatics_genomes-pipeline.drep:v2.sif dRep dereplicate -p 8 \
  ${outdir} -g ${bindir}/*.fa -pa 0.9 -sa 0.95 -nc 0.30 -cm larger --genomeInfo \
  farm_cyanolichen_analysis/dereplicated_genomes/bacteria/drep_genomeinfo.csv -comp 50 -con 5"
