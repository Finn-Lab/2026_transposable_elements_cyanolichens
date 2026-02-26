set -x

sampledir=farm_cyanolichen_analysis/cyanobacteria_analysis/
drepdir=farm_cyanolichen_analysis/cyanobacteria_analysis/drep_mags2refs/
drepbindir=farm_cyanolichen_analysis/cyanobacteria_analysis/ref_genomes
drepqc=${sampledir}/cyanomag_ref_genomeinfo.csv
logdir=farm_cyanolichen_analysis/logs

mkdir -p ${drepbindir}


mkdir -p ${logdir}


bsub  -M 25000 -n 8 -oo ${logdir}/drep_cyanorefs_april9.log "singularity exec \
  singularity-cache/quay.io_microbiome-informatics_genomes-pipeline.drep:v2.sif dRep dereplicate -p 8 \
  ${drepdir} -g ${drepbindir}/*.fa -pa 0.9 -sa 0.95 -nc 0.30 -cm larger --genomeInfo \
  ${drepqc} -comp 50 -con 5"
