set -ex

bindir=farm_cyanolichen_analysis/cyanobacteria_analysis/drep_mags2refs/dereplicated_genomes

taxdir=farm_cyanolichen_analysis/cyanobacteria_analysis/drep_gtdbtk

logdir=farm_cyanolichen_analysis/logs

export GENOMEDB=ref-dbs

bsub -M 100000 -n 2 -oo ${logdir}/cyanomagsdrep_gtdbtk.log "singularity run --bind \
  $GENOMEDB/release214:/refdata \
  singularity-cache/quay.io_microbiome-informatics_genomes-pipeline.gtdb-tk:v2.3.0.sif gtdbtk classify_wf --cpus 2 \
  --genome_dir ${bindir} --out_dir ${taxdir} -x fa --skip_ani_screen"
