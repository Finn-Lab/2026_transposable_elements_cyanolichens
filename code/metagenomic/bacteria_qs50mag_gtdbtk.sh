set -ex

bindir=farm_cyanolichen_analysis/bacteria_qs50

taxdir=farm_cyanolichen_analysis/qs50_bacteria_gtdbtk

logdir=farm_cyanolichen_analysis/logs

export GENOMEDB=ref-dbs/genomes-pipeline  

sbatch --time=1-0 --mem=100G -c 2 -J gtdbtk -o ${logdir}/qs50_bacteria_gtdbtk.log --wrap="singularity run --bind \
  $GENOMEDB/release214:/refdata \
  singularity-cache/quay.io_microbiome-informatics_genomes-pipeline.gtdb-tk:v2.3.0.sif gtdbtk classify_wf --cpus 2 \
  --genome_dir ${bindir} --out_dir ${taxdir} -x fa --skip_ani_screen"
