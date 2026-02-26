set -ex
workdir=farm_cyanolichen_analysis/dereplicated_genomes_metamdbg2bin/bacteria
bindir=${workdir}/dereplicated_genomes
prokkaout=${workdir}/prokka
logdir=farm_cyanolichen_analysis/logs/prokka_metamdbg
mkdir -p ${logdir}

for SAMPLE in ${bindir}/*

do
  sample=`basename $SAMPLE`
  sample=${sample%.fa}

  bsub -M 16000 -n 1 -oo ${logdir}/${sample}_prokka.log \
  "prokka --outdir ${prokkaout}/${sample} --prefix ${sample} ${SAMPLE}"

done