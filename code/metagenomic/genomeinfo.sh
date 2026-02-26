set -x

newdir=farm_cyanolichen_analysis/dereplicated_genomes/bacteria/

for SAMPLE in farm_cyanolichen_analysis/samples/*; do

	sample=`basename $SAMPLE`

	refined_bins=${SAMPLE}/refined_bins/hifiasm_binref/metawrap_50_5_bins.stats

    awk -v ID=${sample} 'BEGIN { FS="\t"; OFS="," } NR==1{ print $1,$2,$3 } NR>1{ print ID "_" $1 ".fa",$2,$3 }' ${refined_bins} > ${newdir}/${sample}_qc.csv

done

awk '(NR == 1) || (FNR > 1)' ${newdir}/*qc.csv > ${newdir}/drep_qc.csv

awk 'BEGIN { FS=","; OFS="," } NR==1{ print "genome",$2,$3 } NR>1{ print $0 }' ${newdir}/drep_qc.csv > ${newdir}/drep_genomeinfo.csv

rm ${newdir}/*qc.csv


