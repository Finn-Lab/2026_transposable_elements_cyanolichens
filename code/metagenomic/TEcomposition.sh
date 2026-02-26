set -ex

algaedir=farm_cyanolichen_analysis/trebouxia_analysis
fungidir=farm_cyanolichen_analysis/metamdbg_fungi_analysis

mkdir -p ${algaedir}/te_bed_gt100_frac
mkdir -p ${fungidir}/te_bed_gt100_frac

mkdir -p ${algaedir}/chromsizes
algchroms=${algaedir}/chromsizes

mkdir -p ${fungidir}/chromsizes
funchroms=${fungidir}/chromsizes

#for SAMPLE in ${algaedir}/ncbi-genomes-2024-04-05/*.fna; do
#	samtools faidx ${SAMPLE}
#done

for SAMPLE in ${fungidir}/dothidio_refs/*.fa; do
	samtools faidx ${SAMPLE}
done

#for SAMPLE in ${algaedir}/ncbi-genomes-2024-04-05/*.fna.fai; do
#	sample=`basename $SAMPLE`
#	sample=${sample%.fna.fai}
#	awk '{print $1 "\t" $2}' ${SAMPLE} > ${algaedir}/chromsizes/${sample}.txt
#done
#
for SAMPLE in ${fungidir}/dothidio_refs/*.fa.fai; do
	sample=`basename $SAMPLE`
	sample=${sample%.fa.fai}
	awk '{print $1 "\t" $2}' ${SAMPLE} > ${fungidir}/chromsizes/${sample}.txt
done

algaete=${algaedir}/earlgrey_refs
fungite=${fungidir}/earlgrey_ncbirefs

mkdir -p ${algaedir}/te_bed_gt100
mkdir -p ${fungidir}/te_bed_gt100

#for SAMPLE in ${algaete}/*; do
#	sample=`basename ${SAMPLE}`
#	sample=${sample%_EarlGrey}
#	bed=${SAMPLE}/${sample}_summaryFiles/${sample}.filteredRepeats.bed
#	if [ -f ${bed} ]; then
#		awk '{if (($3-$2) >= 100) {print $0}}' ${bed} > ${algaedir}/te_bed_gt100/${sample}.filteredRepeats.bed
#	fi
#done

for SAMPLE in ${fungite}/*; do
	sample=`basename ${SAMPLE}`	
	sample=${sample%_EarlGrey}
	bed=${SAMPLE}/${sample}_summaryFiles/${sample}.filteredRepeats.bed

	if [ -f ${bed} ]; then
		awk '{if (($3-$2) >= 100) {print $0}}' ${bed} > ${fungidir}/te_bed_gt100/${sample}.filteredRepeats.bed
	fi
done

#for i in ${algaedir}/te_bed_gt100/* ; do
#  SAMPLE=`basename $i`
#  SAMPLE=${SAMPLE%.filteredRepeats.bed}
#  bedtools genomecov -max 1 -g ${algaedir}/chromsizes/${SAMPLE}.txt \
#    -i ${i} \
#    > ${algaedir}/te_bed_gt100_frac/${SAMPLE}.txt
#done
		
for i in ${fungidir}/te_bed_gt100/* ; do
  SAMPLE=`basename $i`
  SAMPLE=${SAMPLE%.filteredRepeats.bed}

  if [ -f ${fungidir}/chromsizes/${SAMPLE}.txt ]; then
  	if [ ! -f ${fungidir}/te_bed_gt100_frac/${SAMPLE}.txt ]; then
	  	bedtools genomecov -max 1 -g ${fungidir}/chromsizes/${SAMPLE}.txt \
	    	-i ${i} \
	    	> ${fungidir}/te_bed_gt100_frac/${SAMPLE}.txt
	fi
  fi
done
