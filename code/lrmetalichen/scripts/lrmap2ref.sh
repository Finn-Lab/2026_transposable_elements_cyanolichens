#!/bin/bash

usage()
{
cat << EOF
usage: $0 options

Map metagenomic reads against a custom-made bwa reference
and output statistics to infer genome prevalence and abundance.

'bsub' with '-M 250000'

OPTIONS:
   -t      Number of threads [REQUIRED]
   -i      Input first or only read (.fastq or .fastq.gz format) [REQUIRED]
   -r      Reference database (.fa or .fasta) [REQUIRED]
   -o      Output files prefix (include path) [REQUIRED]
   -c      Whether entries in ref database are "contigs" or "complete" genomes [REQUIRED]
EOF
}

# variables
threads=
ref=
reads=
outprefix=
mode=

while getopts “t:m:i:r:o:c:” OPTION
do
     case ${OPTION} in
         t)
             threads=${OPTARG}
             ;;
         i)
             reads=${OPTARG}
             ;;
         r)
             ref=${OPTARG}
             ;;
         o)
             outprefix=${OPTARG}
             ;;
         c)
             mode=${OPTARG}
             ;;
         ?)
             usage
             exit
             ;;

     esac
done

# check arguments
if [[ -z ${threads} ]] || [[ -z ${reads} ]] || [[ -z ${ref} ]] || [[ -z ${outprefix} ]] || [[ -z ${mode} ]]
then
     echo "ERROR : Please supply correct arguments"
     usage
     exit 1
fi

timestamp() {
  date +"%H:%M:%S"
}

readname=$(basename ${outprefix})
refix=$(basename ${ref%%.fna*})
if [ ${threads} -eq 1 ]
then
    threads_sam=1
else
    threads_sam=$((${threads}-1))
fi


# initial mapping and sorting
echo "$(timestamp) Running minimap2 ..."
minimap2 -ax map-hifi ${ref} ${reads} | samtools view  -F 256 -uS | samtools sort  -o ${outprefix}_${refix}_raw.bam
samtools index  ${outprefix}_${refix}_raw.bam
lrmetalichen/scripts/bam_ani-filter.py ${outprefix}_${refix}_raw.bam 90 60 | samtools view -u -o ${outprefix}_${refix}_sorted.bam


#remove tmp files (1)
echo "$(timestamp) [ map2ref pipeline ] Cleaning tmp files: *_raw.bam  ..."
rm -rf ${outprefix}_${refix}_raw.ba*


# parse total count output
echo "$(timestamp)  Parsing total counts results ..."
samtools index  ${outprefix}_${refix}_sorted.bam
samtools idxstats ${outprefix}_${refix}_sorted.bam > ${outprefix}_${refix}_depth.tab
samtools depth ${outprefix}_${refix}_sorted.bam > ${outprefix}_${refix}_depth-pos.tab
lrmetalichen/scripts/parse_bwa-depth.py ${outprefix}_${refix}_depth.tab ${outprefix}_${refix}_depth-pos.tab ${mode} _${refix} > ${outprefix}_${refix}_total.tab


#remove tmp files (2)
echo "$(timestamp) [ map2ref pipeline ] Cleaning tmp files: total counts files  ..."
rm -rf ${outprefix}_${refix}_depth*

# extract unique counts
echo "$(timestamp) [ map2ref pipeline ] Extracting unique counts ..."

samtools view  -q 1 -u ${outprefix}_${refix}_sorted.bam -o ${outprefix}_${refix}_unique_sorted.bam
samtools  index  ${outprefix}_${refix}_unique_sorted.bam


echo "$(timestamp) [ map2ref pipeline ] Cleaning tmp files: *_sorted.bam  ..."
rm -rf ${outprefix}_${refix}_sorted.ba*


# parse unique count output
echo "$(timestamp) [ map2ref pipeline ] Parsing results ..."
samtools idxstats ${outprefix}_${refix}_unique_sorted.bam > ${outprefix}_${refix}_unique_depth.tab
samtools depth ${outprefix}_${refix}_unique_sorted.bam > ${outprefix}_${refix}_unique_depth-pos.tab
lrmetalichen/scripts/parse_bwa-depth.py ${outprefix}_${refix}_unique_depth.tab ${outprefix}_${refix}_unique_depth-pos.tab ${mode} _${refix} > ${outprefix}_${refix}_unique.tab

#clean up tmp files (3)
echo "$(timestamp) [ map2ref pipeline ] Cleaning tmp files: unique_sorted.bams  ..."
rm -rf ${outprefix}_${refix}_unique_sorted.ba* 
rm -rf ${outprefix}_${refix}_unique_depth*

#final step
echo "$(timestamp) [ map2ref pipeline ] Mapping of ${outprefix} in ${mode} FINISHED!!!"
echo "Final results in ${outprefix}_${refix}_total.tab, ${outprefix}_${refix}_unique.tab."

