usage()
{
cat << EOF
usage: $0 options

Bin metagenomic assemblies using CONCOCT

OPTIONS:
      -i  Input metagenomic assembly, either contigs.fasta or scaffolds.fasta [REQUIRED]
      -m  Mapped reads to assembly BAM [REQUIRED]
      -o Output directory for binning[REQUIRED]

EOF
}

#variables
assembly=
mappedreads=
binningdir=

while getopts "i:m:o:h:" OPTION

do

  case ${OPTION} in
    i)
      assembly=${OPTARG}
      ;;
    m)
      mappedreads=${OPTARG}
      ;;
    o)
      binningdir=${OPTARG}
      ;;
    h)
      usage
      exit
      ;;
    ?)
      usage
      exit
      ;;
  esac

done

if [ "${assembly: -3}" == ".gz" ]; then
	echo -e "Unzipping assembly..."
	gunzip ${assembly}
	assembly=${assembly%.gz}
fi

echo -e "Chunking up the assembly..."
cut_up_fasta.py "${assembly}" -c 10000 -o 0 --merge_last -b "${binningdir}"/contigs_10K.bed > "${binningdir}"/contigs_10K.fa

echo -e "Generating coverage table for chunks..."
concoct_coverage_table.py "${binningdir}"/contigs_10K.bed "${mappedreads}" > "${binningdir}"/coverage_table.tsv

echo -e "Running cononcoct..."
concoct --composition_file "${binningdir}"/contigs_10K.fa --coverage_file "${binningdir}"/coverage_table.tsv -b "${binningdir}"/concoct_output/ -t 32 -s 77

echo -e "Merging clusters..."
merge_cutup_clustering.py "${binningdir}"/concoct_output/clustering_gt1000.csv > "${binningdir}"/concoct_output/clustering_merged.csv

mkdir -p "${binningdir}"/concoct_output/fasta_bins
extract_fasta_bins.py "${assembly}" "${binningdir}"/concoct_output/clustering_merged.csv --output_path "${binningdir}"/concoct_output/fasta_bins

if [ "${assembly}: -3}" != ".gz" ]; then
	echo -e "Zipping assembly..."
	gzip ${assembly}
fi

echo -e "Binning complete."
