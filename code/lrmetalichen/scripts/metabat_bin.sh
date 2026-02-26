usage()
{
cat << EOF
usage: $0 options

Bin metagenomic assemblies using metabat2

OPTIONS:
      -i  Input metagenomic assembly, either contigs.fasta or scaffolds.fasta [REQUIRED]
      -m Mapped reads alignment bam [REQUIRED]
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

echo -e "Calculating metabat depth..."
jgi_summarize_bam_contig_depths --outputDepth "${binningdir}"/metabat_depth.txt "${mappedreads}"

echo -e "Performing binning using metaBAT2..."
metabat2 -i "${assembly}" -a "${binningdir}"/metabat_depth.txt -o "${binningdir}"/bin --unbinned

mkdir -p "${binningdir}"/bins

mv "${binningdir}"/*.fa "${binningdir}"/bins

if [ -f "${binningdir}"/bins/*unbinned.fa ]; then
	mv "${binningdir}"/bins/*unbinned.fa "${binningdir}"
fi

if [ -f "${binningdir}"/bins/*lowDepth.fa ]; then
	mv "${binningdir}"/bins/*lowDepth.fa "${binningdir}"
fi

if [ -f "${binningdir}"/tooShort.fa ]; then
	mv "${binningdir}"/bins/*tooShort.fa "${binningdir}"
fi

echo -e "Binning complete."
