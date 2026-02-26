usage()
{
cat << EOF
usage: $0 options

Bin metagenomic assemblies using LRBinner

OPTIONS:
      -r Reads file [REQUIRED]
      -a Assembly file [REQUIRED]
      -o Output directory for binning[REQUIRED]

EOF
}

#variables
reads=
assembly=
binningdir=

while getopts "r:a:o:h:" OPTION

do

  case ${OPTION} in
    r)
      reads=${OPTARG}
      ;;
    a)
      assembly=${OPTARG}
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

if [ "${reads: -3}" == ".gz" ]; then
	echo -e "Unzipping reads..."
	gunzip ${reads}
	reads=${reads%.gz}
fi

mkdir -p ${binningdir}
echo -e "Performing binning using LRBinner..."
./lrbinner.py contigs --reads-path "${reads}" --output "${binningdir}" --contigs "${assembly}" --separate

if [ "${assembly: -3}" != ".gz" ]; then
	echo -e "Zipping assembly..."
	gzip ${assembly}
fi

#if [ "${reads: -3}" != ".gz" ]; then
#	echo -e "Zipping reads..."
#	gzip ${reads}
#fi

echo -e "Binning complete."
