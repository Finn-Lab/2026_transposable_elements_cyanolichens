usage()
{
cat << EOF
usage: $0 options

Map reads to metagenomic assembly.
OPTIONS:
      -r Forward reads [REQUIRED]
      -a Metagenomic assembly [REQUIRED]
      -o Output directory path [REQUIRED]

EOF
}

#variables
reads=
assembly=
bamoutput=

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
      bamoutput=${OPTARG}
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

#if [ "${assembly: -3}" == ".gz" ]; then
#	echo -e "Unzipping assembly..."
#	gunzip ${assembly}
#	assembly=${assembly%.gz}
#fi

if [ "${reads: -3}" == ".gz" ]; then
	echo -e "Unzipping reads..."
	gunzip ${reads}
	reads=${reads%.gz}
fi

echo -e "Mapping reads to assembly with minimap2..."
minimap2 -ax map-hifi ${assembly} ${reads} | samtools sort -o "${bamoutput}"

samtools index "${bamoutput}"

if [ "${reads: -3}" != ".gz" ]; then
	echo -e "Zipping reads..."
	gzip ${reads}
fi

#if [ "${assembly: -3}" != ".gz" ]; then
#	echo -e "Zipping assembly..."
#	gzip ${assembly}
#fi

echo -e "Mapping complete"
