usage()
{
cat << EOF
usage: $0 options

Map reads to metagenomic assembly.
OPTIONS:
      -i Directory containing bin files [REQUIRED]
      -s Sample name [REQUIRED]
      -o Output directory [REQUIRED]
      -b Binning tool [REQUIRED]

EOF
}

#variables
input=
sample=
output=
binner=

while getopts "i:s:o:b:h:" OPTION

do

  case ${OPTION} in
    i)
      input=${OPTARG}
      ;;
    s)
      sample=${OPTARG}
      ;;
    o)
      output=${OPTARG}
      ;;	
    b)
      binner=${OPTARG}
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

if [ ! -f "${output}"/binsummary.csv ]; then
	echo -e "sample,binner,bincount" > "${output}"/binsummary.csv
fi

count=`ls "${input}" | wc -l`
echo "${count}"

echo -e ""${sample}","${binner}","${count}"" >> "${output}"/binsummary.csv
