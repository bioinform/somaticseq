#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long output-dir:,vcf-string:,vcf-out:,out-script:,standalone -n 'MergeTN.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

#echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )

keep_intermediates=0

while true; do
    case "$1" in
        -o | --output-dir )
            case "$2" in
                "") shift 2 ;;
                *)  outdir=$2 ; shift 2 ;;
            esac ;;
            
        --vcf-out )
            case "$2" in
                "") shift 2 ;;
                *)  outvcf=$2 ; shift 2 ;;
            esac ;;

        --vcf-string )
            case "$2" in
                "") shift 2 ;;
                *)  vcf_string=$2 ; shift 2 ;;
            esac ;;

        --out-script )
            case "$2" in
                "") shift 2 ;;
                *)  out_script_name=$2 ; shift 2 ;;
            esac ;;

        --standalone )
            standalone=1 ; shift ;;

        -- ) shift; break ;;
        * ) break ;;
    esac
done

logdir=${outdir}/logs
mkdir -p ${logdir}

if [[ ${out_script_name} ]]
then
    out_script="${out_script_name}"
else
    out_script="${logdir}/concatVcfFiles.${timestamp}.cmd"    
fi

if [[ $standalone ]]
then
    echo "#!/bin/bash" > $out_script
    echo "" >> $out_script
    echo "#$ -o ${logdir}" >> $out_script
    echo "#$ -e ${logdir}" >> $out_script
    echo "#$ -S /bin/bash" >> $out_script
    echo '#$ -l h_vmem=2G' >> $out_script
    echo 'set -e' >> $out_script
fi

echo "" >> $out_script


for file in ${vcf_string}
do
    input_file_string="/mnt/${file} ${input_file_string}"
done

# Merge the BAM files
echo "docker run -v /:/mnt -u $UID --memory 2g --rm lethalfang/vcftools:0.1.15 bash -c \\" >> $out_script
echo "\"vcf-concat \\" >> $out_script
echo "${input_file_string} \\" >> $out_script
echo "> /mnt/${outdir}/${outvcf}\"" >> $out_script
echo "" >> $out_script
