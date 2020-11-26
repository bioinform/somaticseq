#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long out-dir:,out-vcf:,bam:,human-reference:,selector:,dbsnp:,extra-arguments:,action:,MEM:,threads:,out-script:,standalone -n 'gatk_haplotypecaller.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

#echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )
action=echo
MEM=8
threads=12

while true; do
    case "$1" in
    
    -o | --out-dir )
        case "$2" in
            "") shift 2 ;;
            *)  outdir=$2 ; shift 2 ;;
        esac ;;
    
    --out-vcf )
        case "$2" in
            "") shift 2 ;;
            *)  outVcfName=$2 ; shift 2 ;;
        esac ;;

    --bam )
        case "$2" in
            "") shift 2 ;;
            *)  bamFile=$2 ; shift 2 ;;
        esac ;;

    --human-reference )
        case "$2" in
            "") shift 2 ;;
            *)  HUMAN_REFERENCE=$2 ; shift 2 ;;
        esac ;;

    --selector )
        case "$2" in
            "") shift 2 ;;
            *) SELECTOR=$2 ; shift 2 ;;
        esac ;;

    --dbsnp )
        case "$2" in
            "") shift 2 ;;
            *)  dbsnp=$2 ; shift 2 ;;
        esac ;;

    --MEM )
        case "$2" in
            "") shift 2 ;;
            *)  MEM=$2 ; shift 2 ;;
        esac ;;

    --threads )
        case "$2" in
            "") shift 2 ;;
            *)  threads=$2 ; shift 2 ;;
        esac ;;

    --extra-arguments )
        case "$2" in
            "") shift 2 ;;
            *)  extra_arguments=$2 ; shift 2 ;;
        esac ;;

    --out-script )
        case "$2" in
            "") shift 2 ;;
            *)  out_script_name=$2 ; shift 2 ;;
        esac ;;

    --action )
        case "$2" in
            "") shift 2 ;;
            *)  action=$2 ; shift 2 ;;
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
    out_script="${logdir}/HaplotypeCaller.${timestamp}.cmd"    
fi


if [[ $standalone ]]
then
    echo "#!/bin/bash" > $out_script
    echo "" >> $out_script
    echo "#$ -o ${logdir}" >> $out_script
    echo "#$ -e ${logdir}" >> $out_script
    echo "#$ -S /bin/bash" >> $out_script
    echo '#$ -l h_vmem=12G' >> $out_script
    echo "#$ -pe smp ${threads}" >> $out_script
    echo 'set -e' >> $out_script
fi

echo "" >> $out_script
echo 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $out_script
echo "" >> $out_script


if [[ ${SELECTOR} ]]
then
    selector_text="-L /mnt/${SELECTOR}"
fi

dbsnp_text=''
if [[ ${dbsnp} ]]; then
    dbsnp_text="--dbsnp /mnt/${dbsnp}"
fi


echo "docker run --rm -v /:/mnt -u $UID broadinstitute/gatk:4.0.5.2 \\" >> $out_script
echo "java -Xmx${MEM}g -jar /gatk/gatk.jar \\" >> $out_script
echo "HaplotypeCaller \\" >> $out_script
echo "--reference /mnt/${HUMAN_REFERENCE} \\" >> $out_script
echo "--input /mnt/${bamFile} \\" >> $out_script
echo "--native-pair-hmm-threads ${threads} \\" >> $out_script
echo "$selector_text \\" >> $out_script
echo "$dbsnp_text \\" >> $out_script
echo "${extra_arguments} \\" >> $out_script
echo "--output /mnt/${outdir}/${outVcfName}" >> $out_script

echo "" >> $out_script

echo 'echo -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $out_script
