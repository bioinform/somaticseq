#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long out-dir:,out-vcf:,tumor-bam:,normal-bam:,human-reference:,selector:,action: -n 'submit_Strelka.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )
action=echo

while true; do
    case "$1" in
    --out-vcf )
        case "$2" in
            "") shift 2 ;;
            *)  outvcf=$2 ; shift 2 ;;
        esac ;;

    --out-dir )
        case "$2" in
            "") shift 2 ;;
            *)  outdir=$2 ; shift 2 ;;
        esac ;;

    --tumor-bam )
        case "$2" in
            "") shift 2 ;;
            *)  tumor_bam=$2 ; shift 2 ;;
        esac ;;

    --normal-bam )
        case "$2" in
            "") shift 2 ;;
            *)  normal_bam=$2 ; shift 2 ;;
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
        
    --action )
        case "$2" in
            "") shift 2 ;;
            *) action=$2 ; shift 2 ;;
        esac ;;
        

    -- ) shift; break ;;
    * ) break ;;
    esac

done

vcf_prefix=${outvcf%\.vcf}

logdir=${outdir}/logs
mkdir -p ${logdir}

if [[ ${SELECTOR} ]]; then
    region_txt=`cat ${SELECTOR} | awk -F "\t" '{print "--region=" $1 ":" $2+1 "-" $3}' | tr '\n' '\ '`
fi

strelka_script=${outdir}/logs/strelka_${timestamp}.cmd

echo "#!/bin/bash" > $strelka_script
echo "" >> $strelka_script

echo "#$ -o ${logdir}" >> $strelka_script
echo "#$ -e ${logdir}" >> $strelka_script
echo "#$ -S /bin/bash" >> $strelka_script
echo '#$ -l h_vmem=6G' >> $strelka_script
echo 'set -e' >> $strelka_script
echo "" >> $strelka_script

echo 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $strelka_script
echo "" >> $strelka_script

echo "docker run -v /:/mnt -i lethalfang/strelka:2.7.1 \\" >> $strelka_script
echo "/opt/strelka2/bin/configureStrelkaSomaticWorkflow.py \\" >> $strelka_script
echo "--tumorBam=/mnt/${tumor_bam} \\" >> $strelka_script
echo "--normalBam=/mnt/${normal_bam} \\" >> $strelka_script
echo "--referenceFasta=/mnt/${HUMAN_REFERENCE}  \\" >> $strelka_script
echo "--callMemMb=4096 \\" >> $strelka_script
echo "$region_txt \\" >> $strelka_script
echo "--runDir=/mnt/${outdir}/${outvcf%\.vcf}" >> $strelka_script
echo "" >> $strelka_script

echo "docker run -v /:/mnt -i lethalfang/strelka:2.7.1 \\" >> $strelka_script
echo "/mnt/${outdir}/${outvcf%\.vcf}/runWorkflow.py -m local -j 1" >> $strelka_script

echo "" >> $strelka_script
echo 'echo -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $strelka_script

${action} $strelka_script
