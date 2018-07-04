#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long out-dir:,out-vcf:,tumor-bam:,normal-bam:,human-reference:,selector:,extra-config-arguments:,extra-run-arguments:,MEM:,action:,exome -n 'submit_Strelka.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

#echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )
action=echo
MEM=6

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

    --extra-config-arguments )
        case "$2" in
            "") shift 2 ;;
            *) extra_config_arguments=$2 ; shift 2 ;;
        esac ;;

    --extra-run-arguments )
        case "$2" in
            "") shift 2 ;;
            *) extra_run_arguments=$2 ; shift 2 ;;
        esac ;;

    --MEM )
        case "$2" in
            "") shift 2 ;;
            *) MEM=$2 ; shift 2 ;;
        esac ;;
        
    --action )
        case "$2" in
            "") shift 2 ;;
            *) action=$2 ; shift 2 ;;
        esac ;;

    --exome )
        if_exome=1 ; shift ;;

    -- ) shift; break ;;
    * ) break ;;
    esac

done

vcf_prefix=${outvcf%\.vcf}

logdir=${outdir}/logs
mkdir -p ${logdir}

if [[ ${SELECTOR} &&  `cat ${SELECTOR} | wc -l` -lt 50 ]]; then
    region_txt=`cat ${SELECTOR} | awk -F "\t" '{print "--region=" $1 ":" $2+1 "-" $3}' | tr '\n' '\ '`
fi


out_script=${outdir}/logs/strelka_${timestamp}.cmd

echo "#!/bin/bash" > $out_script
echo "" >> $out_script

echo "#$ -o ${logdir}" >> $out_script
echo "#$ -e ${logdir}" >> $out_script
echo "#$ -S /bin/bash" >> $out_script
echo "#$ -l h_vmem=${MEM}G" >> $out_script
echo 'set -e' >> $out_script
echo "" >> $out_script

echo 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $out_script
echo "" >> $out_script


selector_basename=`basename ${SELECTOR}`
if [[ -r ${SELECTOR}.gz && -r ${SELECTOR}.gz.tbi ]]
then
    input_BED=${SELECTOR}.gz
    
else
    echo "docker run -v /:/mnt -u $UID --rm --memory ${MEM}G lethalfang/tabix:1.7 bash -c \"cat /mnt/${SELECTOR} | bgzip > /mnt/${outdir}/${selector_basename}.gz\"" >> $out_script
    echo "docker run -v /:/mnt -u $UID --rm --memory ${MEM}G lethalfang/tabix:1.7 tabix /mnt/${outdir}/${selector_basename}.gz" >> $out_script
    echo "" >> $out_script
    
    input_BED=${outdir}/${selector_basename}.gz
fi

if [[ $if_exome ]]
then
    exome='--exome'
fi

echo "docker run --rm -v /:/mnt -u $UID --memory ${MEM}G lethalfang/strelka:2.9.5 \\" >> $out_script
echo "/opt/strelka/bin/configureStrelkaSomaticWorkflow.py \\" >> $out_script
echo "--tumorBam=/mnt/${tumor_bam} \\" >> $out_script
echo "--normalBam=/mnt/${normal_bam} \\" >> $out_script
echo "--referenceFasta=/mnt/${HUMAN_REFERENCE}  \\" >> $out_script
echo "--callMemMb=$(( 1024 * MEM )) \\" >> $out_script
echo "$region_txt \\" >> $out_script
echo "--callRegions=/mnt/${input_BED} \\" >> $out_script
echo "${exome} ${extra_config_arguments} \\" >> $out_script
echo "--runDir=/mnt/${outdir}/${outvcf%\.vcf}" >> $out_script

echo "" >> $out_script

echo "docker run --rm -v /:/mnt -u $UID --memory ${MEM}G lethalfang/strelka:2.9.5 \\" >> $out_script
echo "/mnt/${outdir}/${outvcf%\.vcf}/runWorkflow.py -m local -j 1 ${extra_run_arguments}" >> $out_script

echo "" >> $out_script
echo 'echo -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $out_script

${action} $out_script
