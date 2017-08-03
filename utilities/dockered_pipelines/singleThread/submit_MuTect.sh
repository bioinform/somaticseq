#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long out-dir:,out-vcf:,tumor-bam:,normal-bam:,human-reference:,selector:,dbsnp:,cosmic:,action: -n 'submit_MuTect.sh'  -- "$@"`

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

    --dbsnp )
        case "$2" in
            "") shift 2 ;;
            *)  dbsnp=$2 ; shift 2 ;;
        esac ;;

    --cosmic )
        case "$2" in
            "") shift 2 ;;
            *)  cosmic=$2 ; shift 2 ;;
        esac ;;

    --action )
        case "$2" in
            "") shift 2 ;;
            *)  action=$2 ; shift 2 ;;
        esac ;;


    -- ) shift; break ;;
    * ) break ;;
    esac

done

logdir=${outdir}/logs
mkdir -p ${logdir}

mutect_script=${outdir}/logs/mutect_${timestamp}.cmd

selector_text=''
if [[ -r $SELECTOR ]]; then
    selector_text="-L /mnt/${SELECTOR}"
fi

dbsnp_text=''
if [[ -r $dbsnp ]]; then
    dbsnp_text="--dbsnp /mnt/${dbsnp}"
fi

cosmic_text=''
if [[ -r $cosmic ]]; then
    cosmic_text="--cosmic /mnt/${cosmic}"
fi


echo "#!/bin/bash" > $mutect_script
echo "" >> $mutect_script

echo "#$ -o ${logdir}" >> $mutect_script
echo "#$ -e ${logdir}" >> $mutect_script
echo "#$ -S /bin/bash" >> $mutect_script
echo '#$ -l h_vmem=10G' >> $mutect_script
echo 'set -e' >> $mutect_script
echo "" >> $mutect_script

echo 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $mutect_script
echo "" >> $mutect_script

echo "docker run -v /:/mnt -i rpbuv011:5000/fangl10/mutect:2014.4-8-gf14306c \\" >> $mutect_script
echo "java -Xmx7g -jar /opt/CancerAnalysisPackage-2014.4-8-gf14306c/SomaticAnalysisTK.jar -T MuTect \\" >> $mutect_script
echo "--reference_sequence /mnt/${HUMAN_REFERENCE} \\" >> $mutect_script
echo "$selector_text \\" >> $mutect_script
echo "--input_file:normal /mnt/${normal_bam} \\" >> $mutect_script
echo "--input_file:tumor /mnt/${tumor_bam} \\" >> $mutect_script
echo "$dbsnp_text \\" >> $mutect_script
echo "$cosmic_text \\" >> $mutect_script
echo "--vcf /mnt/${outdir}/${outvcf} \\" >> $mutect_script
echo "--out /dev/null" >> $mutect_script
echo "" >> $mutect_script

echo "docker run -v /:/mnt -i rpbuv011:5000/fangl10/mutect:2014.4-8-gf14306c \\" >> $mutect_script
echo "java -Xmx7g -jar /opt/CancerAnalysisPackage-2014.4-8-gf14306c/SomaticAnalysisTK.jar -T SomaticIndelDetector \\" >> $mutect_script
echo "--reference_sequence /mnt/${HUMAN_REFERENCE} \\" >> $mutect_script
echo "-L /mnt/${SELECTOR} \\" >> $mutect_script
echo "--input_file:normal /mnt/${normal_bam} \\" >> $mutect_script
echo "--input_file:tumor /mnt/${tumor_bam} \\" >> $mutect_script
echo "--out /mnt/${outdir}/Indel.${outvcf}" >> $mutect_script

echo "" >> $mutect_script
echo 'echo -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $mutect_script

${action} $mutect_script
