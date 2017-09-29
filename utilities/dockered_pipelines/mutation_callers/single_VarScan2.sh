set -e

OPTS=`getopt -o o: --long out-dir:,out-vcf:,in-bam:,human-reference:,selector:,action:,VAF: -n 'submit_VarScan2.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

#echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )
VAF=0.10
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

    --in-bam )
        case "$2" in
            "") shift 2 ;;
            *)  tumor_bam=$2 ; shift 2 ;;
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

    --VAF )
        case "$2" in
            "") shift 2 ;;
            *) VAF=$2 ; shift 2 ;;
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


logdir=${outdir}/logs
mkdir -p ${logdir}

varscan2_script=${outdir}/logs/varscan2_${timestamp}.cmd

selector_text=''
if [[ -r $SELECTOR ]]; then
    selector_text="-l /mnt/${SELECTOR}"
fi

echo "#!/bin/bash" > $varscan2_script
echo "" >> $varscan2_script

echo "#$ -o ${logdir}" >> $varscan2_script
echo "#$ -e ${logdir}" >> $varscan2_script
echo "#$ -S /bin/bash" >> $varscan2_script
echo '#$ -l h_vmem=6G' >> $varscan2_script
echo 'set -e' >> $varscan2_script
echo "" >> $varscan2_script

echo 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $varscan2_script
echo "" >> $varscan2_script

echo "docker run --rm -u $UID -v /:/mnt --memory 8g lethalfang/samtools:0.1.19 bash -c \\" >> $varscan2_script
echo "\"samtools mpileup \\" >> $varscan2_script
echo "-B -q 25 -Q 20 $selector_text -f \\" >> $varscan2_script
echo "/mnt/${HUMAN_REFERENCE} \\" >> $varscan2_script
echo "/mnt/${tumor_bam} \\" >> $varscan2_script
echo "> /mnt/${outdir}/tumor.pileup\"" >> $varscan2_script

echo "" >> $varscan2_script

echo "docker run --rm -u $UID -v /:/mnt --memory 8g djordjeklisic/sbg-varscan2:v1 bash -c \\" >> $varscan2_script
echo "\"java -Xmx6g -jar VarScan2.3.7.jar mpileup2cns \\" >> $varscan2_script
echo "/mnt/${outdir}/tumor.pileup \\" >> $varscan2_script
echo "--variants --min-var-freq $VAF --output-vcf 1" >> $varscan2_script
echo "> /mnt/${outdir}/${outvcf}\"" >> $varscan2_script

echo "" >> $varscan2_script

echo "rm ${outdir}/tumor.pileup" >> $varscan2_script
echo "" >> $varscan2_script

echo 'echo -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $varscan2_script

${action} $varscan2_script
