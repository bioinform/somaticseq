set -e

OPTS=`getopt -o o: --long out-dir:,out-vcf:,tumor-bam:,normal-bam:,human-reference:,selector:,action:,VAF: -n 'submit_VarScan2.sh'  -- "$@"`

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
echo '#$ -l h_vmem=9G' >> $varscan2_script
echo 'set -e' >> $varscan2_script
echo "" >> $varscan2_script

echo 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $varscan2_script
echo "" >> $varscan2_script

echo "docker run --rm -u $UID -v /:/mnt --memory 8g lethalfang/samtools:0.1.19 \\" >> $varscan2_script
echo "samtools mpileup \\" >> $varscan2_script
echo "-B -q 25 -Q 20 $selector_text -f \\" >> $varscan2_script
echo "/mnt/${HUMAN_REFERENCE} \\" >> $varscan2_script
echo "/mnt/${normal_bam} \\" >> $varscan2_script
echo "> ${outdir}/normal.pileup" >> $varscan2_script

echo "" >> $varscan2_script

echo "docker run --rm -u $UID -v /:/mnt --memory 8g lethalfang/samtools:0.1.19 \\" >> $varscan2_script
echo "samtools mpileup \\" >> $varscan2_script
echo "-B -q 25 -Q 20 $selector_text -f \\" >> $varscan2_script
echo "/mnt/${HUMAN_REFERENCE} \\" >> $varscan2_script
echo "/mnt/${tumor_bam} \\" >> $varscan2_script
echo "> ${outdir}/tumor.pileup" >> $varscan2_script

echo "" >> $varscan2_script

echo "docker run --rm -u $UID -v /:/mnt --memory 8g djordjeklisic/sbg-varscan2:v1 \\" >> $varscan2_script
echo "java -Xmx8g -jar VarScan2.3.7.jar somatic \\" >> $varscan2_script
echo "/mnt/${outdir}/normal.pileup \\" >> $varscan2_script
echo "/mnt/${outdir}/tumor.pileup \\" >> $varscan2_script
echo "/mnt/${outdir}/${outvcf%.vcf} --output-vcf 1 --min-var-freq $VAF" >> $varscan2_script

echo "" >> $varscan2_script

echo "docker run --rm -u $UID -v /:/mnt --memory 8g djordjeklisic/sbg-varscan2:v1 \\" >> $varscan2_script
echo "java -Xmx8g -jar VarScan2.3.7.jar processSomatic \\" >> $varscan2_script
echo "/mnt/${outdir}/${outvcf%.vcf}.snp.vcf" >> $varscan2_script

echo "" >> $varscan2_script

echo "docker run --rm -u $UID -v /:/mnt --memory 8g djordjeklisic/sbg-varscan2:v1 \\" >> $varscan2_script
echo "java -Xmx8g -jar VarScan2.3.7.jar somaticFilter \\" >> $varscan2_script
echo "/mnt/${outdir}/${outvcf%.vcf}.snp.Somatic.hc.vcf \\" >> $varscan2_script
echo "-indel-file /mnt/${outdir}/${outvcf%.vcf}.indel.vcf \\" >> $varscan2_script
echo "-output-file /mnt/${outdir}/${outvcf%.vcf}.snp.Somatic.hc.filter.vcf" >> $varscan2_script

echo "" >> $varscan2_script

echo "rm ${outdir}/normal.pileup" >> $varscan2_script
echo "rm ${outdir}/tumor.pileup" >> $varscan2_script
echo "" >> $varscan2_script

echo 'echo -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $varscan2_script

${action} $varscan2_script
