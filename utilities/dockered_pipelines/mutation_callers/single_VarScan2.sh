set -e

OPTS=`getopt -o o: --long out-dir:,out-vcf:,in-bam:,human-reference:,selector:,action:,VAF:,MEM: -n 'submit_VarScan2.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

#echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )
VAF=0.10
action=echo
MEM=7
minMQ=25
minBQ=20

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

    --minMQ )
        case "$2" in
            "") shift 2 ;;
            *) minMQ=$2 ; shift 2 ;;
        esac ;;

    --minBQ )
        case "$2" in
            "") shift 2 ;;
            *) minBQ=$2 ; shift 2 ;;
        esac ;;

    --extra-pileup-arguments )
        case "$2" in
            "") shift 2 ;;
            *) extra_pileup_arguments=$2 ; shift 2 ;;
        esac ;;

    --extra-arguments )
        case "$2" in
            "") shift 2 ;;
            *) extra_arguments=$2 ; shift 2 ;;
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

    -- ) shift; break ;;
    * ) break ;;
    esac

done


logdir=${outdir}/logs
mkdir -p ${logdir}

out_script=${outdir}/logs/varscan2_${timestamp}.cmd

selector_text=''
if [[ -r $SELECTOR ]]; then
    selector_text="-l /mnt/${SELECTOR}"
fi

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

echo "docker run --rm -u $UID -v /:/mnt --memory ${MEM}G lethalfang/samtools:0.1.19 bash -c \\" >> $out_script
echo "\"samtools mpileup \\" >> $out_script
echo "-B -q ${minMQ} -Q ${minBQ} ${extra_pileup_arguments} $selector_text -f \\" >> $out_script
echo "/mnt/${HUMAN_REFERENCE} \\" >> $out_script
echo "/mnt/${tumor_bam} \\" >> $out_script
echo "> /mnt/${outdir}/tumor.pileup\"" >> $out_script

echo "" >> $out_script

echo "docker run --rm -u $UID -v /:/mnt --memory ${MEM}G djordjeklisic/sbg-varscan2:v1 bash -c \\" >> $out_script
echo "\"java -Xmx${MEM}g -jar VarScan2.3.7.jar mpileup2cns \\" >> $out_script
echo "/mnt/${outdir}/tumor.pileup \\" >> $out_script
echo "--variants ${extra_arguments} --min-var-freq $VAF --output-vcf 1 \\" >> $out_script
echo "> /mnt/${outdir}/${outvcf}\"" >> $out_script

echo "" >> $out_script

echo "rm ${outdir}/tumor.pileup" >> $out_script
echo "" >> $out_script

echo 'echo -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $out_script

${action} $out_script
