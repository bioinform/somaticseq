#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long out-dir:,out-vcf:,tumor-bam:,normal-bam:,human-reference:,mapping-quality:,base-quality:,prior:,split:,action: -n 'submit_SomaticSniper.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

#echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )

MQ=25
BQ=15
prior=0.0001
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

    --mapping-quality )
        case "$2" in
            "") shift 2 ;;
            *)  MQ=$2 ; shift 2 ;;
        esac ;;

    --base-quality )
        case "$2" in
            "") shift 2 ;;
            *)  BQ=$2 ; shift 2 ;;
        esac ;;

    --prior-probability )
        case "$2" in
            "") shift 2 ;;
            *)  prior=$2 ; shift 2 ;;
        esac ;;

    --action )
        case "$2" in
            "") shift 2 ;;
            *)  action=$2 ; shift 2 ;;
        esac ;;

    --split )
        case "$2" in
            "") shift 2 ;;
            *)  split=$2 ; shift 2 ;;
        esac ;;

    -- ) shift; break ;;
    * ) break ;;
    esac

done

logdir=${outdir}/logs
mkdir -p ${logdir}

sniper_script=${outdir}/logs/somaticsniper_${timestamp}.cmd

echo "#!/bin/bash" > $sniper_script
echo "" >> $sniper_script

echo "#$ -o ${logdir}" >> $sniper_script
echo "#$ -e ${logdir}" >> $sniper_script
echo "#$ -S /bin/bash" >> $sniper_script
echo '#$ -l h_vmem=6G' >> $sniper_script
echo 'set -e' >> $sniper_script
echo "" >> $sniper_script

echo 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $sniper_script
echo "" >> $sniper_script

echo "docker run --rm -v /:/mnt -u $UID -i lethalfang/somaticsniper:1.0.5.0 \\" >> $sniper_script
echo "/opt/somatic-sniper/build/bin/bam-somaticsniper \\" >> $sniper_script
echo "-q ${MQ} -Q ${BQ} -s ${prior} -F vcf \\" >> $sniper_script
echo "-f /mnt/${HUMAN_REFERENCE} \\" >> $sniper_script
echo "/mnt/${tumor_bam} \\" >> $sniper_script
echo "/mnt/${normal_bam} \\" >> $sniper_script
echo "/mnt/${outdir}/${outvcf}" >> $sniper_script
echo "" >> $sniper_script

if [[ $split ]]
then
    echo "i=1" >> $sniper_script
    echo "while [[ \$i -le $split ]]" >> $sniper_script
    echo "do" >> $sniper_script
    echo "    docker run --rm -v /:/mnt -u $UID -i lethalfang/somaticseq:base-1.0 bedtools intersect -a /mnt/${outdir}/${outvcf} -b /mnt/${outdir}/\${i}/\${i}.bed -header | uniq >  ${outdir}/\${i}/${outvcf}" >> $sniper_script
    echo "    i=\$(( \$i + 1 ))" >> $sniper_script
    echo "done" >> $sniper_script
fi

echo "" >> $sniper_script
echo 'echo -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $sniper_script

${action} $sniper_script
