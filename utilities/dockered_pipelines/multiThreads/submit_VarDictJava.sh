#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long out-dir:,out-vcf:,tumor-bam:,normal-bam:,human-reference:,selector:,action:,VAF -n 'submit_VarDictJava.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )
VAF=0.05
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

VERSION='2.3.1'

logdir=${outdir}/logs
mkdir -p ${logdir}

vardict_script=${outdir}/logs/vardict_${timestamp}.cmd

docker run --rm -v /:/mnt -u $UID -i lethalfang/somaticseq:${VERSION} \
/opt/somaticseq/utilities/split_mergedBed.py \
-infile /mnt/${SELECTOR} -outfile /mnt/${outdir}/split_regions.bed

echo "#!/bin/bash" > $vardict_script
echo "" >> $vardict_script

echo "#$ -o ${logdir}" >> $vardict_script
echo "#$ -e ${logdir}" >> $vardict_script
echo "#$ -S /bin/bash" >> $vardict_script
echo '#$ -l h_vmem=4G' >> $vardict_script
echo 'set -e' >> $vardict_script
echo "" >> $vardict_script

echo 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $vardict_script
echo "" >> $vardict_script

echo "docker run --rm -v /:/mnt -u $UID -i lethalfang/vardictjava:1.5.1 \\" >> $vardict_script
echo "/opt/VarDict-1.5.1/bin/VarDict \\" >> $vardict_script
echo "-G /mnt/${HUMAN_REFERENCE} \\" >> $vardict_script
echo "-f $VAF -h \\" >> $vardict_script
echo "-b '/mnt/${tumor_bam}|/mnt/${normal_bam}' \\" >> $vardict_script
echo "-Q 1 -c 1 -S 2 -E 3 -g 4 /mnt/${outdir}/split_regions.bed \\" >> $vardict_script
echo "> ${outdir}/${timestamp}.var \\" >> $vardict_script
echo "" >> $vardict_script

echo "docker run --rm -v /:/mnt -u $UID -i lethalfang/vardictjava:1.5.1 \\" >> $vardict_script
echo "bash -c \"cat /mnt/${outdir}/${timestamp}.var | awk 'NR!=1' | /opt/VarDict/testsomatic.R | /opt/VarDict/var2vcf_paired.pl -N 'TUMOR|NORMAL' -f $VAF\" \\" >> $vardict_script
echo "> ${outdir}/${outvcf}" >> $vardict_script

echo "" >> $vardict_script
echo 'echo -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $vardict_script

${action} $vardict_script
