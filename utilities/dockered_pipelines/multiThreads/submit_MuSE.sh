#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long out-dir:,out-vcf:,tumor-bam:,normal-bam:,human-reference:,selector:,action:,dbsnp: -n 'submit_MuSE.sh'  -- "$@"`

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

    --dbsnp )
        case "$2" in
            "") shift 2 ;;
            *) dbsnp=$2 ; shift 2 ;;
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

muse_script=${outdir}/logs/muse_${timestamp}.cmd

echo "#!/bin/bash" > $muse_script
echo "" >> $muse_script

echo "#$ -o ${logdir}" >> $muse_script
echo "#$ -e ${logdir}" >> $muse_script
echo "#$ -S /bin/bash" >> $muse_script
echo '#$ -l h_vmem=4G' >> $muse_script
echo 'set -e' >> $muse_script
echo "" >> $muse_script

echo 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $muse_script
echo "" >> $muse_script

echo "cat ${SELECTOR} | awk -F \"\t\" '{print \$1 \"\t\" \$2 \"\t\" \$3}' > ${outdir}/bed_3columns.bed" >> $muse_script
echo "" >> $muse_script

echo "docker run -v /:/mnt -i marghoob/muse:1.0rc_c \\" >> $muse_script
echo "MuSEv1.0rc_submission_c039ffa call \\" >> $muse_script
echo "-O /mnt/${outdir}/MuSE \\" >> $muse_script
echo "-l /mnt/${outdir}/bed_3columns.bed \\" >> $muse_script
echo "-f /mnt/${HUMAN_REFERENCE} \\" >> $muse_script
echo "/mnt/${tumor_bam} \\" >> $muse_script
echo "/mnt/${normal_bam}" >> $muse_script
echo "" >> $muse_script

echo "docker run -v /:/mnt -i marghoob/muse:1.0rc_c \\" >> $muse_script
echo "MuSEv1.0rc_submission_c039ffa sump \\" >> $muse_script
echo "-I /mnt/${outdir}/MuSE.MuSE.txt \\" >> $muse_script
echo "-E -O /mnt/${outdir}/${outvcf} \\" >> $muse_script
echo "-D /mnt/${dbsnp}.gz" >> $muse_script

echo "" >> $muse_script
echo 'echo -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $muse_script

${action} $muse_script
