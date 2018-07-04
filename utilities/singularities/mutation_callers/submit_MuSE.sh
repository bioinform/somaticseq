#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long out-dir:,out-vcf:,tumor-bam:,normal-bam:,human-reference:,selector:,extra-arguments:,action:,MEM:,dbsnp: -n 'submit_MuSE.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

#echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )
VAF=0.05
action=echo
MEM=4

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

out_script=${outdir}/logs/muse_${timestamp}.cmd

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

echo "cat ${SELECTOR} | awk -F \"\t\" '{print \$1 \"\t\" \$2 \"\t\" \$3}' > ${outdir}/bed_3columns.bed" >> $out_script
echo "" >> $out_script

echo "singularity exec --bind /:/mnt   docker://marghoob/muse:1.0rc_c \\" >> $out_script
echo "MuSEv1.0rc_submission_c039ffa call \\" >> $out_script
echo "-O /mnt/${outdir}/MuSE \\" >> $out_script
echo "-l /mnt/${outdir}/bed_3columns.bed \\" >> $out_script
echo "-f /mnt/${HUMAN_REFERENCE} \\" >> $out_script
echo "/mnt/${tumor_bam} \\" >> $out_script
echo "/mnt/${normal_bam}" >> $out_script
echo "" >> $out_script

echo "singularity exec --bind /:/mnt   docker://marghoob/muse:1.0rc_c \\" >> $out_script
echo "MuSEv1.0rc_submission_c039ffa sump \\" >> $out_script
echo "-I /mnt/${outdir}/MuSE.MuSE.txt \\" >> $out_script
echo "${extra_arguments} -O /mnt/${outdir}/${outvcf} \\" >> $out_script
echo "-D /mnt/${dbsnp}.gz" >> $out_script

echo "" >> $out_script
echo 'echo -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $out_script

${action} $out_script
