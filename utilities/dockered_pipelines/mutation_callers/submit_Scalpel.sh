#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long out-dir:,out-vcf:,tumor-bam:,normal-bam:,human-reference:,selector:,action: -n 'submit_Scalpel.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

#echo "$OPTS"
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

scalpel_script=${outdir}/logs/scalpel_${timestamp}.cmd

echo "#!/bin/bash" > $scalpel_script
echo "" >> $scalpel_script

echo "#$ -o ${logdir}" >> $scalpel_script
echo "#$ -e ${logdir}" >> $scalpel_script
echo "#$ -S /bin/bash" >> $scalpel_script
echo '#$ -l h_vmem=14G' >> $scalpel_script
echo 'set -e' >> $scalpel_script
echo "" >> $scalpel_script

echo 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $scalpel_script
echo "" >> $scalpel_script

echo "docker run --rm -v /:/mnt -u $UID --memory 14g -i lethalfang/scalpel:0.5.3 \\" >> $scalpel_script
echo "/opt/scalpel-0.5.3/scalpel-discovery --somatic \\" >> $scalpel_script
echo "--ref /mnt/${HUMAN_REFERENCE} \\" >> $scalpel_script
echo "--bed /mnt/${SELECTOR} \\" >> $scalpel_script
echo "--normal /mnt/${normal_bam} \\" >> $scalpel_script
echo "--tumor /mnt/${tumor_bam} \\" >> $scalpel_script
echo "--window 600 \\" >> $scalpel_script
echo "--dir /mnt/${outdir}/scalpel" >> $scalpel_script
echo "" >> $scalpel_script

echo "docker run --rm -v /:/mnt -u $UID --memory 14g -i lethalfang/scalpel:0.5.3 bash -c \\" >> $scalpel_script
echo "\"/opt/scalpel-0.5.3/scalpel-export --somatic \\" >> $scalpel_script
echo "--db /mnt/${outdir}/scalpel/main/somatic.db.dir \\" >> $scalpel_script
echo "--ref /mnt/${HUMAN_REFERENCE} \\" >> $scalpel_script
echo "--bed /mnt/${SELECTOR} \\" >> $scalpel_script
echo "> /mnt/${outdir}/scalpel/scalpel.vcf\"" >> $scalpel_script
echo "" >> $scalpel_script

echo "cat ${outdir}/scalpel/scalpel.vcf |\\" >> $scalpel_script
echo "docker run --rm -v /:/mnt -u $UID -i lethalfang/scalpel:0.5.3 /opt/vcfsorter.pl /mnt/${HUMAN_REFERENCE%\.fa*}.dict - \\" >> $scalpel_script
echo "> ${outdir}/${outvcf}" >> $scalpel_script

echo "" >> $scalpel_script
echo 'echo -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $scalpel_script

${action} $scalpel_script
