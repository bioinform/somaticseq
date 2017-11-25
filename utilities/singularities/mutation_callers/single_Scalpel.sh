#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long out-dir:,out-vcf:,in-bam:,human-reference:,selector:,extra-discovery-arguments:,extra-export-arguments:,action:,MEM:,two-pass -n 'submit_Scalpel.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

#echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )
VAF=0.05
action=echo
MEM=16

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

    --extra-discovery-arguments )
        case "$2" in
            "") shift 2 ;;
            *)  extra_discovery_arguments=$2 ; shift 2 ;;
        esac ;;

    --extra-export-arguments )
        case "$2" in
            "") shift 2 ;;
            *)  extra_export_arguments=$2 ; shift 2 ;;
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

    --two-pass )
        two_pass=1 ; shift ;;
        
    -- ) shift; break ;;
    * ) break ;;
    esac

done

vcf_prefix=${outvcf%\.vcf}

logdir=${outdir}/logs
mkdir -p ${logdir}

out_script=${outdir}/logs/scalpel_${timestamp}.cmd

echo "#!/bin/bash" > $out_script
echo "" >> $out_script

echo "#$ -o ${logdir}" >> $out_script
echo "#$ -e ${logdir}" >> $out_script
echo "#$ -S /bin/bash" >> $out_script
echo "#$ -l h_vmem=${MEM}G" >> $out_script
echo 'set -e' >> $out_script
echo "" >> $out_script

if [[ $two_pass ]]
then
    two_pass='--two-pass'
fi

echo 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $out_script
echo "" >> $out_script

echo "singularity exec --bind /:/mnt docker://lethalfang/scalpel:0.5.3-2 bash -c \\" >> $out_script
echo "\"/opt/scalpel-0.5.3/scalpel-discovery --single \\" >> $out_script
echo "--ref /mnt/${HUMAN_REFERENCE} \\" >> $out_script
echo "--bed /mnt/${SELECTOR} \\" >> $out_script
echo "--bam /mnt/${tumor_bam} \\" >> $out_script
echo "--window 600 \\" >> $out_script
echo "${extra_discovery_arguments} \\" >> $out_script
echo "--dir /mnt/${outdir}/scalpel && \\" >> $out_script
echo "/opt/scalpel-0.5.3/scalpel-export --single \\" >> $out_script
echo "--db /mnt/${outdir}/scalpel/variants.db.dir \\" >> $out_script
echo "--ref /mnt/${HUMAN_REFERENCE} \\" >> $out_script
echo "--bed /mnt/${SELECTOR} \\" >> $out_script
echo "${extra_export_arguments} \\" >> $out_script
echo "> /mnt/${outdir}/scalpel/scalpel.vcf\"" >> $out_script
echo "" >> $out_script

echo "singularity exec --bind /:/mnt docker://lethalfang/scalpel:0.5.3-2 bash -c \\" >> $out_script
echo "\"cat /mnt/${outdir}/scalpel/scalpel.vcf | /opt/vcfsorter.pl /mnt/${HUMAN_REFERENCE%\.fa*}.dict - \\" >> $out_script
echo "> /mnt/${outdir}/${outvcf}\"" >> $out_script

echo "" >> $out_script
echo 'echo -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $out_script

${action} $out_script
