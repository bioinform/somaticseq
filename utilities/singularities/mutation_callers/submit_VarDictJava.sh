#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long out-dir:,out-vcf:,tumor-bam:,normal-bam:,human-reference:,selector:,extra-arguments:,action:,VAF:,MEM: -n 'submit_VarDictJava.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

#echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )
VAF=0.05
action=echo
MEM=8

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

VERSION=`cat ${MYDIR}/../../../VERSION | sed 's/##SomaticSeq=v//'`

logdir=${outdir}/logs
mkdir -p ${logdir}

out_script=${outdir}/logs/vardict_${timestamp}.cmd

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


total_bases=`cat ${SELECTOR} | awk -F "\t" '{print $3-$2}' | awk '{ sum += $1 } END { print sum }'`
num_lines=`cat ${SELECTOR} | wc -l`

input_bed=${SELECTOR}
if [[ $(( $total_bases / $num_lines )) -gt 50000 ]]
then
    echo "singularity exec --bind /:/mnt docker://lethalfang/somaticseq:${VERSION} \\" >> $out_script
    echo "/opt/somaticseq/utilities/split_mergedBed.py \\" >> $out_script
    echo "-infile /mnt/${SELECTOR} -outfile /mnt/${outdir}/split_regions.bed" >> $out_script
    echo "" >> $out_script
    
    input_bed="${outdir}/split_regions.bed"
fi


echo "singularity exec --bind /:/mnt docker://lethalfang/vardictjava:1.5.1 bash -c \\" >> $out_script
echo "\"/opt/VarDict-1.5.1/bin/VarDict \\" >> $out_script
echo "${extra_arguments} \\" >> $out_script
echo "-G /mnt/${HUMAN_REFERENCE} \\" >> $out_script
echo "-f $VAF -h \\" >> $out_script
echo "-b '/mnt/${tumor_bam}|/mnt/${normal_bam}' \\" >> $out_script
echo "-Q 1 -c 1 -S 2 -E 3 -g 4 /mnt/${input_bed} \\" >> $out_script
echo "> /mnt/${outdir}/${timestamp}.var\"" >> $out_script
echo "" >> $out_script

echo "singularity exec --bind /:/mnt docker://lethalfang/vardictjava:1.5.1 \\" >> $out_script
echo "bash -c \"cat /mnt/${outdir}/${timestamp}.var | awk 'NR!=1' | /opt/VarDict/testsomatic.R | /opt/VarDict/var2vcf_paired.pl -N 'TUMOR|NORMAL' -f $VAF \\" >> $out_script
echo "> /mnt/${outdir}/${outvcf}\"" >> $out_script

echo "" >> $out_script
echo 'echo -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $out_script

${action} $out_script
