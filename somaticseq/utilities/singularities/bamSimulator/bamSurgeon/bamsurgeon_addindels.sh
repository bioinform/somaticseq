#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long output-dir:,genome-reference:,bam-out:,bam-in:,indels:,cnv-file:,min-vaf:,max-vaf:,min-depth:,max-depth:,min-variant-reads:,aligner:,out-script:,seed:,standalone -n 'bamsurgeon_addindels.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

#echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )

seed=$( date +"%Y" )
min_depth=10
max_dpeth=5000
min_var_reads=1
aligner='mem'

while true; do
    case "$1" in
        -o | --output-dir )
            case "$2" in
                "") shift 2 ;;
                *)  outdir=$2 ; shift 2 ;;
            esac ;;
            
        --genome-reference )
            case "$2" in
                "") shift 2 ;;
                *)  HUMAN_REFERENCE=$2 ; shift 2 ;;
            esac ;;

        --selector )
            case "$2" in
                "") shift 2 ;;
                *) SELECTOR=$2 ; shift 2 ;;
            esac ;;

        --bam-out )
            case "$2" in
                "") shift 2 ;;
                *)  outbam=$2 ; shift 2 ;;
            esac ;;

        --bam-in )
            case "$2" in
                "") shift 2 ;;
                *)  inbam=$2 ; shift 2 ;;
            esac ;;

        --indels )
            case "$2" in
                "") shift 2 ;;
                *)  indels=$2 ; shift 2 ;;
            esac ;;

        --cnv-file )
            case "$2" in
                "") shift 2 ;;
                *)  cnvfile=$2 ; shift 2 ;;
            esac ;;

        --min-vaf )
            case "$2" in
                "") shift 2 ;;
                *)  min_vaf=$2 ; shift 2 ;;
            esac ;;

        --max-vaf )
            case "$2" in
                "") shift 2 ;;
                *)  max_vaf=$2 ; shift 2 ;;
            esac ;;

        --min-depth )
            case "$2" in
                "") shift 2 ;;
                *)  min_depth=$2 ; shift 2 ;;
            esac ;;

        --max-depth )
            case "$2" in
                "") shift 2 ;;
                *)  max_depth=$2 ; shift 2 ;;
            esac ;;

        --min-variant-reads )
            case "$2" in
                "") shift 2 ;;
                *)  min_var_reads=$2 ; shift 2 ;;
            esac ;;

        --out-script )
            case "$2" in
                "") shift 2 ;;
                *)  out_script_name=$2 ; shift 2 ;;
            esac ;;

        --aligner )
            case "$2" in
                "") shift 2 ;;
                *)  aligner=$2 ; shift 2 ;;
            esac ;;

        --seed )
            case "$2" in
                "") shift 2 ;;
                *)  seed=$2 ; shift 2 ;;
            esac ;;

        --standalone )
            standalone=1 ; shift ;;

        -- ) shift; break ;;
        * ) break ;;
    esac
done

hg_dict=${HUMAN_REFERENCE%\.fa*}.dict

logdir=${outdir}/logs
mkdir -p ${logdir}

if [[ ${out_script_name} ]]
then
    out_script="${out_script_name}"
else
    out_script="${logdir}/addindels.${timestamp}.cmd"    
fi

if [[ $standalone ]]
then
    echo "#!/bin/bash" > $out_script
    echo "" >> $out_script
    echo "#$ -o ${logdir}" >> $out_script
    echo "#$ -e ${logdir}" >> $out_script
    echo "#$ -S /bin/bash" >> $out_script
    echo '#$ -l h_vmem=8G' >> $out_script
    echo 'set -e' >> $out_script
fi

echo "" >> $out_script

echo "singularity exec --bind /:/mnt  --pwd /mnt/${outdir} docker://lethalfang/bamsurgeon:1.1-3 \\" >> $out_script
echo "/usr/local/bamsurgeon/bin/addindel.py \\" >> $out_script
echo "--snvfrac 0.1 --mutfrac 0.5 --coverdiff 0.9 --procs 1 \\" >> $out_script
echo "--varfile /mnt/${indels} \\" >> $out_script
echo "--bamfile /mnt/${inbam} \\" >> $out_script
echo "--reference /mnt/${HUMAN_REFERENCE} \\" >> $out_script
echo "--cnvfile /mnt/${cnvfile} \\" >> $out_script
echo "--outbam /mnt/${outdir}/unsorted.${outbam} \\" >> $out_script
echo "--picardjar /usr/local/bin/picard.jar \\" >> $out_script
echo "--mindepth $min_depth \\" >> $out_script
echo "--maxdepth $max_depth \\" >> $out_script
echo "--minmutreads $min_var_reads \\" >> $out_script
echo "--seed $seed \\" >> $out_script
echo "--tagreads --force \\" >> $out_script
echo "--aligner "${aligner}"" >> $out_script
echo "" >> $out_script

echo "singularity exec --bind /:/mnt   docker://lethalfang/bamsurgeon:1.1-3 bash -c \\" >> $out_script
echo "\"/usr/local/bamsurgeon/scripts/makevcf_indels.py \\" >> $out_script
echo "/mnt/${outdir}/addindel_logs_unsorted.${outbam} /mnt/${HUMAN_REFERENCE} \\" >> $out_script
echo "| bedtools sort -header -faidx /mnt/${HUMAN_REFERENCE}.fai \\" >> $out_script
echo "> /mnt/${outdir}/synthetic_indels.vcf\"" >> $out_script
echo "" >> $out_script

echo "singularity exec --bind /:/mnt   docker://lethalfang/somaticseq:base-1.1 bash -c \\" >> $out_script
echo "\"java -jar /opt/GATK/GenomeAnalysisTK.jar -T LeftAlignAndTrimVariants \\" >> $out_script
echo "-R /mnt/${HUMAN_REFERENCE} \\" >> $out_script
echo "--variant /mnt/${outdir}/synthetic_indels.vcf \\" >> $out_script
echo "| egrep -v '^[0-9]+ variants|^INFO' > /mnt/${outdir}/synthetic_indels.leftAlign.vcf\"" >> $out_script
echo "" >> $out_script

echo "singularity exec --bind /:/mnt   docker://lethalfang/samtools:1.7 \\" >> $out_script
echo "samtools sort -m 4G --reference /mnt/${HUMAN_REFERENCE} -o /mnt/${outdir}/${outbam} /mnt/${outdir}/unsorted.${outbam}" >> $out_script

echo "singularity exec --bind /:/mnt   docker://lethalfang/samtools:1.7 samtools index /mnt/${outdir}/${outbam}" >> $out_script

echo "" >> $out_script
echo "rm ${outdir}/unsorted.${outbam}" >> $out_script
