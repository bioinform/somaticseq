#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long output-dir:,genome-reference:,bam-out:,bam-in:,svs:,cnv-file:,min-vaf:,max-vaf:,min-depth:,max-depth:,min-variant-reads:,aligner:,out-script:,seed:,standalone -n 'bamsurgeon_addsvs.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

#echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )
seed=$( date +"%Y" )
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

        --svs )
            case "$2" in
                "") shift 2 ;;
                *)  svs=$2 ; shift 2 ;;
            esac ;;

        --cnv-file )
            case "$2" in
                "") shift 2 ;;
                *)  cnvfile=$2 ; shift 2 ;;
            esac ;;

        --aligner )
            case "$2" in
                "") shift 2 ;;
                *)  aligner=$2 ; shift 2 ;;
            esac ;;

        --out-script )
            case "$2" in
                "") shift 2 ;;
                *)  out_script_name=$2 ; shift 2 ;;
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
    out_script="${logdir}/addsvs.${timestamp}.cmd"    
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

echo "docker run -v /:/mnt -u $UID --rm --memory 14g --workdir=/mnt/${outdir} lethalfang/bamsurgeon:1.1-3 \\" >> $out_script
echo "/usr/local/bamsurgeon/bin/addsv.py \\" >> $out_script
echo "--svfrac 0.1 --procs 1 \\" >> $out_script
echo "--varfile /mnt/${svs} \\" >> $out_script
echo "--bamfile /mnt/${inbam} \\" >> $out_script
echo "--reference /mnt/${HUMAN_REFERENCE} \\" >> $out_script
echo "--outbam /mnt/${outdir}/unsorted.${outbam} \\" >> $out_script
echo "--maxlibsize 1000 \\" >> $out_script
echo "--cnvfile /mnt/${cnvfile} \\" >> $out_script
echo "--tagreads \\" >> $out_script
echo "--seed $seed \\" >> $out_script
echo "--aligner "${aligner}"" >> $out_script
echo "" >> $out_script


echo "docker run -v /:/mnt -u $UID --rm lethalfang/bamsurgeon:1.1-3 bash -c \\" >> $out_script
echo "\"/usr/local/bamsurgeon/scripts/makevcf_sv.py -l /mnt/${outdir}/addsv_logs_unsorted.${outbam} \\" >> $out_script
echo "-r /mnt/${HUMAN_REFERENCE} \\" >> $out_script
echo "| bedtools sort -header -faidx /mnt/${HUMAN_REFERENCE}.fai \\" >> $out_script
echo "> /mnt/${outdir}/synthetic_svs.vcf\"" >> $out_script
echo "" >> $out_script


echo "docker run -v /:/mnt -u $UID --rm lethalfang/samtools:1.7 \\" >> $out_script
echo "samtools sort -m 4G --reference /mnt/${HUMAN_REFERENCE} -o /mnt/${outdir}/${outbam} /mnt/${outdir}/unsorted.${outbam}" >> $out_script

echo "docker run -v /:/mnt -u $UID --rm lethalfang/samtools:1.7 samtools index /mnt/${outdir}/${outbam}" >> $out_script

echo "" >> $out_script
echo "rm ${outdir}/unsorted.${outbam}" >> $out_script
