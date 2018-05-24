#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long output-dir:,bam:,out-vcf:,genome-reference:,dbsnp:,hapmap:,omni:,thousandG:,mills:,out-script:,action:,threads:, -n 'bam2vcf.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )

tumor_bam_header='@RG\tID:myPipeline\tLB:myLibrary\tPL:illumina\tSM:TUMOR'
normal_bam_header='@RG\tID:myPipeline\tLB:myLibrary\tPL:illumina\tSM:NORMAL'
MEM=16
threads=24
action=echo

while true; do
    case "$1" in
    
    -o | --output-dir )
        case "$2" in
            "") shift 2 ;;
            *)  outdir=$2 ; shift 2 ;;
        esac ;;

    --bam )
        case "$2" in
            "") shift 2 ;;
            *)  bam=$2 ; shift 2 ;;
        esac ;;

    --out-vcf )
        case "$2" in
            "") shift 2 ;;
            *)  outVcf=$2 ; shift 2 ;;
        esac ;;


    --genome-reference )
        case "$2" in
            "") shift 2 ;;
            *)  GENOME_REFERENCE=$2 ; shift 2 ;;
        esac ;;

    --dbsnp )
        case "$2" in
            "") shift 2 ;;
            *)  dbsnp=$2 ; shift 2 ;;
        esac ;;

    --hapmap )
        case "$2" in
            "") shift 2 ;;
            *)  hapmapFile=$2 ; shift 2 ;;
        esac ;;

    --thousandG )
        case "$2" in
            "") shift 2 ;;
            *)  thousandGFile=$2 ; shift 2 ;;
        esac ;;

    --omni )
        case "$2" in
            "") shift 2 ;;
            *)  omniFile=$2 ; shift 2 ;;
        esac ;;

    --mills )
        case "$2" in
            "") shift 2 ;;
            *)  millsFile=$2 ; shift 2 ;;
        esac ;;

    --threads )
        case "$2" in
            "") shift 2 ;;
            *)  threads=$2 ; shift 2 ;;
        esac ;;

    --MEM )
        case "$2" in
            "") shift 2 ;;
            *)  MEM=$2 ; shift 2 ;;
        esac ;;

    --out-script )
        case "$2" in
            "") shift 2 ;;
            *)  out_script_name=$2 ; shift 2 ;;
        esac ;;

    --action )
        case "$2" in
            "") shift 2 ;;
            *)  action=$2 ; shift 2 ;;
        esac ;;

    -- ) shift; break ;;
    * ) break ;;
    
    esac
done



logdir=${outdir}/logs
mkdir -p ${logdir}


if [[ ${out_script_name} ]]
then
    out_script="${logdir}/${out_script_name}"
else
    out_script="${logdir}/bam2vcf.${timestamp}.cmd"    
fi


echo "#!/bin/bash" > $out_script
echo "" >> $out_script

echo "#$ -o ${logdir}" >> $out_script
echo "#$ -e ${logdir}" >> $out_script
echo "#$ -S /bin/bash" >> $out_script
echo "#$ -l h_vmem=6G" >> $out_script
echo "#$ -pe smp ${threads}" >> $out_script

echo 'set -e' >> $out_script
echo "" >> $out_script

files_to_delete=''

echo 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $out_script
echo "" >> $out_script



$MYDIR/haplotypeCaller.sh \
--out-dir         ${outdir} \
--bam             ${bam} \
--human-reference ${GENOME_REFERENCE} \
--dbsnp           ${dbsnp} \
--out-vcf         preVQSR.${outVcf} \
--threads         ${threads} \
--MEM             ${MEM} \
--out-script      ${out_script}

$MYDIR/VQSR.sh \
--out-dir         ${outdir} \
--in-vcf          ${outdir}/preVQSR.${outVcf} \
--human-reference ${GENOME_REFERENCE} \
--dbsnp           ${dbsnp} \
--hapmap          ${hapmapFile} \
--omni            ${omniFile} \
--thousandG       ${thousandGFile} \
--mills           ${millsFile} \
--out-vcf         ${outVcf} \
--out-script      ${out_script}
