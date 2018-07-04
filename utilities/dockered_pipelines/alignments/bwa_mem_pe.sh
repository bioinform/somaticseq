#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long output-dir:,fq1:,fq2:,ID:,LB:,PL:,SM:,bam-header:,genome-reference:,out-script:,out-bam:,threads:,standalone, -n 'bwa_mem_pe.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

#echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )

ID='unknownIdentifier'
LB='unknown'
PL='illumina'
SM='Alignment'
outBam='Alignment.bam'
threads=1

while true; do
    case "$1" in
        -o | --output-dir )
            case "$2" in
                "") shift 2 ;;
                *)  outdir=$2 ; shift 2 ;;
            esac ;;

        --fq1 )
            case "$2" in
                "") shift 2 ;;
                *)  fq1=$2 ; shift 2 ;;
            esac ;;

        --fq2 )
            case "$2" in
                "") shift 2 ;;
                *)  fq2=$2 ; shift 2 ;;
            esac ;;

        --out-bam )
            case "$2" in
                "") shift 2 ;;
                *)  outBam=$2 ; shift 2 ;;
            esac ;;

        --bam-header )
            case "$2" in
                "") shift 2 ;;
                *)  bamHeader=$2 ; shift 2 ;;
            esac ;;

        --ID )
            case "$2" in
                "") shift 2 ;;
                *)  ID=$2 ; shift 2 ;;
            esac ;;

        --LB )
            case "$2" in
                "") shift 2 ;;
                *)  LB=$2 ; shift 2 ;;
            esac ;;

        --PL )
            case "$2" in
                "") shift 2 ;;
                *)  PL=$2 ; shift 2 ;;
            esac ;;

        --SM )
            case "$2" in
                "") shift 2 ;;
                *)  SM=$2 ; shift 2 ;;
            esac ;;

        --genome-reference )
            case "$2" in
                "") shift 2 ;;
                *)  HUMAN_REFERENCE=$2 ; shift 2 ;;
            esac ;;

        --threads )
            case "$2" in
                "") shift 2 ;;
                *)  threads=$2 ; shift 2 ;;
            esac ;;

        --out-script )
            case "$2" in
                "") shift 2 ;;
                *)  out_script_name=$2 ; shift 2 ;;
            esac ;;

        --standalone )
            standalone=1 ; shift ;;

        -- ) shift; break ;;
        * ) break ;;
    esac
done

logdir=${outdir}/logs
mkdir -p ${logdir}

if [[ ${out_script_name} ]]
then
    out_script="${out_script_name}"
else
    out_script="${logdir}/bwa.mem.${timestamp}.cmd"
fi

if [[ $standalone ]]
then
    echo "#!/bin/bash" > $out_script
    echo "" >> $out_script
    echo "#$ -o ${logdir}" >> $out_script
    echo "#$ -e ${logdir}" >> $out_script
    echo "#$ -S /bin/bash" >> $out_script
    echo '#$ -l h_vmem=16G' >> $out_script
    echo "#$ -pe smp ${threads}" >> $out_script
    echo 'set -e' >> $out_script
fi

echo "" >> $out_script

if [[ ! $bamHeader ]]
then
    bamHeader="@RG\tID:${ID}\tLB:${LB}\tPL:${PL}\tSM:${SM}"
fi

echo "docker run --rm -v /:/mnt -u $UID lethalfang/bwa:0.7.17_samtools bash -c \\" >> $out_script
echo "\"bwa mem \\" >> $out_script
echo "-R '${bamHeader}' \\" >> $out_script
echo "-M -t ${threads} \\" >> $out_script
echo "/mnt/${HUMAN_REFERENCE} \\" >> $out_script
echo "/mnt/${fq1} \\" >> $out_script
echo "/mnt/${fq2} \\" >> $out_script
echo "| samtools view -Sbh - \\" >> $out_script
echo "| samtools sort -m 4G --threads ${threads} -o /mnt/${outdir}/${outBam}\"" >> $out_script

echo "" >> $out_script

echo "docker run --rm -v /:/mnt -u $UID lethalfang/bwa:0.7.17_samtools \\" >> $out_script
echo "samtools index /mnt/${outdir}/${outBam}" >> $out_script
