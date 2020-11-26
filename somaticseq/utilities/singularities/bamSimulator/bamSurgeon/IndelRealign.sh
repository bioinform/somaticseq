#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long output-dir:,tumor-bam:,normal-bam:,genome-reference:,selector:,out-tag:,extra-arguments:,out-script:,standalone, -n 'IndelRealign.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

#echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )

out_tag='JointRealigned'
#extra_arguments='-dt NONE --maxReadsForConsensuses 150000 --maxReadsInMemory 500000 --maxReadsForRealignment 2000000'

while true; do
    case "$1" in
        -o | --output-dir )
            case "$2" in
                "") shift 2 ;;
                *)  outdir=$2 ; shift 2 ;;
            esac ;;

        --tumor-bam )
            case "$2" in
                "") shift 2 ;;
                *)  tbam=$2 ; shift 2 ;;
            esac ;;

        --normal-bam )
            case "$2" in
                "") shift 2 ;;
                *)  nbam=$2 ; shift 2 ;;
            esac ;;

        --genome-reference )
            case "$2" in
                "") shift 2 ;;
                *)  HUMAN_REFERENCE=$2 ; shift 2 ;;
            esac ;;

        --selector )
            case "$2" in
                "") shift 2 ;;
                *)  SELECTOR=$2 ; shift 2 ;;
            esac ;;

        --out-tag )
            case "$2" in
                "") shift 2 ;;
                *)  out_tag=$2 ; shift 2 ;;
            esac ;;

        --extra-arguments )
            case "$2" in
                "") shift 2 ;;
                *)  extra_arguments=$2 ; shift 2 ;;
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
    out_script="${logdir}/indelRealign.${timestamp}.cmd"
fi

if [[ $standalone ]]
then
    echo "#!/bin/bash" > $out_script
    echo "" >> $out_script
    echo "#$ -o ${logdir}" >> $out_script
    echo "#$ -e ${logdir}" >> $out_script
    echo "#$ -S /bin/bash" >> $out_script
    echo '#$ -l h_vmem=10G' >> $out_script
    echo 'set -e' >> $out_script
fi

echo "" >> $out_script

if [[ $SELECTOR ]]
then
    selector_input="-L /mnt/${SELECTOR}"
fi

echo "singularity exec --bind /:/mnt   docker://broadinstitute/gatk3:3.8-0 java -Xmx9g -jar /usr/GenomeAnalysisTK.jar \\" >> $out_script
echo "-T RealignerTargetCreator \\" >> $out_script
echo "-R /mnt/${HUMAN_REFERENCE} \\" >> $out_script
echo "-I /mnt/${tbam} \\" >> $out_script
echo "-I /mnt/${nbam} \\" >> $out_script
echo "$selector_input \\" >> $out_script
echo "-o /mnt/${outdir}/T.N.intervals" >> $out_script
echo "" >> $out_script

echo "singularity exec --bind /:/mnt  --pwd /mnt/${outdir} docker://broadinstitute/gatk3:3.8-0 \\" >> $out_script
echo "java -Xmx9g -jar /usr/GenomeAnalysisTK.jar \\" >> $out_script
echo "-T IndelRealigner \\" >> $out_script
echo "-R /mnt/${HUMAN_REFERENCE} \\" >> $out_script
echo "-I /mnt/${tbam} \\" >> $out_script
echo "-I /mnt/${nbam} \\" >> $out_script
echo "-targetIntervals /mnt/${outdir}/T.N.intervals \\" >> $out_script
echo "${extra_arguments} \\" >> $out_script
echo "-nWayOut .${out_tag}.bam" >> $out_script
echo "" >> $out_script

realigned_normal=${nbam%.bam}.${out_tag}.bam
realigned_tumor=${tbam%.bam}.${out_tag}.bam

echo "mv ${realigned_normal%.bam}.bai ${realigned_normal}.bai" >> $out_script
echo "mv ${realigned_tumor%.bam}.bai ${realigned_tumor}.bai" >> $out_script
