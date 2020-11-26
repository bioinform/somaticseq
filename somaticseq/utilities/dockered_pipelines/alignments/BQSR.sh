#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long output-dir:,in-bam:,out-bam:,genome-reference:,dbsnp:,out-script:,standalone,plot, -n 'BQSR.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

#echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )

while true; do
    case "$1" in
        -o | --output-dir )
            case "$2" in
                "") shift 2 ;;
                *)  outdir=$2 ; shift 2 ;;
            esac ;;
            
        --in-bam )
            case "$2" in
                "") shift 2 ;;
                *)  inBam=$2 ; shift 2 ;;
            esac ;;

        --out-bam )
            case "$2" in
                "") shift 2 ;;
                *)  outBam=$2 ; shift 2 ;;
            esac ;;

        --genome-reference )
            case "$2" in
                "") shift 2 ;;
                *)  HUMAN_REFERENCE=$2 ; shift 2 ;;
            esac ;;

        --dbsnp )
            case "$2" in
                "") shift 2 ;;
                *) dbsnp=$2 ; shift 2 ;;
            esac ;;

        --out-script )
            case "$2" in
                "") shift 2 ;;
                *)  out_script_name=$2 ; shift 2 ;;
            esac ;;

        --standalone )
            standalone=1 ; shift ;;

        --plot )
            plotBQSR=1 ; shift ;;

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
    out_script="${logdir}/BQSR.${timestamp}.cmd"    
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

echo "docker run --rm -v /:/mnt -u $UID broadinstitute/gatk3:3.8-1 \\" >> $out_script
echo "java -Xmx8g -jar /usr/GenomeAnalysisTK.jar \\" >> $out_script
echo "-T BaseRecalibrator \\" >> $out_script
echo "-R /mnt/${HUMAN_REFERENCE} \\" >> $out_script
echo "-I /mnt/${inBam} \\" >> $out_script
echo "-knownSites /mnt/${dbsnp} \\" >> $out_script
echo "-o /mnt/${outdir}/BQSR.${timestamp}.table" >> $out_script

echo "" >> $out_script

if [[ $plotBQSR ]]
then
    echo "docker run --rm -v /:/mnt -u $UID broadinstitute/gatk3:3.8-1 \\" >> $out_script
    echo "java -Xmx8g -jar /usr/GenomeAnalysisTK.jar \\" >> $out_script
    echo "-T BaseRecalibrator \\" >> $out_script
    echo "-R /mnt/${HUMAN_REFERENCE} \\" >> $out_script
    echo "-I /mnt/${inBam} \\" >> $out_script
    echo "-knownSites /mnt/${dbsnp} \\" >> $out_script
    echo "-BQSR /mnt/${outdir}/BQSR.${timestamp}.table \\" >> $out_script
    echo "-o /mnt/${outdir}/post_BQSR.${timestamp}.table" >> $out_script
    
    echo "" >> $out_script
    
    echo "# The GATK docker doesn't have all the R packages to do the plotting. Do it yourself locally: install R packages of ggplot2, gplots, reshape, grid, gsalib" >> $out_script
    echo "#java -Xmx8g -jar /usr/GenomeAnalysisTK.jar \\" >> $out_script
    echo "#-T AnalyzeCovariates \\" >> $out_script
    echo "#-R ${HUMAN_REFERENCE} \\" >> $out_script
    echo "#-before ${outdir}/BQSR.${timestamp}.table \\" >> $out_script
    echo "#-after ${outdir}/post_BQSR.${timestamp}.table \\" >> $out_script
    echo "#-plots ${outdir}/BQSR.plot.${timestamp}.pdf" >> $out_script
    
    echo "" >> $out_script
fi

echo "docker run --rm -v /:/mnt -u $UID broadinstitute/gatk3:3.8-1 \\" >> $out_script
echo "java -Xmx8g -jar /usr/GenomeAnalysisTK.jar \\" >> $out_script
echo "-T PrintReads \\" >> $out_script
echo "-R /mnt/${HUMAN_REFERENCE} \\" >> $out_script
echo "-I /mnt/${inBam} \\" >> $out_script
echo "-BQSR /mnt/${outdir}/BQSR.${timestamp}.table \\" >> $out_script
echo "-o /mnt/${outdir}/${outBam}" >> $out_script

echo "" >> $out_script

echo "mv ${outdir}/${outBam%.bam}.bai ${outdir}/${outBam}.bai" >> $out_script
