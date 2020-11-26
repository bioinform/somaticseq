#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long out-dir:,out-vcf:,in-vcf:,human-reference:,selector:,dbsnp:,hapmap:,omni:,thousandG:,mills:,extra-arguments:,action:,MEM:,threads:,out-script:,standalone -n 'gatk_vqsr.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

#echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )
action=echo
MEM=12

while true; do
    case "$1" in
    
    -o | --out-dir )
        case "$2" in
            "") shift 2 ;;
            *)  outdir=$2 ; shift 2 ;;
        esac ;;
    
    --out-vcf )
        case "$2" in
            "") shift 2 ;;
            *)  outVcfName=$2 ; shift 2 ;;
        esac ;;

    --in-vcf )
        case "$2" in
            "") shift 2 ;;
            *)  inputVcfFile=$2 ; shift 2 ;;
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
            *)  dbsnpFile=$2 ; shift 2 ;;
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

    --MEM )
        case "$2" in
            "") shift 2 ;;
            *)  MEM=$2 ; shift 2 ;;
        esac ;;

    --threads )
        case "$2" in
            "") shift 2 ;;
            *)  threads=$2 ; shift 2 ;;
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

    --action )
        case "$2" in
            "") shift 2 ;;
            *)  action=$2 ; shift 2 ;;
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
    out_script="${logdir}/vqsr.${timestamp}.cmd"
fi


if [[ $standalone ]]
then
    echo "#!/bin/bash" > $out_script
    echo "" >> $out_script
    echo "#$ -o ${logdir}" >> $out_script
    echo "#$ -e ${logdir}" >> $out_script
    echo "#$ -S /bin/bash" >> $out_script
    echo '#$ -l h_vmem=14G' >> $out_script
    echo 'set -e' >> $out_script
fi

echo "" >> $out_script
echo 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $out_script
echo "" >> $out_script


if [[ ${SELECTOR} ]]
then
    selector_text="-L /mnt/${SELECTOR}"
fi

dbsnp_text=''
if [[ ${dbsnp} ]]; then
    dbsnp_text="--dbsnp /mnt/${dbsnp}"
fi


# VariantRecalibrator for SNP
echo "docker run --rm -v /:/mnt -u $UID broadinstitute/gatk:4.0.5.2 \\" >> $out_script
echo "java -Xmx${MEM}g -jar /gatk/gatk.jar VariantRecalibrator \\" >> $out_script
echo "--variant /mnt/${inputVcfFile} \\" >> $out_script
echo "--reference /mnt/${HUMAN_REFERENCE} \\" >> $out_script
echo "--output /mnt/${outdir}/snp.${timestamp}.recal \\" >> $out_script
echo "--tranches-file /mnt/${outdir}/snp.${timestamp}.tranches \\" >> $out_script
echo "-an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \\" >> $out_script
echo "-mode SNP \\" >> $out_script
echo "-resource hapmap,known=false,training=true,truth=true,prior=15.0:/mnt/${hapmapFile} \\" >> $out_script
echo "-resource omni,known=false,training=true,truth=true,prior=12.0:/mnt/${omniFile} \\" >> $out_script
echo "-resource 1000G,known=false,training=true,truth=false,prior=10.0:/mnt/${thousandGFile} \\" >> $out_script
echo "-resource dbsnp,known=true,training=false,truth=false,prior=2.0:/mnt/${dbsnpFile}" >> $out_script
echo "" >> $out_script

# VariantRecalibrator for INDEL
echo "docker run --rm -v /:/mnt -u $UID broadinstitute/gatk:4.0.5.2 \\" >> $out_script
echo "java -Xmx${MEM}g -jar /gatk/gatk.jar VariantRecalibrator \\" >> $out_script
echo "--variant /mnt/${inputVcfFile} \\" >> $out_script
echo "--reference /mnt/${HUMAN_REFERENCE} \\" >> $out_script
echo "--output /mnt/${outdir}/indel.${timestamp}.recal \\" >> $out_script
echo "--tranches-file /mnt/${outdir}/indel.${timestamp}.tranches \\" >> $out_script
echo "-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum \\" >> $out_script
echo "-mode INDEL \\" >> $out_script
echo "-resource mills,known=false,training=true,truth=true,prior=12.0:/mnt/${millsFile} \\" >> $out_script
echo "-resource dbsnp,known=true,training=false,truth=false,prior=2.0:/mnt/${dbsnpFile}" >> $out_script
echo "" >> $out_script


# Apply VQSR for SNP
echo "docker run --rm -v /:/mnt -u $UID broadinstitute/gatk:4.0.5.2 \\" >> $out_script
echo "java -Xmx${MEM}g -jar /gatk/gatk.jar ApplyVQSR \\" >> $out_script
echo "--reference /mnt/${HUMAN_REFERENCE} \\" >> $out_script
echo "--recal-file /mnt/${outdir}/snp.${timestamp}.recal \\" >> $out_script
echo "--tranches-file /mnt/${outdir}/snp.${timestamp}.tranches \\" >> $out_script
echo "--variant /mnt/${inputVcfFile} \\" >> $out_script
echo "--output /mnt/${outdir}/snp.vqsr.${timestamp}.vcf \\" >> $out_script
echo "-mode SNP" >> $out_script
echo "" >> $out_script

# Apply VQSR for INDEL
echo "docker run --rm -v /:/mnt -u $UID broadinstitute/gatk:4.0.5.2 \\" >> $out_script
echo "java -Xmx${MEM}g -jar /gatk/gatk.jar ApplyVQSR \\" >> $out_script
echo "--reference /mnt/${HUMAN_REFERENCE} \\" >> $out_script
echo "--recal-file /mnt/${outdir}/indel.${timestamp}.recal \\" >> $out_script
echo "--tranches-file /mnt/${outdir}/indel.${timestamp}.tranches \\" >> $out_script
echo "--variant /mnt/${inputVcfFile} \\" >> $out_script
echo "--output /mnt/${outdir}/indel.vqsr.${timestamp}.vcf \\" >> $out_script
echo "-mode INDEL" >> $out_script
echo "" >> $out_script


# Combine SNP and INDEL VCF files:
echo "docker run --rm -v /:/mnt -u $UID lethalfang/vcftools:0.1.15 bash -c \\" >> $out_script
echo "\"vcf-concat /mnt/${outdir}/snp.vqsr.${timestamp}.vcf /mnt/${outdir}/indel.vqsr.${timestamp}.vcf | egrep '^#|VQSLOD' | perl /opt/vcfsorter.pl /mnt/${HUMAN_REFERENCE%\.fa*}.dict - > /mnt/${outdir}/${outVcfName}\"" >> $out_script


echo '' >> $out_script
echo "for file in ${outdir}/snp.${timestamp}.recal ${outdir}/snp.${timestamp}.tranches ${outdir}/indel.${timestamp}.recal ${outdir}/indel.${timestamp}.tranches ${outdir}/snp.vqsr.${timestamp}.vcf ${outdir}/indel.vqsr.${timestamp}.vcf ${outdir}/snp.${timestamp}.recal.idx ${outdir}/indel.${timestamp}.recal.idx ${outdir}/snp.vqsr.${timestamp}.vcf.idx ${outdir}/indel.vqsr.${timestamp}.vcf.idx" >> $out_script
echo '    do rm -fv $file' >> $out_script
echo "done" >> $out_script


echo 'echo -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $out_script
