#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long out-dir:,out-vcf:,tumor-bam:,normal-bam:,tumor-name:,normal-name:,human-reference:,selector:,dbsnp:,MEM:,threads:,extra-arguments:,action: -n 'submit_MuTect.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

#echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )
action=echo
MEM=7
threads=4

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

    --tumor-name )
        case "$2" in
            "") shift 2 ;;
            *)  tumor_name=$2 ; shift 2 ;;
        esac ;;

    --normal-name )
        case "$2" in
            "") shift 2 ;;
            *)  normal_name=$2 ; shift 2 ;;
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
            *)  dbsnp=$2 ; shift 2 ;;
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

out_script=${outdir}/logs/mutect2_${timestamp}.cmd

selector_text=''
if [[ -r $SELECTOR ]]; then
    selector_text="--intervals /mnt/${SELECTOR}"
fi

dbsnp_text=''
if [[ -r $dbsnp ]]; then
    dbsnp_text="--dbsnp /mnt/${dbsnp}"
fi


echo "#!/bin/bash" > $out_script
echo "" >> $out_script

echo "#$ -o ${logdir}" >> $out_script
echo "#$ -e ${logdir}" >> $out_script
echo "#$ -S /bin/bash" >> $out_script
echo "#$ -l h_vmem=${MEM}G" >> $out_script
echo "#$ -pe smp ${threads}" >> $out_script

echo 'set -e' >> $out_script
echo "" >> $out_script

echo 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $out_script
echo "" >> $out_script

if [[ $tumor_name ]]
then
    echo "tumor_name=${tumor_name}" >> $out_script
else
    echo "tumor_name=\`docker run --rm -v /:/mnt -u $UID --memory 1g lethalfang/samtools:0.1.19 samtools view -H /mnt/${tumor_bam} | egrep -w '^@RG' | grep -Po 'SM:[^\t$]+' | sed 's/SM://' | uniq | sed -e 's/[[:space:]]*$//'\`" >> $out_script
fi

if [[ $normal_name ]]
then
    echo "normal_name=${normal_name}" >> $out_script
else
    echo "normal_name=\`docker run --rm -v /:/mnt -u $UID --memory 1g lethalfang/samtools:0.1.19 samtools view -H /mnt/${normal_bam} | egrep -w '^@RG' | grep -Po 'SM:[^\t$]+' | sed 's/SM://' | uniq | sed -e 's/[[:space:]]*$//'\`" >> $out_script
    echo "" >> $out_script
fi

echo "" >> $out_script

echo "docker run --rm -v /:/mnt -u $UID --memory $(( MEM * threads ))G broadinstitute/gatk:4.beta.5 \\" >> $out_script
echo "java -Xmx6g -jar gatk.jar Mutect2 \\" >> $out_script
echo "--reference /mnt/${HUMAN_REFERENCE} \\" >> $out_script
echo "$selector_text \\" >> $out_script
echo "--input /mnt/${normal_bam} \\" >> $out_script
echo "--input /mnt/${tumor_bam} \\" >> $out_script
echo "--normalSampleName \${normal_name} \\" >> $out_script
echo "--tumorSampleName \${tumor_name} \\" >> $out_script
echo "$dbsnp_text \\" >> $out_script
echo "${extra_arguments} \\" >> $out_script
echo "--output /mnt/${outdir}/unfiltered.${outvcf}" >> $out_script
echo "" >> $out_script

echo "docker run --rm -v /:/mnt -u $UID --memory 6g broadinstitute/gatk:4.beta.5 \\" >> $out_script
echo "java -Xmx4g -jar gatk.jar FilterMutectCalls \\" >> $out_script
echo "--variant /mnt/${outdir}/unfiltered.${outvcf} \\" >> $out_script
echo "--output /mnt/${outdir}/${outvcf}" >> $out_script

echo "" >> $out_script
echo 'echo -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $out_script

${action} $out_script
