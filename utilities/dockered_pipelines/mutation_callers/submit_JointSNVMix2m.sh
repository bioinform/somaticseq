#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long out-dir:,out-vcf:,tumor-bam:,normal-bam:,human-reference:,action:,skip-size:,convergence-threshold:,split: -n 'submit_JointSNVMix2.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

#echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )
convergence_threshold=0.01
skip_size=10000
action=echo

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

    --skip-size )
        case "$2" in
            "") shift 2 ;;
            *) skip_size=$2 ; shift 2 ;;
        esac ;;

    --convergence-threshold )
        case "$2" in
            "") shift 2 ;;
            *) convergence_threshold=$2 ; shift 2 ;;
        esac ;;

    --split )
        case "$2" in
            "") shift 2 ;;
            *)  split=$2 ; shift 2 ;;
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

logdir=${outdir}/logs
mkdir -p ${logdir}

jsm_script=${outdir}/logs/jsm2_${timestamp}.cmd

echo "#!/bin/bash" > $jsm_script
echo "" >> $jsm_script

echo "#$ -o ${logdir}" > $jsm_script
echo "#$ -e ${logdir}" >> $jsm_script
echo "#$ -S /bin/bash" >> $jsm_script
echo '#$ -l h_vmem=8G' >> $jsm_script
echo 'set -e' >> $jsm_script
echo "" >> $jsm_script

echo 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $jsm_script
echo "" >> $jsm_script

echo "docker run --rm -v /:/mnt -u $UID lethalfang/jointsnvmix2:0.7.5 \\" >> $jsm_script
echo "/opt/JointSNVMix-0.7.5/build/scripts-2.7/jsm.py train joint_snv_mix_two \\" >> $jsm_script
echo "--convergence_threshold $convergence_threshold \\" >> $jsm_script
echo "--skip_size $skip_size \\" >> $jsm_script
echo "/mnt/${HUMAN_REFERENCE} \\" >> $jsm_script
echo "/mnt/${normal_bam} \\" >> $jsm_script
echo "/mnt/${tumor_bam} \\" >> $jsm_script
echo "/opt/JointSNVMix-0.7.5/config/joint_priors.cfg \\" >> $jsm_script
echo "/opt/JointSNVMix-0.7.5/config/joint_params.cfg \\" >> $jsm_script
echo "/mnt/${outdir}/jsm.parameter.cfg" >> $jsm_script
echo "" >> $jsm_script

echo "echo -e '##fileformat=VCFv4.1' > ${outdir}/JointSNVMix2.vcf" >> $jsm_script
echo "echo -e '##INFO=<ID=AAAB,Number=1,Type=Float,Description=\"Probability of Joint Genotype AA in Normal and AB in Tumor\">' >> ${outdir}/${outvcf}" >> $jsm_script
echo "echo -e '##INFO=<ID=AABB,Number=1,Type=Float,Description=\"Probability of Joint Genotype AA in Normal and BB in Tumor\">' >> ${outdir}/${outvcf}" >> $jsm_script
echo "echo -e '##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases (reads1)\">' >> ${outdir}/${outvcf}" >> $jsm_script
echo "echo -e '##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases (reads2)\">' >> ${outdir}/${outvcf}" >> $jsm_script
echo "echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR' >> ${outdir}/${outvcf}" >> $jsm_script
echo "" >> $jsm_script

echo "docker run --rm -v /:/mnt -u $UID lethalfang/jointsnvmix2:0.7.5 bash -c \\"  >> $jsm_script
echo \""/opt/JointSNVMix-0.7.5/build/scripts-2.7/jsm.py classify joint_snv_mix_two \\" >> $jsm_script
echo "/mnt/${HUMAN_REFERENCE} \\" >> $jsm_script
echo "/mnt/${normal_bam} \\" >> $jsm_script
echo "/mnt/${tumor_bam} \\" >> $jsm_script
echo "/mnt/${outdir}/jsm.parameter.cfg \\" >> $jsm_script
echo "/dev/stdout | awk -F \\\"\t\\\" 'NR!=1 && \\\$4!=\\\"N\\\" && \\\$10+\\\$11>=0.95' | \\" >> $jsm_script
echo "awk -F \\\"\t\\\" '{print \\\$1 \\\"\t\\\" \\\$2 \\\"\t.\t\\\" \\\$3 \\\"\t\\\" \\\$4 \\\"\t.\t.\tAAAB=\\\" \\\$10 \\\";AABB=\\\" \\\$11 \\\"\tRD:AD\t\\\" \\\$5 \\\":\\\" \\\$6 \\\"\t\\\" \\\$7 \\\":\\\" \\\$8}' \\" >> $jsm_script
echo "| /opt/vcfsorter.pl /mnt/${HUMAN_REFERENCE%\.fa*}.dict - >> /mnt/${outdir}/${outvcf}\"" >> $jsm_script
echo "" >> $jsm_script

if [[ $split ]]
then
    echo "i=1" >> $jsm_script
    echo "while [[ \$i -le $split ]]" >> $jsm_script
    echo "do" >> $jsm_script
    echo "    docker run --rm -v /:/mnt -u $UID lethalfang/somaticseq:base-1.0 bedtools intersect -a /mnt/${outdir}/${outvcf} -b /mnt/${outdir}/\${i}/\${i}.bed -header | uniq >  ${outdir}/\${i}/${outvcf}" >> $jsm_script
    echo "    i=\$(( \$i + 1 ))" >> $jsm_script
    echo "done" >> $jsm_script
fi

echo "" >> $jsm_script
echo 'echo -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $jsm_script

${action} $jsm_script
