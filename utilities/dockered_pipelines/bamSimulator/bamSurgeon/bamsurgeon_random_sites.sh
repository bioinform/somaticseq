#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long output-dir:,genome-reference:,selector:,num-snvs:,num-indels:,num-svs:,max-indel-len:,min-vaf:,max-vaf:,left-beta:,right-beta:,out-script:,seed:,standalone -n 'bamsurgeon_split_BAM.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

#echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )

seed=$( date +"%Y" )
num_snvs=500
num_indels=100
num_svs=50
min_vaf=0.05
max_vaf=0.5
left_beta=2
right_beta=2
max_indel_len=18

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

        --num-snvs )
            case "$2" in
                "") shift 2 ;;
                *)  num_snvs=$2 ; shift 2 ;;
            esac ;;

        --num-indels )
            case "$2" in
                "") shift 2 ;;
                *)  num_indels=$2 ; shift 2 ;;
            esac ;;

        --num-svs )
            case "$2" in
                "") shift 2 ;;
                *)  num_svs=$2 ; shift 2 ;;
            esac ;;

        --max-indel-len )
            case "$2" in
                "") shift 2 ;;
                *)  max_indel_len=$2 ; shift 2 ;;
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

        --left-beta )
            case "$2" in
                "") shift 2 ;;
                *)  left_beta=$2 ; shift 2 ;;
            esac ;;

        --right-beta )
            case "$2" in
                "") shift 2 ;;
                *)  right_beta=$2 ; shift 2 ;;
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
    out_script="${logdir}/pickSites.${timestamp}.cmd"    
fi

if [[ $standalone ]]
then
    echo "#!/bin/bash" > $out_script
    echo "" >> $out_script
    echo "#$ -o ${logdir}" >> $out_script
    echo "#$ -e ${logdir}" >> $out_script
    echo "#$ -S /bin/bash" >> $out_script
    echo '#$ -l h_vmem=4G' >> $out_script
    echo 'set -e' >> $out_script
fi


echo "" >> $out_script

#1) Generate mutation sites and VAF's
echo "docker run -v /:/mnt -u $UID --rm lethalfang/bamsurgeon:1.1-3 bash -c \\" >> $out_script
echo "\"/usr/local/bamsurgeon/scripts/randomsites.py \\" >> $out_script
echo "--genome /mnt/${HUMAN_REFERENCE} \\" >> $out_script
echo "--seed $seed \\" >> $out_script
echo "--bed /mnt/${SELECTOR} \\" >> $out_script
echo "--numpicks ${num_snvs} \\" >> $out_script
echo "--minvaf $min_vaf \\" >> $out_script
echo "--maxvaf $max_vaf \\" >> $out_script
echo "--vafbeta1 $left_beta \\" >> $out_script
echo "--vafbeta2 $right_beta \\" >> $out_script
echo "--avoidN snv \\" >> $out_script
echo "| bedtools sort -header -faidx /mnt/${HUMAN_REFERENCE}.fai > /mnt/${outdir}/random_sSNV.bed\"" >> $out_script
echo "" >> $out_script


echo "docker run -v /:/mnt -u $UID --rm lethalfang/bamsurgeon:1.1-3 bash -c \\" >> $out_script
echo "\"/usr/local/bamsurgeon/scripts/randomsites.py \\" >> $out_script
echo "--genome /mnt/${HUMAN_REFERENCE} \\" >> $out_script
echo "--seed $seed \\" >> $out_script
echo "--bed /mnt/${SELECTOR} \\" >> $out_script
echo "--numpicks $num_indels \\" >> $out_script
echo "--minvaf $min_vaf \\" >> $out_script
echo "--maxvaf $max_vaf \\" >> $out_script
echo "--vafbeta1 $left_beta \\" >> $out_script
echo "--vafbeta2 $right_beta \\" >> $out_script
echo "--avoidN indel --maxlen ${max_indel_len}  \\" >> $out_script
echo "| vcfsorter.pl /mnt/${hg_dict} - > /mnt/${outdir}/random_sINDEL.bed\"" >> $out_script
echo "" >> $out_script


echo "docker run -v /:/mnt -u $UID --rm lethalfang/bamsurgeon:1.1-3 bash -c \\" >> $out_script
echo "\"/usr/local/bamsurgeon/scripts/randomsites.py \\" >> $out_script
echo "--genome /mnt/${HUMAN_REFERENCE} \\" >> $out_script
echo "--seed $seed \\" >> $out_script
echo "--bed /mnt/${SELECTOR} \\" >> $out_script
echo "--numpicks $num_svs \\" >> $out_script
echo "--minvaf $min_vaf \\" >> $out_script
echo "--maxvaf $max_vaf \\" >> $out_script
echo "--vafbeta1 $left_beta \\" >> $out_script
echo "--vafbeta2 $right_beta \\" >> $out_script
echo "sv --cnvfile /mnt/${outdir}/cnvfile.bed \\" >> $out_script
echo "| vcfsorter.pl /mnt/${hg_dict} - > /mnt/${outdir}/random_sSV.bed\"" >> $out_script
echo "" >> $out_script


echo "docker run -v /:/mnt -u $UID --rm lethalfang/bedtools:2.26.0 bash -c \\" >> $out_script
echo "\"bedtools sort -header -faidx /mnt/${HUMAN_REFERENCE}.fai -i /mnt/${outdir}/cnvfile.bed \\" >> $out_script
echo "> /mnt/${outdir}/sorted.cnvfile.bed\"" >> $out_script
echo "" >> $out_script


echo "docker run -v /:/mnt -u $UID --rm lethalfang/tabix:1.7 bgzip -f /mnt/${outdir}/sorted.cnvfile.bed" >> $out_script
echo "docker run -v /:/mnt -u $UID --rm lethalfang/tabix:1.7 tabix -f /mnt/${outdir}/sorted.cnvfile.bed.gz" >> $out_script
