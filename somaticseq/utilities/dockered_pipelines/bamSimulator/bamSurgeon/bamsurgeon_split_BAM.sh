#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long output-dir:,genome-reference:,bam-out1:,bam-out2:,bam-in:,split-proportion:,down-sample:,seed:,out-script:,clean-bam,standalone -n 'bamsurgeon_split_BAM.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

#echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )
seed=$( date +"%Y" )
proportion=0.5
down_sample=1

while true; do
    case "$1" in
        -o | --output-dir )
            case "$2" in
                "") shift 2 ;;
                *)  outdir=$2 ; shift 2 ;;
            esac ;;

        --bam-in )
            case "$2" in
                "") shift 2 ;;
                *)  inbam=$2 ; shift 2 ;;
            esac ;;

        --bam-out1 )
            case "$2" in
                "") shift 2 ;;
                *)  outbam1=$2 ; shift 2 ;;
            esac ;;

        --bam-out2 )
            case "$2" in
                "") shift 2 ;;
                *)  outbam2=$2 ; shift 2 ;;
            esac ;;

        --genome-reference )
            case "$2" in
                "") shift 2 ;;
                *)  HUMAN_REFERENCE=$2 ; shift 2 ;;
            esac ;;

        --split-proportion )
            case "$2" in
                "") shift 2 ;;
                *)  proportion=$2 ; shift 2 ;;
            esac ;;

        --down-sample )
            case "$2" in
                "") shift 2 ;;
                *)  down_sample=$2 ; shift 2 ;;
            esac ;;

        --seed )
            case "$2" in
                "") shift 2 ;;
                *)  seed=$2 ; shift 2 ;;
            esac ;;

        --out-script )
            case "$2" in
                "") shift 2 ;;
                *)  out_script_name=$2 ; shift 2 ;;
            esac ;;

        --clean-bam )
            clean_bam=1 ; shift ;;

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
    out_script="${logdir}/splitBams.${timestamp}.cmd"    
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


# Then you can split
echo "docker run -v /:/mnt -u $UID --rm --memory 8g lethalfang/bamsurgeon:1.1-3 \\" >> $out_script
echo "/usr/local/bamsurgeon/scripts/sortedBamSplit.py \\" >> $out_script
echo "--bam /mnt/${inbam} \\" >> $out_script
echo "--proportion ${proportion} \\" >> $out_script
echo "--downsample ${down_sample} \\" >> $out_script
echo "--pick1 /mnt/${outdir}/${outbam1} \\" >> $out_script
echo "--pick2 /mnt/${outdir}/${outbam2} \\" >> $out_script
echo "--seed ${seed}" >> $out_script
echo "" >> $out_script

echo "docker run -v /:/mnt -u $UID --rm --memory 8g lethalfang/samtools:1.7 samtools index /mnt/${outdir}/${outbam1}" >> $out_script
echo "docker run -v /:/mnt -u $UID --rm --memory 8g lethalfang/samtools:1.7 samtools index /mnt/${outdir}/${outbam2}" >> $out_script
echo "" >> $out_script
