#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long output-dir:,bam-out:,bam-in:,out-SM:,out-script:,standalone -n 'Reheader_SM.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

#echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )

keep_intermediates=0
outSM='TN_Merged'

while true; do
    case "$1" in
        -o | --output-dir )
            case "$2" in
                "") shift 2 ;;
                *)  outdir=$2 ; shift 2 ;;
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

        --out-SM )
            case "$2" in
                "") shift 2 ;;
                *)  outSM=$2 ; shift 2 ;;
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
    out_script="${logdir}/reheader.${timestamp}.cmd"    
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

# Uniform sample and read group names in the merged file
echo "singularity exec --bind /:/mnt   docker://lethalfang/bamsurgeon:1.1-3 \\" >> $out_script
echo "java -Xmx6g -jar /usr/local/bin/picard.jar AddOrReplaceReadGroups \\" >> $out_script
echo "I=/mnt/${outdir}/${inbam} \\" >> $out_script
echo "RGID=BAMSurgeon \\" >> $out_script
echo "RGLB=TNMerged \\" >> $out_script
echo "RGPL=illumina \\" >> $out_script
echo "RGPU=BAMSurgeon \\" >> $out_script
echo "RGSM=${outSM} \\" >> $out_script
echo "CREATE_INDEX=true \\" >> $out_script
echo "O=/mnt/${outdir}/${outbam}" >> $out_script
echo "" >> $out_script

# Remove temp files
echo "mv ${outdir}/${outbam%.bam}.bai ${outdir}/${outbam}.bai" >> $out_script
echo "" >> $out_script
