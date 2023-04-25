#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long out-dir:,bam:,in-vcf:,sample-name:,canvas-reference:,extra-arguments:,action:,MEM:,threads:,out-script:,standalone -n 'canvas.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

#echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )
action=echo
MEM=8
threads=12
SAMPLE_NAME='Canvas'

while true; do
    case "$1" in
    
    -o | --out-dir )
        case "$2" in
            "") shift 2 ;;
            *)  outdir=$2 ; shift 2 ;;
        esac ;;
    
    --bam )
        case "$2" in
            "") shift 2 ;;
            *)  bamFile=$2 ; shift 2 ;;
        esac ;;

    --in-vcf )
        case "$2" in
            "") shift 2 ;;
            *)  inVcf=$2 ; shift 2 ;;
        esac ;;

    --sample-name )
        case "$2" in
            "") shift 2 ;;
            *)  SAMPLE_NAME=$2 ; shift 2 ;;
        esac ;;

    --canvas-reference )
        case "$2" in
            "") shift 2 ;;
            *)  CANVAS_REFERENCE=$2 ; shift 2 ;;
        esac ;;

    --genome-reference-dir )
        case "$2" in
            "") shift 2 ;;
            *)  GENOMIC_REFERENCE_DIR=$2 ; shift 2 ;;
        esac ;;

    --filter-bed )
        case "$2" in
            "") shift 2 ;;
            *)  filterBed=$2 ; shift 2 ;;
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
    out_script="${logdir}/canvas.${timestamp}.cmd"    
fi


if [[ $standalone ]]
then
    echo "#!/bin/bash" > $out_script
    echo "" >> $out_script
    echo "#$ -o ${logdir}" >> $out_script
    echo "#$ -e ${logdir}" >> $out_script
    echo "#$ -S /bin/bash" >> $out_script
    echo '#$ -l h_vmem=12G' >> $out_script
    echo "#$ -pe smp ${threads}" >> $out_script
    echo 'set -e' >> $out_script
fi

echo "" >> $out_script
echo 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $out_script
echo "" >> $out_script


# Cannot yet control number of threads to invoke
# First, create Canvas-ready fasta
# docker run --rm -v /:/mnt -u $UID lethalfang/canvas:1.35.1 bash -c 'export COMPlus_gcAllowVeryLargeObjects=1 && dotnet /opt/Canvas/Tools/FlagUniqueKmers/FlagUniqueKmers.dll /mnt/sc1/groups/bfx-red/data/datainsights/SEQC2_Resources/GRCh38.d1.vd1.fa /mnt/sc1/groups/bfx-red/data/datainsights/SEQC2_Resources/GRCh38.d1.vd1.Canvas-ready.fasta'

# Run Canvas
echo "docker run -u $UID --rm -v /:/mnt lethalfang/canvas:1.35.1 \\" >> $out_script
echo "dotnet /opt/Canvas/Canvas.dll Germline-WGS \\" >> $out_script
echo "--bam                 /mnt/${bamFile} \\" >> $out_script
echo "--sample-b-allele-vcf /mnt/${inVcf} \\" >> $out_script
echo "--sample-name         ${SAMPLE_NAME} \\" >> $out_script
echo "--reference           /mnt/${CANVAS_REFERENCE} \\" >> $out_script
echo "--genome-folder       /mnt/${GENOMIC_REFERENCE_DIR} \\" >> $out_script
echo "--filter-bed          /mnt/${filterBed} \\" >> $out_script
echo "--output              /mnt/${out_dir} " >> $out_script

echo "" >> $out_script
echo 'echo -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $out_script
