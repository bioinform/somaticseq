#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long out-dir:,bam:,genome-reference:,extra-arguments:,action:,MEM:,threads:,out-script:,standalone -n 'manta.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

#echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )
action=echo
MEM=8
threads=12

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

    --genome-reference )
        case "$2" in
            "") shift 2 ;;
            *)  GENOME_REFERENCE=$2 ; shift 2 ;;
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

    --threads )
        case "$2" in
            "") shift 2 ;;
            *)  threads=$2 ; shift 2 ;;
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
    echo '#$ -l h_vmem=4G' >> $out_script
    echo "#$ -pe smp ${threads}" >> $out_script
    echo 'set -e' >> $out_script
fi

echo "" >> $out_script
echo 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $out_script
echo "" >> $out_script


echo "docker run -v /:/mnt -u $UID --rm lethalfang/manta:1.4.0 \\" >> $out_script
echo "/opt/manta/bin/configManta.py \\" >> $out_script
echo "--bam            /mnt/${bamFile} \\" >> $out_script
echo "--referenceFasta /mnt/${GENOME_REFERENCE} \\" >> $out_script
echo "--runDir         /mnt/${outdir}" >> $out_script

echo "" >> $out_script

echo "docker run -v /:/mnt -u $UID --rm lethalfang/manta:1.4.0 /mnt/${outdir}/runWorkflow.py -m local -j $thread" >> $out_script
