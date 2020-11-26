#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long output-dir:,bam-in:,bam-out:,selector:,out-script:,standalone -n 'split_BAM_by_BED.sh'  -- "$@"`

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
            
        --bam-in )
            case "$2" in
                "") shift 2 ;;
                *)  inbam=$2 ; shift 2 ;;
            esac ;;

        --bam-out )
            case "$2" in
                "") shift 2 ;;
                *)  outbam=$2 ; shift 2 ;;
            esac ;;

        --selector )
            case "$2" in
                "") shift 2 ;;
                *)  SELECTOR=$2 ; shift 2 ;;
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
    out_script="${logdir}/splitByBed.${timestamp}.cmd"    
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


echo "docker run --rm -v /:/mnt -u $UID --memory 4g lethalfang/samtools:1.7 bash -c \\" >> $out_script
echo "\"samtools view /mnt/${inbam} -L /mnt/${SELECTOR} -Sbh \\" >> $out_script
echo "> /mnt/${outdir}/${outbam}\"" >> $out_script

echo "" >> $out_script

echo "docker run --rm -v /:/mnt -u $UID --memory 4g lethalfang/samtools:1.7 \\" >> $out_script
echo "samtools index /mnt/${outdir}/${outbam}" >> $out_script
