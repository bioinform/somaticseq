#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long output-dir:,in-bam:,out-bam:,out-script:,standalone, -n 'markdup.sh'  -- "$@"`

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
    out_script="${logdir}/markdup.${timestamp}.cmd"
fi

if [ ! $outBam ]
then
    outBam=${inBam%.bam}.markdup.bam
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

echo "mkdir -p ${outdir}/${timestamp}.temp" >> $out_script
echo "" >> $out_script

echo "docker run --rm -v /:/mnt -u $UID lethalfang/picard:2.18.4 \\" >> $out_script
echo "java -Xmx16g -jar /opt/picard.jar MarkDuplicates \\" >> $out_script
echo "I=/mnt/${inBam} \\" >> $out_script
echo "M=/mnt/${outdir}/${outBam%.bam} \\" >> $out_script
echo "CREATE_INDEX=true \\" >> $out_script
echo "ASSUME_SORTED=true \\" >> $out_script
echo "TMP_DIR=/mnt/${outdir}/${timestamp}.temp \\" >> $out_script
#echo "MINIMUM_DISTANCE=1000 \\" >> $out_script
echo "O=/mnt/${outdir}/${outBam}" >> $out_script

echo "" >> $out_script

echo "mv ${outdir}/${outBam%.bam}.bai ${outdir}/${outBam}.bai" >> $out_script
echo "rm -r ${outdir}/${timestamp}.temp" >> $out_script
