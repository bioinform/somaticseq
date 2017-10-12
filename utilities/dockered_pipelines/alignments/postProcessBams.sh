#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long output-dir:,genome-reference:,selector:,dbsnp:,known-indel:,tumor-in-bam:,tumor-out-bam:,normal-in-bam:,normal-out-bam:,MEM:,out-script:,action:,threads:,markdup:,indel-realign,bqsr, -n 'postProcessBams.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )

MEM=14
threads=1
action=echo

while true; do
    case "$1" in
    
    -o | --output-dir )
        case "$2" in
            "") shift 2 ;;
            *)  outdir=$2 ; shift 2 ;;
        esac ;;

    --tumor-in-bam )
        case "$2" in
            "") shift 2 ;;
            *)  t_inbam=$2 ; shift 2 ;;
        esac ;;

    --tumor-out-bam )
        case "$2" in
            "") shift 2 ;;
            *)  t_outbam=$2 ; shift 2 ;;
        esac ;;

    --normal-in-bam )
        case "$2" in
            "") shift 2 ;;
            *)  n_inbam=$2 ; shift 2 ;;
        esac ;;

    --normal-out-bam )
        case "$2" in
            "") shift 2 ;;
            *)  n_outbam=$2 ; shift 2 ;;
        esac ;;

    --genome-reference )
        case "$2" in
            "") shift 2 ;;
            *)  GENOME_REFERENCE=$2 ; shift 2 ;;
        esac ;;

    --selector )
        case "$2" in
            "") shift 2 ;;
            *)  SELECTOR=$2 ; shift 2 ;;
        esac ;;

    --dbsnp )
        case "$2" in
            "") shift 2 ;;
            *)  dbsnp=$2 ; shift 2 ;;
        esac ;;

    --known-indel )
        case "$2" in
            "") shift 2 ;;
            *)  known_indel=$2 ; shift 2 ;;
        esac ;;

    --threads )
        case "$2" in
            "") shift 2 ;;
            *)  threads=$2 ; shift 2 ;;
        esac ;;

    --MEM )
        case "$2" in
            "") shift 2 ;;
            *)  MEM=$2 ; shift 2 ;;
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

    --indel-realign )
        indel_realign=1 ; shift ;;

    --bqsr )
        bqsr=1 ; shift ;;

    -- ) shift; break ;;
    * ) break ;;
    
    esac
done

VERSION='latest'

hg_dict=${GENOME_REFERENCE%\.fa*}.dict

logdir=${outdir}/logs
mkdir -p ${logdir}


if [[ $threads -ge 1 ]]
then

    if [[ $SELECTOR ]]
    then
        cp $SELECTOR ${outdir}/genome.bed
    else
        cat ${GENOME_REFERENCE}.fai | awk -F "\t" '{print $1 "\t0\t" $2}' > ${outdir}/genome.bed
    fi
        
    docker run --rm -v /:/mnt -u $UID lethalfang/somaticseq:${VERSION} \
    /opt/somaticseq/utilities/split_Bed_into_equal_regions.py \
    -infile /mnt/${outdir}/genome.bed -num $threads -outfiles /mnt/${outdir}/bed
fi




ith_thread=1
while [[ $ith_thread -le $threads ]]
do

    mkdir -p ${outdir}/${ith_thread}/logs
    mv ${outdir}/${ith_thread}.bed ${outdir}/${ith_thread}


    if [[ ${out_script_name} ]]
    then
        out_script="${outdir}/${ith_thread}/logs/${out_script_name}"
    else
        out_script="${outdir}/${ith_thread}/logs/postProcessBams.${timestamp}.cmd"    
    fi
    
    
    echo "#!/bin/bash" > $out_script
    echo "" >> $out_script
    
    echo "#$ -o ${logdir}" >> $out_script
    echo "#$ -e ${logdir}" >> $out_script
    echo "#$ -S /bin/bash" >> $out_script
    echo "#$ -l h_vmem=${MEM}G" >> $out_script
    #echo "#$ -pe smp 1" >> $out_script
    
    echo 'set -e' >> $out_script
    echo "" >> $out_script
    
    files_to_delete=''
    
    echo 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $out_script
    echo "" >> $out_script

    if [[ $indel_realign ]]
    then
    
        if [[ ${t_inbam} && ${n_inbam} ]]
        then
        
            $MYDIR/jointIndelRealign.sh \
            --output-dir       ${outdir}/${ith_thread} \
            --normal-bam       ${n_inbam} \
            --tumor-bam        ${t_inbam} \
            --selector         ${outdir}/${ith_thread}/${ith_thread}.bed \
            --genome-reference ${GENOME_REFERENCE} \
            --out-script       ${out_script}
            
            latest_normal_bam=${n_inbam%.bam}.jointRealigned.bam
            latest_tumor_bam=${t_inbam%.bam}.jointRealigned.bam
    
        elif [[ ${t_inbam} ]]
        then
            
            $MYDIR/singleIndelRealign.sh \
            --output-dir       ${outdir}/${ith_thread} \
            --tumor-bam        ${t_inbam} \
            --selector         ${outdir}/${ith_thread}/${ith_thread}.bed \            
            --genome-reference ${GENOME_REFERENCE} \
            --out-script       ${out_script}
            
            latest_tumor_bam=${t_inbam%.bam}.indelRealigned.bam
        
        elif [[ ${n_inbam} ]]
        then
        
            $MYDIR/singleIndelRealign.sh \
            --output-dir       ${outdir}/${ith_thread} \
            --tumor-bam        ${n_inbam} \
            --selector         ${outdir}/${ith_thread}/${ith_thread}.bed \            
            --genome-reference ${GENOME_REFERENCE} \
            --out-script       ${out_script}
            
            latest_normal_bam=${n_inbam%.bam}.indelRealigned.bam
        fi
        
    fi
    
    
    if [[ $bqsr ]]
    then
    
        if [[ ${latest_tumor_bam} ]]
        then
            $MYDIR/BQSR.sh \
            --output-dir       ${outdir}/${ith_thread} \
            --in-bam           ${latest_tumor_bam} \
            --out-bam          ${t_outbam} \
            --genome-reference ${GENOME_REFERENCE} \
            --dbsnp            ${dbsnp} \
            --out-script       ${out_script}
        fi
        
        if [[ ${latest_normal_bam} ]]
        then
            $MYDIR/BQSR.sh \
            --output-dir       ${outdir}/${ith_thread} \
            --in-bam           ${latest_normal_bam} \
            --out-bam          ${n_outbam} \
            --genome-reference ${GENOME_REFERENCE} \
            --dbsnp            ${dbsnp} \
            --out-script       ${out_script}
        fi
    
    fi

        
    ith_thread=$(( $ith_thread + 1))

done
