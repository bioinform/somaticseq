#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long output-dir:,tumor-fq1:,tumor-fq2:,tumor-bam-header:,normal-fq1:,normal-fq2:,normal-bam-header:,genome-reference:,dbsnp:,known-indel:,tumor-in-bam:,tumor-out-bam:,normal-in-bam:,normal-out-bam:,MEM:,out-script:,action:,threads:,bwa,pre-realign-markdup,indel-realign,bqsr, -n 'fastq2bam_pipeline_singleThread.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )

tumor_bam_header='@RG\tID:myPipeline\tLB:myLibrary\tPL:illumina\tSM:TUMOR'
normal_bam_header='@RG\tID:myPipeline\tLB:myLibrary\tPL:illumina\tSM:NORMAL'
MEM=16
threads=1
action=echo
t_outbam=tumor
n_outbam=normal

while true; do
    case "$1" in
    
    -o | --output-dir )
        case "$2" in
            "") shift 2 ;;
            *)  outdir=$2 ; shift 2 ;;
        esac ;;

    --tumor-fq1 )
        case "$2" in
            "") shift 2 ;;
            *)  t_fq1=$2 ; shift 2 ;;
        esac ;;

    --tumor-fq2 )
        case "$2" in
            "") shift 2 ;;
            *)  t_fq2=$2 ; shift 2 ;;
        esac ;;

    --tumor-bam-header )
        case "$2" in
            "") shift 2 ;;
            *)  tumor_bam_header=$2 ; shift 2 ;;
        esac ;;

    --normal-fq1 )
        case "$2" in
            "") shift 2 ;;
            *)  n_fq1=$2 ; shift 2 ;;
        esac ;;

    --normal-fq2 )
        case "$2" in
            "") shift 2 ;;
            *)  n_fq2=$2 ; shift 2 ;;
        esac ;;

    --normal-bam-header )
        case "$2" in
            "") shift 2 ;;
            *)  normal_bam_header=$2 ; shift 2 ;;
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

    --bwa )
        bwa=1 ; shift ;;

    --pre-realign-markdup )
        pre_realign_markdup=1 ; shift ;;

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


if [[ ${out_script_name} ]]
then
    out_script="${out_script_name}"
else
    out_script="${logdir}/fastq2bam.${timestamp}.cmd"    
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

files_to_delete=''

echo 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $out_script
echo "" >> $out_script


if [[ $bwa ]]
then

    if [[ ${t_fq1} && ${t_fq2} ]]
    then
        $MYDIR/bwa_mem_pe.sh \
        --fq1              ${t_fq1} \
        --fq2              ${t_fq2} \
        --genome-reference ${GENOME_REFERENCE} \
        --output-dir       ${outdir} \
        --out-bam          ${t_outbam}.sorted.bam \
        --bam-header       "${tumor_bam_header}" \
        --threads          ${threads} \
        --out-script       ${out_script}

        latest_tumor_bam="${outdir}/${t_outbam}.sorted.bam"
    fi

    if [[ ${n_fq1} && ${n_fq2} ]]
    then
        $MYDIR/bwa_mem_pe.sh \
        --fq1              ${n_fq1} \
        --fq2              ${n_fq2} \
        --genome-reference ${GENOME_REFERENCE} \
        --output-dir       ${outdir} \
        --out-bam          ${n_outbam}.sorted.bam \
        --bam-header       "${normal_bam_header}" \
        --threads          ${threads} \
        --out-script       ${out_script}
        
        latest_normal_bam="${outdir}/${n_outbam}.sorted.bam"
    fi

fi


if [[ $pre_realign_markdup ]]
then

    if [[ ${latest_tumor_bam} ]]
    then
        $MYDIR/markdup.sh \
        --output-dir ${outdir} \
        --in-bam     ${latest_tumor_bam} \
        --out-bam    ${t_outbam}.markdup.bam \
        --out-script ${out_script}
    
        latest_tumor_bam="${outdir}/${t_outbam}.markdup.bam"
    fi
    
  
    if [[ ${latest_normal_bam} ]]
    then
    
        $MYDIR/markdup.sh \
        --output-dir ${outdir} \
        --in-bam     ${latest_normal_bam} \
        --out-bam    ${n_outbam}.markdup.bam \
        --out-script ${out_script}
    
        latest_normal_bam="${outdir}/${n_outbam}.markdup.bam"
    fi
fi



if [[ $indel_realign ]]
then

    if [[ ${latest_tumor_bam} && ${latest_normal_bam} ]]
    then
    
        $MYDIR/jointIndelRealign.sh \
        --output-dir       ${outdir} \
        --normal-bam       ${latest_normal_bam} \
        --tumor-bam        ${latest_tumor_bam} \
        --genome-reference ${GENOME_REFERENCE} \
        --threads          ${threads} \        
        --out-script       ${out_script}
        
        latest_normal_bam=${latest_normal_bam%.bam}.jointRealigned.bam
        latest_tumor_bam=${latest_tumor_bam%.bam}.jointRealigned.bam

    elif [[ ${latest_tumor_bam} ]]
    then
        
        $MYDIR/singleIndelRealign.sh \
        --output-dir       ${outdir} \
        --tumor-bam        ${latest_tumor_bam} \
        --genome-reference ${GENOME_REFERENCE} \
        --threads          ${threads} \
        --out-script       ${out_script}
        
        latest_tumor_bam=${latest_tumor_bam%.bam}.indelRealigned.bam
    
    elif [[ ${latest_normal_bam} ]]
    then
    
        $MYDIR/singleIndelRealign.sh \
        --output-dir       ${outdir} \
        --tumor-bam        ${latest_normal_bam} \
        --genome-reference ${GENOME_REFERENCE} \
        --threads          ${threads} \        
        --out-script       ${out_script}
        
        latest_normal_bam=${latest_normal_bam%.bam}.indelRealigned.bam
    fi
    
fi


if [[ $bqsr ]]
then

    if [[ ${latest_tumor_bam} ]]
    then
        $MYDIR/BQSR.sh \
        --output-dir       ${outdir} \
        --in-bam           ${latest_tumor_bam} \
        --out-bam          ${t_outbam}.bam \
        --genome-reference ${GENOME_REFERENCE} \
        --dbsnp            ${dbsnp} \
        --out-script       ${out_script}
    fi
    
    if [[ ${latest_normal_bam} ]]
    then
        $MYDIR/BQSR.sh \
        --output-dir       ${outdir} \
        --in-bam           ${latest_normal_bam} \
        --out-bam          ${n_outbam}.bam \
        --genome-reference ${GENOME_REFERENCE} \
        --dbsnp            ${dbsnp} \
        --out-script       ${out_script}
    fi

fi
