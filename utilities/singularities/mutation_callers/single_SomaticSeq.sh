#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long out-dir:,in-bam:,human-reference:,selector:,exclude:,dbsnp:,cosmic:,MEM:,action:,mutect2:,varscan:,vardict:,lofreq:,scalpel:,strelka:,ada-r-script:,classifier-snv:,classifier-indel:,truth-snv:,truth-indel:,extra-arguments: -n 'submit_SomaticSeq.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

#echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )
action=echo
MEM=6

while true; do
    case "$1" in

    --out-dir )
        case "$2" in
            "") shift 2 ;;
            *)  outdir=$2 ; shift 2 ;;
        esac ;;

    --in-bam )
        case "$2" in
            "") shift 2 ;;
            *)  tumor_bam=$2 ; shift 2 ;;
        esac ;;

    --human-reference )
        case "$2" in
            "") shift 2 ;;
            *)  HUMAN_REFERENCE=$2 ; shift 2 ;;
        esac ;;

    --selector )
        case "$2" in
            "") shift 2 ;;
            *)  SELECTOR=$2 ; shift 2 ;;
        esac ;;

    --exclude )
        case "$2" in
            "") shift 2 ;;
            *)  EXCLUSION=$2 ; shift 2 ;;
        esac ;;

    --dbsnp )
        case "$2" in
            "") shift 2 ;;
            *)  dbsnp=$2 ; shift 2 ;;
        esac ;;

    --cosmic )
        case "$2" in
            "") shift 2 ;;
            *)  cosmic=$2 ; shift 2 ;;
        esac ;;

    --MEM )
        case "$2" in
            "") shift 2 ;;
            *) MEM=$2 ; shift 2 ;;
        esac ;;

    --action )
        case "$2" in
            "") shift 2 ;;
            *)  action=$2 ; shift 2 ;;
        esac ;;

    --mutect2 )
        case "$2" in
            "") shift 2 ;;
            *)  mutect2_vcf=$2 ; shift 2 ;;
        esac ;;

    --varscan )
        case "$2" in
            "") shift 2 ;;
            *)  varscan_vcf=$2 ; shift 2 ;;
        esac ;;

    --vardict )
        case "$2" in
            "") shift 2 ;;
            *)  vardict_vcf=$2 ; shift 2 ;;
        esac ;;

    --lofreq )
        case "$2" in
            "") shift 2 ;;
            *)  lofreq_vcf=$2 ; shift 2 ;;
        esac ;;


    --scalpel )
        case "$2" in
            "") shift 2 ;;
            *)  scalpel_vcf=$2 ; shift 2 ;;
        esac ;;

    --strelka )
        case "$2" in
            "") shift 2 ;;
            *)  strelka_vcf=$2 ; shift 2 ;
        esac ;;

    --ada-r-script )
        case "$2" in
            "") shift 2 ;;
            *)  ada_r_script=$2 ; shift 2 ;;
        esac ;;
        
    --classifier-snv )
        case "$2" in
            "") shift 2 ;;
            *)  classifier_snv=$2 ; shift 2 ;;
        esac ;;
        
    --classifier-indel )
        case "$2" in
            "") shift 2 ;;
            *)  classifier_indel=$2 ; shift 2 ;;
        esac ;;
        
    --truth-snv )
        case "$2" in
            "") shift 2 ;;
            *)  truth_snv=$2 ; shift 2 ;;
        esac ;;
        
    --truth-indel )
        case "$2" in
            "") shift 2 ;;
            *)  truth_indel=$2 ; shift 2 ;;
        esac ;;

    --extra-arguments )
        case "$2" in
            "") shift 2 ;;
            *) extra_arguments=$2 ; shift 2 ;;
        esac ;;

    -- ) shift; break ;;
    * ) break ;;
    esac

done

# VERSION='2.4.0'
VERSION=`cat ${MYDIR}/../../../VERSION | sed 's/##SomaticSeq=v//'`

logdir=${outdir}/logs
mkdir -p ${logdir}

out_script=${outdir}/logs/sseq_${timestamp}.cmd

selector_text=''
if [[ -r $SELECTOR ]]; then
    selector_text="--inclusion-region /mnt/${SELECTOR}"
fi

exclusion_text=''
if [[ -r $EXCLUSION ]]; then
    exclusion_text="--exclusion-region /mnt/${EXCLUSION}"
fi


dbsnp_text=''
if [[ -r $dbsnp ]]; then
    dbsnp_text="--dbsnp /mnt/${dbsnp}"
fi

cosmic_text=''
if [[ -r $cosmic ]]; then
    cosmic_text="--cosmic /mnt/${cosmic}"
fi

# VCF inputs
if [[ $mutect2_vcf ]];  then mutect2_text="--mutect2 /mnt/${mutect2_vcf}";  fi
if [[ $varscan_vcf ]];  then varscan_text="--varscan /mnt/${varscan_vcf}";  fi
if [[ $vardict_vcf ]];  then vardict_text="--vardict /mnt/${vardict_vcf}";  fi
if [[ $lofreq_vcf ]];   then lofreq_text="--lofreq /mnt/${lofreq_vcf}";     fi
if [[ $scalpel_vcf ]];  then scalpel_text="--scalpel /mnt/${scalpel_vcf}";  fi
if [[ $strelka_vcf ]];  then strelka_text="--strelka /mnt/${strelka_vcf}";  fi

# SomaticSeq modes:
if [[ $classifier_snv ]];   then classifier_snv_text="--classifier-snv ${classifier_snv}";       fi
if [[ $classifier_indel ]]; then classifier_indel_text="--classifier-indel ${classifier_indel}"; fi
if [[ $truth_snv ]];        then truth_snv_text="--truth-snv ${truth_snv}"                     ; fi
if [[ $truth_indel ]];      then truth_indel_text="--truth-indel ${truth_indel}"               ; fi
if [[ $ada_r_script ]];     then ada_r_script_text="--ada-r-script ${ada_r_script}"            ; fi

echo "#!/bin/bash" > $out_script
echo "" >> $out_script

echo "#$ -o ${logdir}" >> $out_script
echo "#$ -e ${logdir}" >> $out_script
echo "#$ -S /bin/bash" >> $out_script
echo "#$ -l h_vmem=${MEM}G" >> $out_script # ML may require lots of RAM
echo 'set -e' >> $out_script
echo "" >> $out_script

echo 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $out_script
echo "" >> $out_script

echo "singularity exec --bind /:/mnt docker://lethalfang/somaticseq:${VERSION} \\" >> $out_script
echo "/opt/somaticseq/ssSomaticSeq.Wrapper.sh \\" >> $out_script
echo "--output-dir /mnt/${outdir} \\" >> $out_script
echo "--genome-reference /mnt/${HUMAN_REFERENCE} \\" >> $out_script
echo "--in-bam /mnt/${tumor_bam} \\" >> $out_script
echo "$mutect2_text \\" >> $out_script
echo "$varscan_text \\" >> $out_script
echo "$vardict_text \\" >> $out_script
echo "$lofreq_text \\" >> $out_script
echo "$scalpel_text \\" >> $out_script
echo "$strelka_text \\" >> $out_script
echo "$selector_text \\" >> $out_script
echo "$exclusion_text \\" >> $out_script
echo "$cosmic_text \\" >> $out_script
echo "$dbsnp_text \\" >> $out_script
echo "$classifier_snv_text \\" >> $out_script
echo "$classifier_indel_text \\" >> $out_script
echo "$truth_snv_text \\" >> $out_script
echo "$truth_indel_text \\" >> $out_script
echo "$ada_r_script_text \\" >> $out_script
echo "${extra_arguments} \\" >> $out_script
echo "--gatk /opt/GATK/GenomeAnalysisTK.jar" >> $out_script

echo "" >> $out_script
echo 'echo -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $out_script

$action $out_script
