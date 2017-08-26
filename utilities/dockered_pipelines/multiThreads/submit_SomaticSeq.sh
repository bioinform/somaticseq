#!/bin/bash
# Use getopt instead of getopts for long options

set -e

OPTS=`getopt -o o: --long out-dir:,tumor-bam:,normal-bam:,human-reference:,selector:,dbsnp:,cosmic:,action:,mutect:,indelocator:,mutect2:,varscan-snv:,varscan-indel:,jsm:,sniper:,vardict:,muse:,lofreq-snv:,lofreq-indel:,scalpel:,strelka-snv:,strelka-indel:,ada-r-script:,classifier-snv:,classifier-indel:,truth-snv:,truth-indel: -n 'submit_SomaticSeq.sh'  -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

echo "$OPTS"
eval set -- "$OPTS"

MYDIR="$( cd "$( dirname "$0" )" && pwd )"

timestamp=$( date +"%Y-%m-%d_%H-%M-%S_%N" )
action=echo

while true; do
    case "$1" in

    --out-dir )
        case "$2" in
            "") shift 2 ;;
            *)  outdir=$2 ; shift 2 ;;
        esac ;;

    --tumor-bam )
        case "$2" in
            "") shift 2 ;;
            *)  tumor_bam=$2 ; shift 2 ;;
        esac ;;

    --normal-bam )
        case "$2" in
            "") shift 2 ;;
            *)  normal_bam=$2 ; shift 2 ;;
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

    --action )
        case "$2" in
            "") shift 2 ;;
            *)  action=$2 ; shift 2 ;;
        esac ;;

    --mutect )
        case "$2" in
            "") shift 2 ;;
            *)  mutect_vcf=$2 ; shift 2 ;;
        esac ;;

    --indelocator )
        case "$2" in
            "") shift 2 ;;
            *)  indelocator_vcf=$2 ; shift 2 ;;
        esac ;;

    --mutect2 )
        case "$2" in
            "") shift 2 ;;
            *)  mutect2_vcf=$2 ; shift 2 ;;
        esac ;;

    --varscan-snv )
        case "$2" in
            "") shift 2 ;;
            *)  varscan_snv_vcf=$2 ; shift 2 ;;
        esac ;;

    --varscan-indel )
        case "$2" in
            "") shift 2 ;;
            *)  varscan_indel_vcf=$2 ; shift 2 ;;
        esac ;;

    --jsm )
        case "$2" in
            "") shift 2 ;;
            *)  jsm_vcf=$2 ; shift 2 ;;
        esac ;;

    --sniper )
        case "$2" in
            "") shift 2 ;;
            *)  sniper_vcf=$2 ; shift 2 ;;
        esac ;;

    --vardict )
        case "$2" in
            "") shift 2 ;;
            *)  vardict_vcf=$2 ; shift 2 ;;
        esac ;;

    --muse )
        case "$2" in
            "") shift 2 ;;
             *)  muse_vcf=$2 ; shift 2 ;;
        esac ;;

    --lofreq-snv )
        case "$2" in
            "") shift 2 ;;
            *)  lofreq_snv_vcf=$2 ; shift 2 ;;
        esac ;;

    --lofreq-indel )
        case "$2" in
            "") shift 2 ;;
            *)  lofreq_indel_vcf=$2 ; shift 2 ;;
        esac ;;

    --scalpel )
        case "$2" in
            "") shift 2 ;;
            *)  scalpel_vcf=$2 ; shift 2 ;;
        esac ;;

    --strelka-snv )
        case "$2" in
            "") shift 2 ;;
            *)  strelka_snv_vcf=$2 ; shift 2 ;
        esac ;;

    --strelka-indel )
        case "$2" in
            "") shift 2 ;;
            *)  strelka_indel_vcf=$2 ; shift 2 ;;
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

    -- ) shift; break ;;
    * ) break ;;
    esac

done

VERSION='2.3.2'

logdir=${outdir}/logs
mkdir -p ${logdir}

sseq_script=${outdir}/logs/sseq_${timestamp}.cmd

selector_text=''
if [[ -r $SELECTOR ]]; then
    selector_text="--inclusion-region /mnt/${SELECTOR}"
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
if [[ $mutect_vcf ]];        then mutect_text="--mutect /mnt/${mutect_vcf}";                      fi
if [[ $indelocator_vcf ]];   then indelocator_text="--indelocator /mnt/${indelocator_vcf}";       fi
if [[ $mutect2_vcf ]];       then mutect2_text="--mutect2 /mnt/${mutect2_vcf}";                   fi
if [[ $varscan_snv_vcf ]];   then varscan_snv_text="--varscan-snv /mnt/${varscan_snv_vcf}";       fi
if [[ $varscan_indel_vcf ]]; then varscan_indel_text="--varscan-indel /mnt/${varscan_indel_vcf}"; fi
if [[ $jsm_vcf ]];           then jsm_text="--jsm /mnt/${jsm_vcf}";                               fi
if [[ $sniper_vcf ]];        then sniper_text="--sniper /mnt/${sniper_vcf}";                      fi
if [[ $vardict_vcf ]];       then vardict_text="--vardict /mnt/${vardict_vcf}";                   fi
if [[ $muse_vcf ]];          then muse_text="--muse /mnt/${muse_vcf}";                            fi
if [[ $lofreq_snv_vcf ]];    then lofreq_snv_text="--lofreq-snv /mnt/${lofreq_snv_vcf}";          fi
if [[ $lofreq_indel_vcf ]];  then lofreq_indel_text="--lofreq-indel /mnt/${lofreq_indel_vcf}";    fi
if [[ $scalpel_vcf ]];       then scalpel_text="--scalpel /mnt/${scalpel_vcf}";                   fi
if [[ $strelka_snv_vcf ]];   then strelka_snv_text="--strelka-snv /mnt/${strelka_snv_vcf}";       fi
if [[ $strelka_indel_vcf ]]; then strelka_indel_text="--strelka-indel /mnt/${strelka_indel_vcf}"; fi

# SomaticSeq modes:
if [[ $classifier_snv ]];   then classifier_snv_text="--classifier-snv ${classifier_snv}";       fi
if [[ $classifier_indel ]]; then classifier_indel_text="--classifier-indel ${classifier_indel}"; fi
if [[ $truth_snv ]];        then truth_snv_text="--truth-snv ${truth_snv}"                     ; fi
if [[ $truth_indel ]];      then truth_indel_text="--truth-indel ${truth_indel}"               ; fi
if [[ $ada_r_script ]];     then ada_r_script_text="--ada-r-script ${ada_r_script}"            ; fi

echo "#!/bin/bash" > $sseq_script
echo "" >> $sseq_script

echo "#$ -o ${logdir}" >> $sseq_script
echo "#$ -e ${logdir}" >> $sseq_script
echo "#$ -S /bin/bash" >> $sseq_script
#echo '#$ -l h_vmem=6G' >> $sseq_script # ML may require lots of RAM
echo 'set -e' >> $sseq_script
echo "" >> $sseq_script

echo 'echo -e "Start at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $sseq_script
echo "" >> $sseq_script

echo "docker pull lethalfang/somaticseq:${VERSION}" >> $sseq_script
echo "" >> $sseq_script

echo "docker run --rm -v /:/mnt -u $UID -i lethalfang/somaticseq:${VERSION} \\" >> $sseq_script
echo "/opt/somaticseq/SomaticSeq.Wrapper.sh \\" >> $sseq_script
echo "--output-dir       /mnt/${outdir} \\" >> $sseq_script
echo "--genome-reference /mnt/${HUMAN_REFERENCE} \\" >> $sseq_script
echo "--tumor-bam        /mnt/${tumor_bam} \\" >> $sseq_script
echo "--normal-bam       /mnt/${normal_bam} \\" >> $sseq_script
echo "$mutect_text \\" >> $sseq_script
echo "$indelocator_text \\" >> $sseq_script
echo "$mutect2_text \\" >> $sseq_script
echo "$varscan_snv_text \\" >> $sseq_script
echo "$varscan_indel_text \\" >> $sseq_script
echo "$jsm_text \\" >> $sseq_script
echo "$sniper_text \\" >> $sseq_script
echo "$vardict_text \\" >> $sseq_script
echo "$muse_text \\" >> $sseq_script
echo "$lofreq_snv_text \\" >> $sseq_script
echo "$lofreq_indel_text \\" >> $sseq_script
echo "$scalpel_text \\" >> $sseq_script
echo "$strelka_snv_text \\" >> $sseq_script
echo "$strelka_indel_text \\" >> $sseq_script
echo "$selector_text \\" >> $sseq_script
echo "$cosmic_text \\" >> $sseq_script
echo "$dbsnp_text \\" >> $sseq_script
echo "$classifier_snv_text \\" >> $sseq_script
echo "$classifier_indel_text \\" >> $sseq_script
echo "$truth_snv_text \\" >> $sseq_script
echo "$truth_indel_text \\" >> $sseq_script
echo "$ada_r_script_text \\" >> $sseq_script
echo "--gatk /opt/GATK/GenomeAnalysisTK.jar" >> $sseq_script

echo "" >> $sseq_script
echo 'echo -e "Done at `date +"%Y/%m/%d %H:%M:%S"`" 1>&2' >> $sseq_script

$action $sseq_script
