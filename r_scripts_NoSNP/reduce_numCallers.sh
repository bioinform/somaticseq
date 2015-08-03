#!/bin/bash
#MVJSD

set -e

MYDIR="$( cd "$( dirname "$0" )" && pwd )"
export PATH=$MYDIR:/home/ltfang/apps/Bina_SomaticMerge:/net/kodiak/volumes/lake/shared/opt/python3/bin:$PATH

fai_file='/home/ltfang/references/human_g1k_v37_decoy.fasta.fai'

first_caller='CGA'
second_caller='VarScan2'
third_caller='JointSNVMix2'
forth_caller='SomaticSniper'
fifth_caller='VarDict'

while getopts "o:i:M:I:V:J:S:D:g:T:N:t:n:" opt
do
    case $opt in
        i)
            vcfin=$OPTARG;;
	M)
	    mutect_vcf=$OPTARG;;
	V)
	    varscan_vcf=$OPTARG;;
	J)
	    jsm_vcf=$OPTARG;;
	S)
	    sniper_vcf=$OPTARG;;
	D)
	    vardict_vcf=$OPTARG;;
	g)
	    fai_file=$OPTARG;;
	T)
	    haploT=$OPTARG;;
	N)
	    haploN=$OPTARG;;
	t)
	    samT=$OPTARG;;
	n)
	    samN=$OPTARG;;
    esac
done


for caller1 in $first_caller ''
do

    for caller2 in $second_caller ''
    do

        for caller3 in $third_caller ''
        do

            for caller4 in $forth_caller ''
            do

                for caller5 in $fifth_caller ''
                do

                    if [[ $caller1 == $first_caller ]]; then
                        c1=.${first_caller}
                    else
                        c1=''
                    fi

                    if [[ $caller2 == $second_caller ]]; then
                        c2=.${second_caller}
                        d2="-varscan $varscan_vcf"
                    else
                        c2=''
                        d2=''
                    fi

                    if [[ $caller3 == $third_caller ]]; then
                        c3=.${third_caller}
                        d3="-jsm $jsm_vcf"
                    else
                        c3=''
                        d3=''
                    fi

                    if [[ $caller4 == $forth_caller ]]; then
                        c4=.${forth_caller}
                        d4="-sniper $sniper_vcf"
                    else
                        c4=''
                        d4=''
                    fi

                    if [[ $caller5 == $fifth_caller ]]; then
                        c5=.${fifth_caller}
                        d5="-vardict $vardict_vcf"
                    else
                        c5=''
                        d5=''
                    fi


                    # COMMAND HERE
                    python3 /home/ltfang/programming/NGS/reduce_numCallers.py -infile ${vcfin} -outfile reduced${c1}${c2}${c3}${c4}${c5}.vcf -tools $caller1 $caller2 $caller3 $caller4 $caller5
                    python3 /home/ltfang/apps/Bina_SomaticMerge/SSeq_merged.vcf2tsv.py -myvcf reduced${c1}${c2}${c3}${c4}${c5}.vcf $d2 $d3 $d4 $d5 -haploN $haploN -haploT $haploT -samN $samN -samT $samT -fai $fai_file -outfile reduced${c1}${c2}${c3}${c4}${c5}.tsv
                    /home/ltfang/shared_delta/data/published/SomaticSeq/somaticseq/r_scripts_NoSNP/cross_validate_reduced.sh reduced${c1}${c2}${c3}${c4}${c5}.tsv
                    # COMMAND HERE


                done
            done
        done
    done
done
