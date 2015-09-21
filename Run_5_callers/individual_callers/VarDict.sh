#!/bin/bash

# Usage: ./pipeline_02_somatic_callers.sh -o /path/to/results [-q]

# Just echo the cmd file, unless with -q flag then qsub it.

whattodo=echo
nthreads=20
path_vardict='/home/ltfang/apps/downloads/vardict'
hg_ref='/home/ltfang/references/human_g1k_v37_decoy.fasta'


# ---------- VarDict ---------- #
# VarDict takes FOREVER, so I'd like to really parallelize it instead of simply chromosome-by-chromosome:
wgsbed='/home/ltfang/oncopipe_cmds/snps/whole_genome_per5K.bed'

while getopts "o:n:t:s:g:V:q" opt
do
    case $opt in
        o)
            out_dir=$OPTARG;;
        q)
            whattodo=qsub;;
        n)
            normal_bam=$OPTARG;;
        t)
            tumor_bam=$OPTARG;;
	s)
	    nthreads=$OPTARG;;
	V)
	    path_vardict=$OPTARG;;
	g)
	    hg_ref=$OPTARG;;
	w)
	    wgsbed=$OPTARG;;
    esac
done


hg_dict=${hg_ref/.fa*/.dict}

timestamp=$( date +"%Y-%m-%d_%H-%M-%S" )

vardict_dir=${out_dir}/vardict_${timestamp}
mkdir ${vardict_dir}
log_dir=${vardict_dir}/logs
mkdir ${log_dir}



if ! [[ -d ${normal_bam} || -d ${tumor_bam} || -d ${log_dir} ]]
then
    echo "Missing ${normal_bam}, ${tumor_bam}, or ${log_dir}"
    exit 1
fi


nrows=$(wc -l $wgsbed | awk '{print $1}')
each_row=$(( $nrows / $nthreads ))

i=1
while [ $i -le $(( $nthreads )) ]
do

    # ---------- VarDict ---------- #
    out_script=No_${i}.vardict.${timestamp}.cmd

    echo '#!/bin/bash' > ${out_script}
    echo '' >> ${out_script}

    echo '#$ -S /bin/bash' >> ${out_script}
    echo "#\$ -o ${log_dir}" >> ${out_script}
    echo "#\$ -e ${log_dir}" >> ${out_script}
    echo ''  >> ${out_script}

    echo "export PATH=${path_samtools}:${path_vardict}:$PATH" >> ${out_script}
    echo "AF_THR=0.01  # minimum allele frequency" >> ${out_script}
    echo '' >> ${out_script}

    echo 'echo -e "Start VarDict at `date +"%Y/%m/%d %H:%M:%S"`\n" 1>&2' >> ${out_script}

    echo "cd ${vardict_dir}" >> $out_script


    if [ $i -lt $nthreads ]; then
        cat $wgsbed | awk -v start=$(( ($i-1) * $each_row + 1 )) -v end=$(( $i * $each_row )) 'NR>=start && NR<=end' > ${vardict_dir}/region.${i}.bed
    else
        cat $wgsbed | awk -v start=$(( ($i-1) * $each_row + 1 ))  'NR>=start' > ${vardict_dir}/region.${i}.bed
    fi

    echo '' >> $out_script

    echo "vardict -G ${hg_ref} -f \$AF_THR -h -b '${tumor_bam}|${normal_bam}' -C -Q 1 -c 1 -S 2 -E 3 -g 4 region.${i}.bed > region.${i}.var" >> $out_script
    echo "cat region.${i}.var | awk 'NR!=1' | testsomatic.R | var2vcf_paired.pl -N 'TUMOR|NORMAL' -f \$AF_THR > paired.${i}.vcf" >> $out_script

    echo '' >> $out_script

    echo 'echo -e "DONE VarDict at `date +"%Y/%m/%d %H:%M:%S"`\n" 1>&2' >> ${out_script}

    ${whattodo} ${out_script}

    i=$(( $i+1 ))  # Work on the next region

done
