#!/bin/bash

# Usage: ./pipeline_02_somatic_callers.sh -o /path/to/results [-q]

# Just echo the cmd file, unless with -q flag then qsub it. 

whattodo=echo
somaticsniper='/home/ltfang/apps/somatic-sniper/build/bin/bam-somaticsniper'

while getopts "o:N:T:e:q" opt
do
    case $opt in
        o)
            out_dir=$OPTARG;;
        q)
            whattodo=qsub;;
        N)
            normal_stuff=$OPTARG;;
        T)
            tumor_stuff=$OPTARG;;
	e)
	    somaticsniper=$OPTARG;;
    esac
done




hg_ref='/home/ltfang/references/genome.GRCh37.fa'
timestamp=$( date +"%Y-%m-%d_%H-%M-%S" )

sniper_dir=${out_dir}/somaticsniper_${timestamp}
log_dir=${sniper_dir}/logs

mkdir ${sniper_dir} ${log_dir}



for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT
do


    if [[ -e ${normal_bam} ]]
    then
        normal_bam_file=${normal_bam}
    else
        normal_bam_file=${normal_stuff}/${c}.bam
    fi


    if [[ -e ${tumor_bam} ]]
    then
        tumor_bam_file=${tumor_bam}
    else
        tumor_bam_file=${tumor_stuff}/${c}.bam
    fi

    # ---------- Somatic Sniper ---------- #
    out_script=chr${c}.sniper.${timestamp}.cmd

    echo '#!/bin/bash' > ${out_script}
    echo '' >> ${out_script}

    echo '#$ -S /bin/bash' >> ${out_script}
    echo "#\$ -o ${log_dir}" >> ${out_script}
    echo "#\$ -e ${log_dir}" >> ${out_script}
    echo ''  >> ${out_script}

    echo 'echo -e "Start running SomaticSniper at `date +"%Y/%m/%d %H:%M:%S"`\n" 1>&2' >> ${out_script}

    echo "${somaticsniper} -q 25 -Q 15 -s 0.0001 \\" >> ${out_script}
    echo "-f ${hg_ref} \\" >> ${out_script}
    echo "${tumor_bam_file}  \\" >> ${out_script}
    echo "${normal_bam_file} \\" >> ${out_script}
    echo "-F vcf ${sniper_dir}/${c}.vcf" >> ${out_script}
    echo "" >> ${out_script}

    echo 'echo -e "SomaticSniper Done at `date +"%Y/%m/%d %H:%M:%S"`\n" 1>&2' >> ${out_script}

    ${whattodo} ${out_script}


done
