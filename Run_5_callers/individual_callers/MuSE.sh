#!/bin/bash

# Usage: ./pipeline_02_somatic_callers.sh -o /path/to/results [-q]
# Just echo the cmd file, unless with -q flag then qsub it.

whattodo=echo
hg_ref='/home/ltfang/references/genome.GRCh37.fa'

muse_dir='/home/ltfang/apps/downloads/MuSE'
muse_exe='MuSEv1.0rc_submission_c039ffa'
muse=${muse_dir}/${muse_exe}


while getopts "o:N:T:g:e:q" opt
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
	g)
	    hg_ref=$OPTARG;;
	e)
	    muse=$OPTARG;;
    esac
done


hg_fai=${hg_ref}.fai

timestamp=$( date +"%Y-%m-%d_%H-%M-%S" )



muse_dir=${out_dir}/muse_${timestamp}
log_dir=${muse_dir}/logs

mkdir ${muse_dir} ${log_dir}


for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT
do

    out_script=chr${c}.muse.${timestamp}.cmd

    cat ${hg_fai} | awk -F "\t" -v chrom=${c} '$1==chrom' | awk -F "\t" '{print $1 "\t1\t" $2}' > ${muse_dir}/${c}.bed

    echo '#!/bin/bash' > ${out_script}
    echo '' >> ${out_script}

    echo '#$ -S /bin/bash' >> ${out_script}
    echo "#\$ -o ${log_dir}" >> ${out_script}
    echo "#\$ -e ${log_dir}" >> ${out_script}
    echo ''  >> ${out_script}

    echo 'echo -e "Start running MuSE at `date +"%Y/%m/%d %H:%M:%S"`\n" 1>&2' >> ${out_script}

    echo "${muse} call -O ${muse_dir}/${c} \\" >> ${out_script}
    echo "-l ${muse_dir}/${c}.bed \\" >> ${out_script}
    echo "-f ${hg_ref} ${tumor_bam_file} ${normal_bam_file}" >> ${out_script}
    echo "" >> ${out_script}

    echo 'echo -e "MuSE Done at `date +"%Y/%m/%d %H:%M:%S"`\n" 1>&2' >> ${out_script}

    ${whattodo} ${out_script}


done
