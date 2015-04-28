#!/bin/bash

# Usage: ./pipeline_01_split_bams.sh -n normal.bam -t tumor.bam -o /path/to/results [-q]
# Optional: -N Normal_Split_Bam_Directory -T Tumor_Split_Bam_Directory, if you want overwrite the default locations. 


# Just echo the cmd file, unless with -q flag then qsub it.
whattodo=echo

while getopts "o:n:t:N:T:q" opt
do
    case $opt in
        o)
            out_dir=$OPTARG;;
        n)
            normal_bam=$OPTARG
            normal_stuff=${out_dir}/normal_stuff
            mkdir ${normal_stuff};;
        t)
            tumor_bam=$OPTARG
            tumor_stuff=${out_dir}/tumor_stuff
            mkdir ${tumor_stuff};;
        N)
            normal_stuff=$OPTARG
            normal_bam='DUMMY.NORMAL'
            ln -s ${normal_stuff} ${out_dir}/normal_stuff;;
        T)
            tumor_stuff=$OPTARG
            tumor_bam='DUMMY.TUMOR'
            ln -s ${tumor_stuff} ${out_dir}/tumor_stuff;;
        q)
            whattodo=qsub;;
    esac
done


hg_ref='/home/ltfang/references/human_g1k_v37_decoy.fasta'


if ! [[ -d ${out_dir} || -d ${hg_ref} ]];
then
    echo "Missing ${out_dir}, or ${hg_ref}"
    exit 1
fi


timestamp=$( date +"%Y-%m-%d_%H-%M-%S" )

log_dir=${out_dir}/logs
mkdir ${log_dir}


#--- LOCATION OF PROGRAMS ------
samtools='/home/ltfang/apps/samtools/samtools'
bgzip='/home/ltfang/apps/tabix/bgzip'
somaticsniper='/home/ltfang/apps/somatic-sniper/build/bin/bam-somaticsniper'
varscan='/home/ltfang/apps/varscan.jar'
jointsnvmix='/net/kodiak/volumes/lake/shared/opt/JointSNVMix-0.7.5/build/scripts-2.7/jsm.py'
snvconfigs='/net/kodiak/volumes/lake/shared/opt/JointSNVMix-0.7.5/config'
mutect='/home/ltfang/apps/CancerAnalysisPackage-2014.1-13-g6b71cb4/SomaticAnalysisTK.jar'



for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT
do

    tag='N'
    stuff=${normal_stuff}

    for bam_file in ${normal_bam} ${tumor_bam}
    do


        out_script=chr${c}.${tag}.split.${timestamp}.cmd

        echo '#!/bin/bash' > ${out_script}
        echo '' >> ${out_script}

        echo '#$ -S /bin/bash' >> ${out_script}
        echo "#\$ -o ${log_dir}" >> ${out_script}
        echo "#\$ -e ${log_dir}" >> ${out_script}
        echo ''  >> ${out_script}

        # Split the bam files, and then create pileups:
        echo 'echo -e "Start splitting bam files at `date +"%Y/%m/%d %H:%M:%S"`\n" 1>&2' >> ${out_script}
        echo "" >> ${out_script}

        echo "${samtools} view -bh ${bam_file} ${c} > ${stuff}/${c}.bam" >> ${out_script}

        echo ''  >> ${out_script}

        echo "if ! [[ -e ${stuff}/${c}.bam.bai || -e ${stuff}/${c}.bai ]]" >> ${out_script}
        echo "then" >> ${out_script}
        echo "    ${samtools} index ${stuff}/${c}.bam" >> ${out_script}
        echo "    ln ${stuff}/${c}.bam.bai ${stuff}/${c}.bai" >> ${out_script}
        echo "elif ! [ -e ${stuff}/${c}.bai ]" >> ${out_script}
        echo "then" >> ${out_script}
        echo "    ln ${stuff}/${c}.bam.bai ${stuff}/${c}.bai" >> ${out_script}
        echo "elif ! [ -e ${stuff}/${c}.bam.bai ]" >> ${out_script}
        echo "then" >> ${out_script}
        echo "    ln ${stuff}/${c}.bai ${stuff}/${c}.bam.bai" >> ${out_script}
        echo "fi" >> ${out_script}
        echo "" >> ${out_script}

        echo 'echo -e "Done splitting and indexing bam at `date +"%Y/%m/%d %H:%M:%S"`\n" 1>&2' >> ${out_script}

        ${whattodo} ${out_script}

        tag='T'
        stuff=${tumor_stuff}

    done

done
