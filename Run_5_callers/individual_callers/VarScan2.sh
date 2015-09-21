#!/bin/bash

# Usage: ./pipeline_02_somatic_callers.sh -o /path/to/results [-q]

# Just echo the cmd file, unless with -q flag then qsub it.

#--- LOCATION OF PROGRAMS ------
varscan='/home/ltfang/apps/varscan.jar'
samtools='samtools'
hg_ref='/home/ltfang/references/genome.GRCh37.fa'

whattodo=echo


while getopts "o:g:N:T:n:t:e:q" opt
do
    case $opt in
        o)
            out_dir=$OPTARG;;
	g)
	    hg_ref=$OPTARG;;
        q)
            whattodo=qsub;;
        N)
            normal_stuff=$OPTARG;;
        T)
            tumor_stuff=$OPTARG;;
	n)
	    normal_bam=$OPTARG;;
	t)
	    tumor_bam=$OPTARG;;
	e)
	    varscan=$OPTARG;;
    esac
done


timestamp=$( date +"%Y-%m-%d_%H-%M-%S" )

varscan_dir=${out_dir}/varscan2_${timestamp}
log_dir=${varscan_dir}/logs

mkdir ${varscan_dir} ${log_dir}


for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 X 15 16 17 18 19 20 21 22 Y MT
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


    # ---------- VarScan2 ---------- #
    out_script=chr${c}.varscan2.${timestamp}.cmd

    echo '#!/bin/bash' > ${out_script}
    echo '' >> ${out_script}

    echo '#$ -S /bin/bash' >> ${out_script}
    echo "#\$ -o ${log_dir}" >> ${out_script}
    echo "#\$ -e ${log_dir}" >> ${out_script}
    echo ''  >> ${out_script}

    echo 'echo -e "Start running VarScan2 at `date +"%Y/%m/%d %H:%M:%S"`\n" 1>&2'  >> ${out_script}
    echo "" >> ${out_script}

    echo "mkfifo ${normal_stuff}/${c}.pileup.fifo ${tumor_stuff}/${c}.pileup.fifo" >> ${out_script}
    echo "" >> ${out_script}


    echo "${samtools} mpileup -B -q 25 -Q 20 -r ${c} -f \\" >> ${out_script}
    echo "${hg_ref} \\" >> ${out_script}
    echo "${normal_bam_file} \\" >> ${out_script}
    echo "> ${normal_stuff}/${c}.pileup.fifo &" >> ${out_script}
    echo "" >> ${out_script}

    echo "${samtools} mpileup -B -q 25 -Q 20 -r ${c} -f \\" >> ${out_script}
    echo "${hg_ref} \\" >> ${out_script}
    echo "${tumor_bam_file} \\" >> ${out_script}
    echo "> ${tumor_stuff}/${c}.pileup.fifo &" >> ${out_script}
    echo "" >> ${out_script}



    echo "java -jar ${varscan} somatic \\" >> ${out_script}
    echo "${normal_stuff}/${c}.pileup.fifo \\" >> ${out_script}
    echo "${tumor_stuff}/${c}.pileup.fifo \\" >> ${out_script}
    echo "${varscan_dir}/${c} --output-vcf 1" >> ${out_script}
    echo "" >> ${out_script}

    echo "rm ${normal_stuff}/${c}.pileup.fifo ${tumor_stuff}/${c}.pileup.fifo" >> ${out_script}
    echo "" >> ${out_script}

    echo 'echo -e "VarScan2 Done at `date +"%Y/%m/%d %H:%M:%S"`\n" 1>&2' >> ${out_script}

    ${whattodo} ${out_script}

done
