#!/bin/bash

# Usage: ./pipeline_02_somatic_callers.sh -o /path/to/results [-q]

# Just echo the cmd file, unless with -q flag then qsub it.

whattodo=echo
pon='/home/ltfang/apps/CancerAnalysisPackage-2013.2-18-g8207e53/reference/refseq_exome_10bp_hg19_300_1kg_normal_panel.vcf'
cosmic='/home/ltfang/references/cosmic.b37.vcf'
dbsnp='/home/ltfang/references/dbsnp_138.b37.vcf'

hg_ref='/home/ltfang/references/genome.GRCh37.fa'
mutect='/home/ltfang/apps/CancerAnalysisPackage-2014.1-13-g6b71cb4/SomaticAnalysisTK.jar'

while getopts "o:n:t:g:e:d:c:p:q" opt
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
	g)
	    hg_ref=$OPTARG;;
	e)
	    mutect=$OPTARG;;
	d)
	    dbsnp=$OPTARG;;
	c)
	    cosmic=$OPTARG;;
	p)
	    pon=$OPTARG;;
    esac
done

timestamp=$( date +"%Y-%m-%d_%H-%M-%S" )


mutect_dir=${out_dir}/mutect_${timestamp}
mutect_snp_dir=${mutect_dir}/snp
mutect_indel_dir=${mutect_dir}/indel

log_dir=${mutect_dir}/logs


mkdir ${mutect_dir} ${mutect_snp_dir} ${mutect_indel_dir} ${log_dir}


for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 X 15 16 17 18 19 20 21 22 Y MT
do


    # ---------- MuTect ---------- #
    out_script=chr${c}.mutect.${timestamp}.cmd

    echo '#!/bin/bash' > ${out_script}
    echo '' >> ${out_script}

    echo '#$ -S /bin/bash' >> ${out_script}
    echo "#\$ -o ${log_dir}" >> ${out_script}
    echo "#\$ -e ${log_dir}" >> ${out_script}
    echo ''  >> ${out_script}

    echo 'echo -e "Start MuTect SNP at `date +"%Y/%m/%d %H:%M:%S"`\n" 1>&2' >> ${out_script}

    echo "java -Xmx4g -jar ${mutect} \\" >> ${out_script}
    echo "--analysis_type MuTect \\" >> ${out_script}
    echo "--reference_sequence ${hg_ref} \\" >> ${out_script}
    echo "-L ${c} \\" >> ${out_script}
    echo "--input_file:normal ${normal_bam_file} \\" >> ${out_script}
    echo "--input_file:tumor  ${tumor_bam_file} \\" >> ${out_script}
    echo "--normal_panel ${pon} \\" >> ${out_script}
    echo "--dbsnp ${dbsnp} \\" >> ${out_script}
    echo "--cosmic ${cosmic} \\" >> ${out_script}
    echo "--out /dev/null \\" >> ${out_script}
    echo "--vcf ${mutect_snp_dir}/${c}.vcf" >> ${out_script}
    echo "" >> ${out_script}

    echo 'echo -e "MuTect SNP done at `date +"%Y/%m/%d %H:%M:%S"`\n" 1>&2' >> ${out_script}

    ${whattodo} ${out_script}



    # ---------- Indelocator ---------- #
    out_script=chr${c}.indelocator.${timestamp}.cmd

    echo '#!/bin/bash' > ${out_script}
    echo '' >> ${out_script}

    echo '#$ -S /bin/bash' >> ${out_script}
    echo "#\$ -o ${log_dir}" >> ${out_script}
    echo "#\$ -e ${log_dir}" >> ${out_script}
    echo ''  >> ${out_script}


    echo 'echo -e "Start Indelocator at `date +"%Y/%m/%d %H:%M:%S"`\n" 1>&2'  >> ${out_script}

    echo "java -Xmx4g -jar ${mutect} \\" >> ${out_script}
    echo "--analysis_type SomaticIndelDetector \\" >> ${out_script}
    echo "--reference_sequence ${hg_ref} \\" >> ${out_script}
    echo "-L ${c} \\" >> ${out_script}
    echo "--input_file:normal ${normal_bam_file} \\" >> ${out_script}
    echo "--input_file:tumor  ${tumor_bam_file} \\" >> ${out_script}
    echo "--out ${mutect_indel_dir}/${c}.vcf" >> ${out_script}
    echo "" >> ${out_script}

    echo 'echo -e "Indelocator done at `date +"%Y/%m/%d %H:%M:%S"`\n" 1>&2' >> ${out_script}

    ${whattodo} ${out_script}

done
