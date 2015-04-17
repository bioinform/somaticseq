#!/bin/bash

# Usage: ./pipeline_02_somatic_callers.sh -o /path/to/results [-q]

# Just echo the cmd file, unless with -q flag then qsub it. 

whattodo=echo

while getopts "o:N:T:n:t:q" opt
do
    case $opt in
        o)
            out_dir=$OPTARG
            normal_stuff=${out_dir}/normal_stuff
            tumor_stuff=${out_dir}/tumor_stuff;;

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
    esac
done


hg_ref='/home/ltfang/references/human_g1k_v37_decoy.fasta'
hg_dict=${hg_ref/.fa*/.dict}
pon='/home/ltfang/apps/CancerAnalysisPackage-2013.2-18-g8207e53/reference/refseq_exome_10bp_hg19_300_1kg_normal_panel.vcf'
cosmic='/home/ltfang/references/cosmic.b37.vcf'
dbsnp='/home/ltfang/references/dbsnp_138.b37.vcf'


timestamp=$( date +"%Y-%m-%d_%H-%M-%S" )

log_dir=${out_dir}/logs

if ! [[ -d ${normal_stuff} || -d ${tumor_stuff} || -d ${log_dir} ]]
then
    echo "Missing ${normal_stuff}, ${tumor_stuff}, or ${log_dir}"
    exit 1
fi



sniper_dir=${out_dir}/somaticsniper
mkdir ${sniper_dir}

varscan_dir=${out_dir}/varscan2
mkdir ${varscan_dir}

snvmix_dir=${out_dir}/jointsnvmix2
mkdir ${snvmix_dir}

mutect_dir=${out_dir}/mutect
mutect_snp_dir=${mutect_dir}/snp
mutect_indel_dir=${mutect_dir}/indel
mkdir ${mutect_dir} ${mutect_snp_dir} ${mutect_indel_dir}

vardict_dir=${out_dir}/vardict
mkdir ${vardict_dir}


#--- LOCATION OF PROGRAMS ------
samtools='/home/ltfang/apps/samtools/samtools'
bgzip='/home/ltfang/apps/tabix/bgzip'
somaticsniper='/home/ltfang/apps/somatic-sniper/build/bin/bam-somaticsniper'
varscan='/home/ltfang/apps/varscan.jar'
jointsnvmix='/net/kodiak/volumes/lake/shared/opt/JointSNVMix-0.7.5/build/scripts-2.7/jsm.py'
snvconfigs='/net/kodiak/volumes/lake/shared/opt/JointSNVMix-0.7.5/config'
mutect='/home/ltfang/apps/CancerAnalysisPackage-2014.1-13-g6b71cb4/SomaticAnalysisTK.jar'


#--- LOCATION OF PATHS ------
path_samtools='/home/ltfang/apps/samtools'
path_vardict='/home/ltfang/apps/VarDict'


for c in 1 2 3 4 5 6 7 8 9 10 X 11 12 13 14 15 16 17 18 19 20 21 22 Y MT
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
    echo "--input_file:normal ${normal_stuff}/${c}.bam \\" >> ${out_script}
    echo "--input_file:tumor  ${tumor_stuff}/${c}.bam \\" >> ${out_script}
    echo "--normal_panel ${pon} \\" >> ${out_script}
    echo "--cosmic ${cosmic} \\" >> ${out_script}
    echo "--dbsnp ${dbsnp} \\" >> ${out_script}
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
    echo "--input_file:normal ${normal_stuff}/${c}.bam \\" >> ${out_script}
    echo "--input_file:tumor  ${tumor_stuff}/${c}.bam \\" >> ${out_script}
    echo "--out ${mutect_indel_dir}/${c}.vcf" >> ${out_script}
    echo "" >> ${out_script}

    echo 'echo -e "Indelocator done at `date +"%Y/%m/%d %H:%M:%S"`\n" 1>&2' >> ${out_script}

    ${whattodo} ${out_script}




    # ---------- JointSNVMix2 ---------- #
    # Train first:
    out_script=chr${c}.snvmix2.${timestamp}.cmd

    echo '#!/bin/bash' > ${out_script}
    echo '' >> ${out_script}

    echo '#$ -S /bin/bash' >> ${out_script}
    echo "#\$ -o ${log_dir}" >> ${out_script}
    echo "#\$ -e ${log_dir}" >> ${out_script}
    echo ''  >> ${out_script}

    echo 'echo -e "Start Training JointSNVMix2 at `date +"%Y/%m/%d %H:%M:%S"`\n" 1>&2'  >> $out_script

    echo "python ${jointsnvmix} train joint_snv_mix_two \\"  >> ${out_script}
    echo "--convergence_threshold 0.01 \\" >> ${out_script}
    echo "${hg_ref} \\"  >> ${out_script}
    echo "${normal_stuff}/${c}.bam \\"  >> ${out_script}
    echo "${tumor_stuff}/${c}.bam \\"  >> ${out_script}
    echo "${snvconfigs}/joint_priors.cfg \\"  >> ${out_script}
    echo "${snvconfigs}/joint_params.cfg \\"  >> ${out_script}
    echo "${tumor_stuff}/${c}.parameter.cfg" >> ${out_script}
    echo ""  >> ${out_script}

    # Classify second:
    # Run classify:
    echo 'echo -e "Start Classifying JointSNVMix2 at `date +"%Y/%m/%d %H:%M:%S"`\n" 1>&2'  >> $out_script

    echo "echo -e '##fileformat=VCFv4.1' > ${snvmix_dir}/${c}.vcf" >> ${out_script}
    echo "echo -e '##INFO=<ID=AAAB,Number=1,Type=Float,Description=\"Probability of Joint Genotype AA in Normal and AB in Tumor\">' >> ${snvmix_dir}/${c}.vcf" >> ${out_script}
    echo "echo -e '##INFO=<ID=AABB,Number=1,Type=Float,Description=\"Probability of Joint Genotype AA in Normal and BB in Tumor\">' >> ${snvmix_dir}/${c}.vcf" >> ${out_script}
    echo "echo -e '##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases (reads1)\">' >> ${snvmix_dir}/${c}.vcf" >> ${out_script}
    echo "echo -e '##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases (reads2)\">' >> ${snvmix_dir}/${c}.vcf" >> ${out_script}
    echo "echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR' >> ${snvmix_dir}/${c}.vcf" >> ${out_script}
    echo '' >> ${out_script}

    echo "python ${jointsnvmix} classify joint_snv_mix_two \\" >> ${out_script}
    echo "${hg_ref} \\" >> ${out_script}
    echo "${normal_stuff}/${c}.bam \\" >> ${out_script}
    echo "${tumor_stuff}/${c}.bam \\" >> ${out_script}
    echo "${tumor_stuff}/${c}.parameter.cfg \\" >> ${out_script}
    echo "/dev/stdout | awk -F \"\t\" 'NR!=1 && \$4!=\"N\" && \$10+\$11>=0.95' | \\" >> ${out_script}
    echo "awk -F \"\t\" '{print \$1 \"\t\" \$2 \"\t.\t\" \$3 \"\t\" \$4 \"\t.\t.\tAAAB=\" \$10 \";AABB=\" \$11 \"\tRD:AD\t\" \$5 \":\" \$6 \"\t\" \$7 \":\" \$8}' | \\" >> ${out_script}
    echo "sort -k1,1V -k2,2n >> ${snvmix_dir}/${c}.vcf" >> ${out_script}
    echo "" >> ${out_script}

    echo 'echo -e "JointJNVMix2 Done at `date +"%Y/%m/%d %H:%M:%S"`\n" 1>&2' >> ${out_script}

    ${whattodo} ${out_script}





    # ---------- SomaticSniper ---------- #
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
    echo "${tumor_stuff}/${c}.bam \\" >> ${out_script}
    echo "${normal_stuff}/${c}.bam \\" >> ${out_script}
    echo "-F vcf ${sniper_dir}/${c}.vcf" >> ${out_script}
    echo "" >> ${out_script}

    echo 'echo -e "SomaticSniper Done at `date +"%Y/%m/%d %H:%M:%S"`\n" 1>&2' >> ${out_script}

    ${whattodo} ${out_script}




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


    echo "${samtools} mpileup -B -q 25 -Q 20 -f \\" >> ${out_script}
    echo "${hg_ref} \\" >> ${out_script}
    echo "${normal_stuff}/${c}.bam \\" >> ${out_script}
    echo "> ${normal_stuff}/${c}.pileup.fifo &" >> ${out_script}
    echo "" >> ${out_script}

    echo "${samtools} mpileup -B -q 25 -Q 20 -f \\" >> ${out_script}
    echo "${hg_ref} \\" >> ${out_script}
    echo "${tumor_stuff}/${c}.bam \\" >> ${out_script}
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





# ---------- VarDict ---------- #
# VarDict takes FOREVER, so I'd like to really parallelize it instead of simply chromosome-by-chromosome:


wgsbed='/home/ltfang/oncopipe_cmds/snps/whole_genome_per5K.bed'


for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT
do


    # ---------- Indelocator ---------- #
    out_script=chr${c}.vardict.${timestamp}.cmd

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


    cat $wgsbed | egrep "^${c}\s" > ${vardict_dir}/${c}.bed

    echo '' >> $out_script

    echo "vardict -G ${hg_ref} -f \$AF_THR -h -b '${tumor_stuff}/${c}.bam|${normal_stuff}/${c}.bam' -z -F -C -c 1 -S 2 -E 3 -g 4 ${c}.bed > ${c}.var" >> $out_script
    echo "cat ${c}.var | awk 'NR!=1' | testsomatic.R | var2vcf_somatic.pl -N 'TUMOR|NORMAL' -f \$AF_THR > ${c}.vcf" >> $out_script

    echo '' >> $out_script



    echo "" >> ${out_script}
    echo 'echo -e "VarDict done at `date +"%Y/%m/%d %H:%M:%S"`\n" 1>&2' >> ${out_script}


    ${whattodo} ${out_script}


done
