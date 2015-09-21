#!/bin/bash

# Usage: ./pipeline_02_somatic_callers.sh -o /path/to/results [-q]

# Just echo the cmd file, unless with -q flag then qsub it.

whattodo=echo

hg_ref='/home/ltfang/references/human_g1k_v37_decoy.fasta'
jointsnvmix='/net/kodiak/volumes/lake/shared/opt/JointSNVMix-0.7.5/build/scripts-2.7/jsm.py'
snvconfigs='/net/kodiak/volumes/lake/shared/opt/JointSNVMix-0.7.5/config'

while getopts "o:N:T:g:j:c:q" opt
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
	    jointsnvmix=$OPTARG;;
	c)
	    snvconfigs=$OPTARG;;
    esac
done

timestamp=$( date +"%Y-%m-%d_%H-%M-%S" )

if ! [[ -d ${normal_stuff} || -d ${tumor_stuff} || -d ${log_dir} ]]
then
    echo "Missing ${normal_stuff}, ${tumor_stuff}, or ${log_dir}"
    exit 1
fi


snvmix_dir=${out_dir}/jsm2_${timestamp}
mkdir ${snvmix_dir}

log_dir=${snvmix_dir}/logs
mkdir $log_dir


for c in 1 2 3 4 5 6 7 8 9 10 X 11 12 13 14 15 16 17 18 19 20 21 22 Y MT
do

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

done
