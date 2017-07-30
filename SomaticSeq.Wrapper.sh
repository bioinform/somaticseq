	$MYDIR/utilities/modify_MuTect.py -type snp -infile ${mutect_vcf} -outfile ${merged_dir}/mutect.snp.vcf -nbam ${nbam} -tbam ${tbam}
	$MYDIR/utilities/modify_MuTect.py -type indel -infile ${indelocator_vcf} -outfile ${merged_dir}/indelocator.vcf -nbam ${nbam} -tbam ${tbam}
	$MYDIR/utilities/modify_VJSD.py -method SomaticSniper -infile ${sniper_vcf} -outfile ${merged_dir}/somaticsniper.vcf
	$MYDIR/utilities/modify_VJSD.py -method JointSNVMix2  -infile ${jsm_vcf} -outfile ${merged_dir}/jsm.vcf
	$MYDIR/utilities/modify_VJSD.py -method VarScan2 -infile ${varscan_vcf} -outfile ${merged_dir}/varscan2.snp.vcf
	$MYDIR/utilities/modify_VJSD.py -method MuSE  -infile ${muse_vcf} -outfile ${merged_dir}/muse.vcf
	$MYDIR/utilities/modify_VJSD.py -method VarScan2 -infile ${varscan_indel_vcf} -outfile ${merged_dir}/varscan2.indel.vcf
	$MYDIR/utilities/modify_VJSD.py -method VarDict -infile ${vardict_vcf} -outfile ${merged_dir}/vardict.vcf -filter paired
