#!/usr/bin/env python3

import os, re, subprocess
import somaticseq.vcfModifier.copy_TextFile as copy_TextFile
import somaticseq.vcfModifier.splitVcf as splitVcf
import somaticseq.vcfModifier.getUniqueVcfPositions as getUniqueVcfPositions
from somaticseq.vcfModifier.vcfIntersector import *



# Combine individual VCF output into a simple combined VCF file, for single-sample callers
def combineSingle(outdir, ref, bam, inclusion=None, exclusion=None, mutect=None, mutect2=None, varscan=None, vardict=None, lofreq=None, scalpel=None, strelka=None, keep_intermediates=False):
    
    hg_dict = re.sub(r'\.fa(sta)?$', '.dict', ref)
    
    intermediate_files  = set()
    snv_intermediates   = []
    indel_intermediates = []

    intermediate_vcfs = {'MuTect2'  :{'snv': None, 'indel': None}, \
                         'VarScan2' :{'snv': None, 'indel': None}, \
                         'VarDict'  :{'snv': None, 'indel': None}, \
                         'LoFreq'   :{'snv': None, 'indel': None}, \
                         'Strelka'  :{'snv': None, 'indel': None}, }

    if mutect:
        
        import somaticseq.vcfModifier.modify_MuTect as mod_mutect
        
        mutect_in = bed_intersector(mutect, os.sep.join(( outdir, 'intersect.mutect1.vcf' )), inclusion, exclusion)
        intermediate_files.add(mutect_in)
        
        snv_mutect_out = os.sep.join(( outdir, 'snv.mutect1.vcf' ))
        mod_mutect.convert(mutect_in, snv_mutect_out, bam)
        
        intermediate_files.add(snv_mutect_out)
        snv_intermediates.append(snv_mutect_out)

    if mutect2:
        import somaticseq.vcfModifier.modify_ssMuTect2 as mod_mutect2
                
        mutect2_in = bed_intersector(mutect2, os.sep.join(( outdir, 'intersect.mutect2.vcf' )), inclusion, exclusion)
        intermediate_files.add(mutect2_in)
        
        snv_mutect_out   = os.sep.join(( outdir, 'snv.mutect2.vcf' ))
        indel_mutect_out = os.sep.join(( outdir, 'indel.mutect2.vcf' ))
        mod_mutect2.convert(mutect2_in, snv_mutect_out, indel_mutect_out)
        
        for file_i in snv_mutect_out, indel_mutect_out:
            intermediate_files.add( file_i )
            
        snv_intermediates.append(snv_mutect_out)
        indel_intermediates.append(indel_mutect_out)
        intermediate_vcfs['MuTect2']['snv']   = snv_mutect_out
        intermediate_vcfs['MuTect2']['indel'] = indel_mutect_out
        
    if varscan:
        import somaticseq.vcfModifier.modify_VarScan2 as mod_varscan2

        varscan_in = bed_intersector(varscan, os.sep.join(( outdir, 'intersect.varscan.vcf' )), inclusion, exclusion)
        intermediate_files.add(varscan_in)

        snv_temp          = os.sep.join(( outdir, 'snv.varscan.temp.vcf' ))
        indel_temp        = os.sep.join(( outdir, 'indel.varscan.temp.vcf' ))
        snv_varscan_out   = os.sep.join(( outdir, 'snv.varscan.vcf' ))
        indel_varscan_out = os.sep.join(( outdir, 'indel.varscan.vcf' ))

        splitVcf.split_into_snv_and_indel(varscan_in, snv_temp, indel_temp)
        mod_varscan2.convert(snv_temp,   snv_varscan_out)
        mod_varscan2.convert(indel_temp, indel_varscan_out)

        for file_i in snv_temp, indel_temp, snv_varscan_out, indel_varscan_out:
            intermediate_files.add( file_i )
        
        snv_intermediates.append(snv_varscan_out)
        indel_intermediates.append(indel_varscan_out)
        intermediate_vcfs['VarScan2']['snv']   = snv_varscan_out
        intermediate_vcfs['VarScan2']['indel'] = indel_varscan_out
        
    if vardict:
        import somaticseq.vcfModifier.modify_VarDict as mod_vardict
        
        # If the VarDict VCF file has line that clash with bedtools
        cleaned_vardict = os.sep.join(( outdir, 'cleaned.vardict.vcf' ))
        cleaned_vardict = remove_vcf_illegal_lines(vardict, cleaned_vardict)
        if cleaned_vardict:
            intermediate_files.add( cleaned_vardict )
        else:
            cleaned_vardict = vardict
        
        vardict_in = bed_intersector(cleaned_vardict, os.sep.join(( outdir, 'intersect.vardict.vcf' )), inclusion, exclusion)
        intermediate_files.add(vardict_in)

        snv_vardict_out   = os.sep.join(( outdir, 'snv.vardict.vcf' ))
        indel_vardict_out = os.sep.join(( outdir, 'indel.vardict.vcf'))
        mod_vardict.convert(vardict_in, snv_vardict_out, indel_vardict_out)
        
        sorted_snv_vardict_out   = os.sep.join(( outdir, 'snv.sort.vardict.vcf'))
        sorted_indel_vardict_out = os.sep.join(( outdir, 'indel.sort.vardict.vcf'))
        
        vcfsorter(ref, snv_vardict_out,   sorted_snv_vardict_out)
        vcfsorter(ref, indel_vardict_out, sorted_indel_vardict_out)
                
        for file_i in snv_vardict_out, indel_vardict_out, sorted_snv_vardict_out, sorted_indel_vardict_out:
            intermediate_files.add( file_i )
        
        snv_intermediates.append(sorted_snv_vardict_out)
        indel_intermediates.append(sorted_indel_vardict_out)
        intermediate_vcfs['VarDict']['snv']   = sorted_snv_vardict_out
        intermediate_vcfs['VarDict']['indel'] = sorted_indel_vardict_out

    if lofreq:
        
        lofreq_in = bed_intersector(lofreq, os.sep.join(( outdir, 'intersect.lofreq.vcf' )), inclusion, exclusion)
        intermediate_files.add(lofreq_in)
        
        snv_lofreq_out   = os.sep.join(( outdir, 'snv.lofreq.vcf' ))
        indel_lofreq_out = os.sep.join(( outdir, 'indel.lofreq.vcf' ))

        splitVcf.split_into_snv_and_indel(lofreq_in, snv_lofreq_out, indel_lofreq_out)

        for file_i in snv_lofreq_out, indel_lofreq_out:
            intermediate_files.add( file_i )
        
        snv_intermediates.append(snv_lofreq_out)
        indel_intermediates.append(indel_lofreq_out)
        intermediate_vcfs['LoFreq']['snv']   = snv_lofreq_out
        intermediate_vcfs['LoFreq']['indel'] = indel_lofreq_out

    if scalpel:
        
        scalpel_in = bed_intersector(scalpel, os.sep.join(( outdir, 'intersect.scalpel.vcf' )), inclusion, exclusion)
        intermediate_files.add(scalpel_in)

        scalpel_out = os.sep.join(( outdir, 'indel.scalpel.vcf' ))
        copy_TextFile.copy(scalpel_in, scalpel_out)
        intermediate_files.add(scalpel_out)
        indel_intermediates.append(scalpel_out)
        
    if strelka:
        import somaticseq.vcfModifier.modify_ssStrelka as mod_strelka
        
        strelka_in = bed_intersector(strelka, os.sep.join(( outdir, 'intersect.strelka.vcf' )), inclusion, exclusion)
        intermediate_files.add(strelka_in)

        snv_strelka_out   = os.sep.join(( outdir, 'snv.strelka.vcf' ))
        indel_strelka_out = os.sep.join(( outdir, 'indel.strelka.vcf' ))
        mod_strelka.convert(strelka_in, snv_strelka_out, indel_strelka_out)
        
        for file_i in snv_strelka_out, indel_strelka_out:
            intermediate_files.add( file_i )
        
        snv_intermediates.append(snv_strelka_out)
        indel_intermediates.append(indel_strelka_out)
        intermediate_vcfs['Strelka']['snv']   = snv_strelka_out
        intermediate_vcfs['Strelka']['indel'] = indel_strelka_out


    # Combine SNV/INDEL variant candidates
    snv_combined   = os.sep.join(( outdir, 'unsorted.CombineVariants.snv.vcf' ))
    indel_combined = os.sep.join(( outdir, 'unsorted.CombineVariants.indel.vcf' ))
    
    getUniqueVcfPositions.combine(snv_intermediates, snv_combined)
    getUniqueVcfPositions.combine(indel_intermediates, indel_combined)
    for file_i in snv_combined, indel_combined:
        intermediate_files.add( file_i )
    
    # Sort them:
    snv_combined_sorted = os.sep.join(( outdir, 'CombineVariants.snv.vcf' ))
    indel_combined_sorted = os.sep.join(( outdir, 'CombineVariants.indel.vcf' ))
    
    vcfsorter(ref, snv_combined,   snv_combined_sorted)
    vcfsorter(ref, indel_combined, indel_combined_sorted)
    

    if not keep_intermediates:
        for file_i in intermediate_files:
            subprocess.call( ('rm', '-v', file_i ) )
    
    return snv_combined_sorted, indel_combined_sorted, intermediate_vcfs, intermediate_files






# Combine individual VCF output into a simple combined VCF file, for paired sample callers
def combinePaired(outdir, ref, tbam, nbam, inclusion=None, exclusion=None, mutect=None, indelocator=None, mutect2=None, varscan_snv=None, varscan_indel=None, jsm=None, sniper=None, vardict=None, muse=None, lofreq_snv=None, lofreq_indel=None, scalpel=None, strelka_snv=None, strelka_indel=None, tnscope=None, platypus=None, keep_intermediates=False):
    
    hg_dict = re.sub(r'\.fa(sta)?$', '.dict', ref)
    
    intermediate_files  = set()
    snv_intermediates   = []
    indel_intermediates = []
    
    intermediate_vcfs = {'MuTect2':  {'snv': None, 'indel': None}, \
                         'VarDict':  {'snv': None, 'indel': None}, \
                         'TNscope':  {'snv': None, 'indel': None}, \
                         'Platypus': {'snv': None, 'indel': None} }
    
    # Modify direct VCF outputs for merging:
    if mutect or indelocator:
        
        import somaticseq.vcfModifier.modify_MuTect as mod_mutect

        if mutect:
            
            mutect_in = bed_intersector(mutect, os.sep.join(( outdir, 'intersect.mutect1.vcf' )), inclusion, exclusion)
            intermediate_files.add(mutect_in)
            
            snv_mutect_out = os.sep.join(( outdir, 'snv.mutect1.vcf' ))
            mod_mutect.convert(mutect_in, snv_mutect_out, tbam, nbam)
            
            intermediate_files.add(snv_mutect_out)
            snv_intermediates.append(snv_mutect_out)
        
        if indelocator:
            
            indelocator_in = bed_intersector(indelocator, os.sep.join(( outdir, 'intersect.indelocator.vcf' )), inclusion, exclusion)
            intermediate_files.add(indelocator_in)
            
            indel_indelocator_out = os.sep.join(( outdir, 'indel.indelocator.vcf'))
            mod_mutect.convert(indelocator_in, indel_indelocator_out, tbam, nbam)
            
            intermediate_files.add(indel_indelocator_out)
            indel_intermediates.append(indel_indelocator_out)

    
    if mutect2:
        
        import somaticseq.vcfModifier.modify_MuTect2 as mod_mutect2
                
        mutect2_in = bed_intersector(mutect2, os.sep.join(( outdir, 'intersect.mutect2.vcf' )), inclusion, exclusion)
        intermediate_files.add(mutect2_in)

        snv_mutect_out   = os.sep.join(( outdir, 'snv.mutect2.vcf'))
        indel_mutect_out = os.sep.join(( outdir, 'indel.mutect2.vcf'))
        mod_mutect2.convert(mutect2_in, snv_mutect_out, indel_mutect_out, False)
        
        for file_i in snv_mutect_out, indel_mutect_out:
            intermediate_files.add( file_i )
        
        snv_intermediates.append(snv_mutect_out)
        indel_intermediates.append(indel_mutect_out)
        
        intermediate_vcfs['MuTect2']['snv']   = snv_mutect_out
        intermediate_vcfs['MuTect2']['indel'] = indel_mutect_out
    
    if varscan_snv or varscan_indel:
        
        import somaticseq.vcfModifier.modify_VarScan2 as mod_varscan2
        
        if varscan_snv:
            
            varscan_in = bed_intersector(varscan_snv, os.sep.join(( outdir, 'intersect.varscan.snv.vcf' )), inclusion, exclusion)
            intermediate_files.add(varscan_in)
            
            snv_varscan_out = os.sep.join(( outdir, 'snv.varscan.vcf'))
            mod_varscan2.convert(varscan_in, snv_varscan_out)
            
            intermediate_files.add(snv_varscan_out)
            snv_intermediates.append(snv_varscan_out)
            
        if varscan_indel:

            varscan_in = bed_intersector(varscan_indel, os.sep.join(( outdir, 'intersect.varscan.indel.vcf' )), inclusion, exclusion)
            intermediate_files.add(varscan_in)
            
            indel_varscan_out = os.sep.join(( outdir, 'indel.varscan.vcf' ))
            mod_varscan2.convert(varscan_in, indel_varscan_out)
            
            intermediate_files.add(indel_varscan_out)
            indel_intermediates.append(indel_varscan_out)
    
    if jsm:
        import somaticseq.vcfModifier.modify_JointSNVMix2 as mod_jsm

        jsm_in = bed_intersector(jsm, os.sep.join(( outdir, 'intersect.jsm.vcf' )), inclusion, exclusion)
        intermediate_files.add(jsm_in)
        
        jsm_out = os.sep.join(( outdir, 'snv.jsm.vcf' ))
        mod_jsm.convert(jsm_in, jsm_out)
        
        intermediate_files.add(jsm_out)
        snv_intermediates.append(jsm_out)
    
    if sniper:
        import somaticseq.vcfModifier.modify_SomaticSniper as mod_sniper
        
        sniper_in = bed_intersector(sniper, os.sep.join(( outdir, 'intersect.sniper.vcf' )), inclusion, exclusion)
        intermediate_files.add(sniper_in)
        
        sniper_out = os.sep.join(( outdir, 'snv.somaticsniper.vcf' ))
        mod_sniper.convert(sniper_in, sniper_out)
        
        intermediate_files.add(sniper_out)
        snv_intermediates.append(sniper_out)
        
    if vardict:
        import somaticseq.vcfModifier.modify_VarDict as mod_vardict

        # If the VarDict VCF file has line that clash with bedtools
        cleaned_vardict = os.sep.join(( outdir, 'cleaned.vardict.vcf' ))
        cleaned_vardict = remove_vcf_illegal_lines(vardict, cleaned_vardict)
        if cleaned_vardict:
            intermediate_files.add( cleaned_vardict )
        else:
            cleaned_vardict = vardict

        vardict_in = bed_intersector(cleaned_vardict, os.sep.join(( outdir, 'intersect.vardict.vcf' )), inclusion, exclusion)
        intermediate_files.add(vardict_in)
        
        snv_vardict_out   = os.sep.join(( outdir, 'snv.vardict.vcf' ))
        indel_vardict_out = os.sep.join(( outdir, 'indel.vardict.vcf' ))
        mod_vardict.convert(vardict_in, snv_vardict_out, indel_vardict_out)
        
        sorted_snv_vardict_out   = os.sep.join(( outdir, 'snv.sort.vardict.vcf' ))
        sorted_indel_vardict_out = os.sep.join(( outdir, 'indel.sort.vardict.vcf' ))
        
        vcfsorter(ref, snv_vardict_out,   sorted_snv_vardict_out)
        vcfsorter(ref, indel_vardict_out, sorted_indel_vardict_out)
        
        for file_i in snv_vardict_out, indel_vardict_out, sorted_snv_vardict_out, sorted_indel_vardict_out:
            intermediate_files.add(file_i)
        
        snv_intermediates.append(sorted_snv_vardict_out)
        indel_intermediates.append(sorted_indel_vardict_out)
        intermediate_vcfs['VarDict']['snv']   = sorted_snv_vardict_out
        intermediate_vcfs['VarDict']['indel'] = sorted_indel_vardict_out
        
    if muse:

        muse_in = bed_intersector(muse, os.sep.join(( outdir, 'intersect.muse.vcf' )), inclusion, exclusion)
        intermediate_files.add(muse_in)
        
        muse_out = os.sep.join(( outdir, 'snv.muse.vcf' ))
        copy_TextFile.copy(muse_in, muse_out)
        
        intermediate_files.add(muse_out)
        snv_intermediates.append(muse_out)
        
    if lofreq_snv:
        
        lofreq_in = bed_intersector(lofreq_snv, os.sep.join(( outdir, 'intersect.lofreq.snv.vcf' )), inclusion, exclusion)
        intermediate_files.add(lofreq_in)


        snv_lofreq_out = os.sep.join(( outdir, 'snv.lofreq.vcf' ))
        copy_TextFile.copy(lofreq_in, snv_lofreq_out)
        
        intermediate_files.add(snv_lofreq_out)
        snv_intermediates.append(snv_lofreq_out)

    if lofreq_indel:
        
        lofreq_in = bed_intersector(lofreq_indel, os.sep.join(( outdir, 'intersect.lofreq.indel.vcf' )), inclusion, exclusion)
        intermediate_files.add(lofreq_in)

        indel_lofreq_out = os.sep.join(( outdir, 'indel.lofreq.vcf' ))
        copy_TextFile.copy(lofreq_in, indel_lofreq_out)
        
        intermediate_files.add(indel_lofreq_out)
        indel_intermediates.append(indel_lofreq_out)
        
    if scalpel:
        
        scalpel_in = bed_intersector(scalpel, os.sep.join(( outdir, 'intersect.scalpel.vcf' )), inclusion, exclusion)
        intermediate_files.add(scalpel_in)
        
        scalpel_out = os.sep.join(( outdir, 'indel.scalpel.vcf' ))
        copy_TextFile.copy(scalpel_in, scalpel_out)
        
        intermediate_files.add(scalpel_out)
        indel_intermediates.append(scalpel_out)
    
    if strelka_snv or strelka_indel:
        
        import somaticseq.vcfModifier.modify_Strelka as mod_strelka

        if strelka_snv:
            
            strelka_in = bed_intersector(strelka_snv, os.sep.join(( outdir, 'intersect.strelka.snv.vcf' )), inclusion, exclusion)
            intermediate_files.add(strelka_in)
            
            snv_strelka_out = os.sep.join(( outdir, 'snv.strelka.vcf' ))
            mod_strelka.convert(strelka_in, snv_strelka_out)
            
            intermediate_files.add(snv_strelka_out)
            snv_intermediates.append(snv_strelka_out)

        if strelka_indel:
            
            strelka_in = bed_intersector(strelka_indel, os.sep.join(( outdir, 'intersect.strelka.indel.vcf' )), inclusion, exclusion)
            intermediate_files.add(strelka_in)
                        
            indel_strelka_out = os.sep.join(( outdir, 'indel.strelka.vcf' ))
            mod_strelka.convert(strelka_in, indel_strelka_out)
            
            intermediate_files.add(indel_strelka_out)
            indel_intermediates.append(indel_strelka_out)
            
    if tnscope:

        import somaticseq.vcfModifier.modify_MuTect2 as mod_mutect2
        
        tnscope_in = bed_intersector(tnscope, os.sep.join(( outdir, 'intersect.tnscope.vcf' )), inclusion, exclusion)
        intermediate_files.add(tnscope_in)
        
        snv_tnscope_out   = os.sep.join(( outdir, 'snv.tnscope.vcf' ))
        indel_tnscope_out = os.sep.join(( outdir, 'indel.tnscope.vcf' ))
        mod_mutect2.convert(tnscope_in, snv_tnscope_out, indel_tnscope_out, True)
        
        for file_i in snv_tnscope_out, indel_tnscope_out:
            intermediate_files.add( file_i )
        
        snv_intermediates.append(snv_tnscope_out)
        indel_intermediates.append(indel_tnscope_out)
        intermediate_vcfs['TNscope']['snv']   = snv_tnscope_out
        intermediate_vcfs['TNscope']['indel'] = indel_tnscope_out
    
    if platypus:
        
        platypus_in = bed_intersector(platypus, os.sep.join(( outdir, 'intersect.platypus.vcf' )), inclusion, exclusion)
        intermediate_files.add(platypus_in)
        
        snv_platypus_out   = os.sep.join(( outdir, 'snv.platypus.vcf' ))
        indel_platypus_out = os.sep.join(( outdir, 'indel.platypus.vcf' ))

        splitVcf.split_into_snv_and_indel(platypus_in, snv_platypus_out, indel_platypus_out)

        for file_i in snv_platypus_out, indel_platypus_out:
            intermediate_files.add( file_i )
        
        snv_intermediates.append(snv_platypus_out)
        indel_intermediates.append(indel_platypus_out)
        intermediate_vcfs['Platypus']['snv']   = snv_platypus_out
        intermediate_vcfs['Platypus']['indel'] = indel_platypus_out

    
    
    # Combine SNV/INDEL variant candidates
    snv_combined   = os.sep.join(( outdir, 'unsorted.CombineVariants.snv.vcf' ))
    indel_combined = os.sep.join(( outdir, 'unsorted.CombineVariants.indel.vcf' ))
    
    getUniqueVcfPositions.combine(snv_intermediates, snv_combined)
    getUniqueVcfPositions.combine(indel_intermediates, indel_combined)
    
    for file_i in snv_combined, indel_combined:
        intermediate_files.add( file_i )
    
    # Sort them:
    snv_combined_sorted = os.sep.join(( outdir, 'CombineVariants.snv.vcf' ))
    indel_combined_sorted = os.sep.join(( outdir, 'CombineVariants.indel.vcf' ))
    
    vcfsorter(ref, snv_combined,   snv_combined_sorted)
    vcfsorter(ref, indel_combined, indel_combined_sorted)
    
    if not keep_intermediates:
        for file_i in intermediate_files:
            subprocess.call( ('rm', '-v', file_i ) )
    
    return snv_combined_sorted, indel_combined_sorted, intermediate_vcfs, intermediate_files

