#!/usr/bin/env python3

import sys, os, argparse, gzip, re, subprocess

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomicFileHandler.genomic_file_handlers as genome
import vcfModifier.copy_TextFile as copy_TextFile
import vcfModifier.getUniqueVcfPositions as getUniqueVcfPositions
from vcfModifier.vcfIntersector import *





# Combine individual VCF output into a simple combined VCF file, for single-sample callers
def combineSingle(outdir, ref, bam, inclusion=None, exclusion=None, mutect=None, mutect2=None, varscan=None, vardict=None, lofreq=None, scalpel=None, strelka=None, keep_intermediates=False):
    
    import vcfModifier.splitVcf as splitVcf
    
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
        
        import vcfModifier.modify_MuTect as mod_mutect
            
        if exclusion:
            mutect_ex = bed_exclude(mutect, exclusion, outdir + os.sep + 'snv.mutect1.ex.vcf')
            intermediate_files.add(mutect_ex)
        else:
            mutect_ex = mutect
        
        if inclusion:
            mutect_in = bed_include(mutect_ex, inclusion, outdir + os.sep + 'snv.mutect1.in.vcf')
            intermediate_files.add(mutect_in)
        else:
            mutect_in = mutect_ex
        
        snv_mutect_out = outdir + os.sep + 'snv.mutect1.vcf'
        mod_mutect.convert(mutect_in, snv_mutect_out, bam)
        
        intermediate_files.add(snv_mutect_out)
        snv_intermediates.append(snv_mutect_out)

    if mutect2:
        import vcfModifier.modify_ssMuTect2 as mod_mutect2
        
        if exclusion:
            mutect2_ex = bed_exclude(mutect2, exclusion, outdir + os.sep + 'mutect.ex.vcf')
            intermediate_files.add(mutect2_ex)
        else:
            mutect2_ex = mutect2
        
        if inclusion:
            mutect2_in = bed_include(mutect2_ex, inclusion, outdir + os.sep + 'mutect.in.vcf')
            intermediate_files.add(mutect2_in)
        else:
            mutect2_in = mutect2_ex
        
        snv_mutect_out   = outdir + os.sep + 'snv.mutect.vcf'
        indel_mutect_out = outdir + os.sep + 'indel.mutect.vcf'
        mod_mutect2.convert(mutect2_in, snv_mutect_out, indel_mutect_out)
        
        for file_i in snv_mutect_out, indel_mutect_out:
            intermediate_files.add( file_i )
            
        snv_intermediates.append(snv_mutect_out)
        indel_intermediates.append(indel_mutect_out)
        intermediate_vcfs['MuTect2']['snv']   = snv_mutect_out
        intermediate_vcfs['MuTect2']['indel'] = indel_mutect_out
        
    if varscan:
        import vcfModifier.modify_VarScan2 as mod_varscan2

        if exclusion:
            varscan_ex = bed_exclude(varscan, exclusion, outdir + os.sep + 'varscan.ex.vcf')
            intermediate_files.add(varscan_ex)
        else:
            varscan_ex = varscan
        
        if inclusion:
            varscan_in = bed_include(varscan_ex, inclusion, outdir + os.sep + 'varscan.in.vcf')
            intermediate_files.add(varscan_in)
        else:
            varscan_in = varscan_ex

        snv_temp          = outdir + os.sep + 'snv.varscan.temp.vcf'
        indel_temp        = outdir + os.sep + 'indel.varscan.temp.vcf'
        snv_varscan_out   = outdir + os.sep + 'snv.varscan.vcf'
        indel_varscan_out = outdir + os.sep + 'indel.varscan.vcf'

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
        import vcfModifier.modify_VarDict as mod_vardict
        
        if exclusion:
            vardict_ex = bed_exclude(vardict, exclusion, outdir + os.sep + 'vardict.ex.vcf')
            intermediate_files.add(vardict_ex)
        else:
            vardict_ex = vardict
        
        if inclusion:
            vardict_in = bed_include(vardict_ex, inclusion, outdir + os.sep + 'vardict.in.vcf')
            intermediate_files.add(vardict_in)
        else:
            vardict_in = vardict_ex

        snv_vardict_out   = outdir + os.sep + 'snv.vardict.vcf'
        indel_vardict_out = outdir + os.sep + 'indel.vardict.vcf'
        mod_vardict.convert(vardict_in, snv_vardict_out, indel_vardict_out)
        
        sorted_snv_vardict_out   = outdir + os.sep + 'snv.sort.vardict.vcf'
        sorted_indel_vardict_out = outdir + os.sep + 'indel.sort.vardict.vcf'
        
        vcfsorter(ref, snv_vardict_out,   sorted_snv_vardict_out)
        vcfsorter(ref, indel_vardict_out, sorted_indel_vardict_out)
                
        for file_i in snv_vardict_out, indel_vardict_out, sorted_snv_vardict_out, sorted_indel_vardict_out:
            intermediate_files.add( file_i )
        
        snv_intermediates.append(sorted_snv_vardict_out)
        indel_intermediates.append(sorted_indel_vardict_out)
        intermediate_vcfs['VarDict']['snv']   = sorted_snv_vardict_out
        intermediate_vcfs['VarDict']['indel'] = sorted_indel_vardict_out

    if lofreq:
        
        if exclusion:
            lofreq_ex = bed_exclude(lofreq, exclusion, outdir + os.sep + 'lofreq.ex.vcf')
            intermediate_files.add(lofreq_ex)
        else:
            lofreq_ex = lofreq
        
        if inclusion:
            lofreq_in = bed_include(lofreq_ex, inclusion, outdir + os.sep + 'lofreq.in.vcf')
            intermediate_files.add(lofreq_in)
        else:
            lofreq_in = lofreq_ex
        
        snv_lofreq_out   = outdir + os.sep + 'snv.lofreq.vcf'
        indel_lofreq_out = outdir + os.sep + 'indel.lofreq.vcf'

        splitVcf.split_into_snv_and_indel(lofreq_in, snv_lofreq_out, indel_lofreq_out)

        for file_i in snv_lofreq_out, indel_lofreq_out:
            intermediate_files.add( file_i )
        
        snv_intermediates.append(snv_lofreq_out)
        indel_intermediates.append(indel_lofreq_out)
        intermediate_vcfs['LoFreq']['snv']   = snv_varscan_out
        intermediate_vcfs['LoFreq']['indel'] = indel_varscan_out

    if scalpel:
        
        if exclusion:
            scalpel_ex = bed_exclude(scalpel, exclusion, outdir + os.sep + 'scalpel.ex.vcf')
            intermediate_files.add(scalpel_ex)
        else:
            scalpel_ex = scalpel
        
        if inclusion:
            scalpel_in = bed_include(scalpel_ex, inclusion, outdir + os.sep + 'scalpel.in.vcf')
            intermediate_files.add(scalpel_in)
        else:
            scalpel_in = scalpel_ex

        scalpel_out = outdir + os.sep + 'indel.scalpel.vcf'
        copy_TextFile.copy(scalpel_in, scalpel_out)
        intermediate_files.add(scalpel_out)
        indel_intermediates.append(scalpel_out)
        
    if strelka:
        import vcfModifier.modify_ssStrelka as mod_strelka
        
        if exclusion:
            strelka_ex = bed_exclude(strelka, exclusion, outdir + os.sep + 'strelka.ex.vcf')
            intermediate_files.add(strelka_ex)
        else:
            strelka_ex = strelka
        
        if inclusion:
            strelka_in = bed_include(strelka_ex, inclusion, outdir + os.sep + 'strelka.in.vcf')
            intermediate_files.add(strelka_in)
        else:
            strelka_in = strelka_ex

        snv_strelka_out   = outdir + os.sep + 'snv.strelka.vcf'
        indel_strelka_out = outdir + os.sep + 'indel.strelka.vcf'
        mod_strelka.convert(strelka_in, snv_strelka_out, indel_strelka_out)
        
        for file_i in snv_strelka_out, indel_strelka_out:
            intermediate_files.add( file_i )
        
        snv_intermediates.append(snv_strelka_out)
        indel_intermediates.append(indel_strelka_out)
        intermediate_vcfs['Strelka']['snv']   = snv_strelka_out
        intermediate_vcfs['Strelka']['indel'] = indel_strelka_out


    # Combine SNV/INDEL variant candidates
    snv_combined   = outdir + os.sep + 'unsorted.CombineVariants.snv.vcf'
    indel_combined = outdir + os.sep + 'unsorted.CombineVariants.indel.vcf'
    
    getUniqueVcfPositions.combine(snv_intermediates, snv_combined)
    getUniqueVcfPositions.combine(indel_intermediates, indel_combined)
    for file_i in snv_combined, indel_combined:
        intermediate_files.add( file_i )
    
    # Sort them:
    snv_combined_sorted = outdir + os.sep + 'CombineVariants.snv.vcf'
    indel_combined_sorted = outdir + os.sep + 'CombineVariants.indel.vcf'
    
    vcfsorter(ref, snv_combined,   snv_combined_sorted)
    vcfsorter(ref, indel_combined, indel_combined_sorted)
    

    if not keep_intermediates:
        for file_i in intermediate_files:
            subprocess.call( ('rm', '-v', file_i ) )
    
    return snv_combined_sorted, indel_combined_sorted, intermediate_vcfs, intermediate_files






# Combine individual VCF output into a simple combined VCF file, for paired sample callers
def combinePaired(outdir, ref, tbam, nbam, inclusion=None, exclusion=None, mutect=None, indelocator=None, mutect2=None, varscan_snv=None, varscan_indel=None, jsm=None, sniper=None, vardict=None, muse=None, lofreq_snv=None, lofreq_indel=None, scalpel=None, strelka_snv=None, strelka_indel=None, tnscope=None, keep_intermediates=False):
    
    hg_dict = re.sub(r'\.fa(sta)?$', '.dict', ref)
    
    intermediate_files  = set()
    snv_intermediates   = []
    indel_intermediates = []
    
    intermediate_vcfs = {'MuTect2':{'snv': None, 'indel': None}, \
                         'VarDict':{'snv': None, 'indel': None}, \
                         'TNscope':{'snv': None, 'indel': None}, }
    
    # Modify direct VCF outputs for merging:
    if mutect or indelocator:
        
        import vcfModifier.modify_MuTect as mod_mutect

        if mutect:
            
            if exclusion:
                mutect_ex = bed_exclude(mutect, exclusion, outdir + os.sep + 'snv.mutect1.ex.vcf')
                intermediate_files.add(mutect_ex)
            else:
                mutect_ex = mutect
            
            if inclusion:
                mutect_in = bed_include(mutect_ex, inclusion, outdir + os.sep + 'snv.mutect1.in.vcf')
                intermediate_files.add(mutect_in)
            else:
                mutect_in = mutect_ex
            
            snv_mutect_out = outdir + os.sep + 'snv.mutect1.vcf'
            mod_mutect.convert(mutect_in, snv_mutect_out, tbam, nbam)
            
            intermediate_files.add(snv_mutect_out)
            snv_intermediates.append(snv_mutect_out)
        
        if indelocator:

            if exclusion:
                indelocator_ex = bed_exclude(indelocator, exclusion, outdir + os.sep + 'indel.indelocator.ex.vcf')
                intermediate_files.add(indelocator_ex)
            else:
                indelocator_ex = indelocator
            
            if inclusion:
                indelocator_in = bed_include(indelocator_ex, inclusion, outdir + os.sep + 'indel.indelocator.in.vcf')
                intermediate_files.add(indelocator_in)
            else:
                indelocator_in = indelocator_ex
            
            indel_indelocator_out = outdir + os.sep + 'indel.indelocator.vcf'
            mod_mutect.convert(indelocator_in, indel_indelocator_out, tbam, nbam)
            
            intermediate_files.add(indel_indelocator_out)
            indel_intermediates.append(indel_indelocator_out)

    
    if mutect2:
        
        import vcfModifier.modify_MuTect2 as mod_mutect2
        
        if exclusion:
            mutect2_ex = bed_exclude(mutect2, exclusion, outdir + os.sep + 'mutect.ex.vcf')
            intermediate_files.add(mutect2_ex)
        else:
            mutect2_ex = mutect2
        
        if inclusion:
            mutect2_in = bed_include(mutect2_ex, inclusion, outdir + os.sep + 'mutect.in.vcf')
            intermediate_files.add(mutect2_in)
        else:
            mutect2_in = mutect2_ex
        
        snv_mutect_out   = outdir + os.sep + 'snv.mutect.vcf'
        indel_mutect_out = outdir + os.sep + 'indel.mutect.vcf'
        mod_mutect2.convert(mutect2_in, snv_mutect_out, indel_mutect_out, False)
        
        for file_i in snv_mutect_out, indel_mutect_out:
            intermediate_files.add( file_i )
        
        snv_intermediates.append(snv_mutect_out)
        indel_intermediates.append(indel_mutect_out)
        
        intermediate_vcfs['MuTect2']['snv']   = snv_mutect_out
        intermediate_vcfs['MuTect2']['indel'] = indel_mutect_out
    
    if varscan_snv or varscan_indel:
                
        import vcfModifier.modify_VarScan2 as mod_varscan2
        
        if varscan_snv:
            
            if exclusion:
                varscan_ex = bed_exclude(varscan_snv, exclusion, outdir + os.sep + 'snv.varscan.ex.vcf')
                intermediate_files.add(varscan_ex)
            else:
                varscan_ex = varscan_snv
            
            if inclusion:
                varscan_in = bed_include(varscan_ex, inclusion, outdir + os.sep + 'snv.varscan.in.vcf')
                intermediate_files.add(varscan_in)
            else:
                varscan_in = varscan_ex
            
            snv_varscan_out = outdir + os.sep + 'snv.varscan.vcf'
            mod_varscan2.convert(varscan_in, snv_varscan_out)
            
            intermediate_files.add(snv_varscan_out)
            snv_intermediates.append(snv_varscan_out)
            
        if varscan_indel:

            if exclusion:
                varscan_ex = bed_exclude(varscan_indel, exclusion, outdir + os.sep + 'indel.varscan.ex.vcf')
                intermediate_files.add(varscan_ex)
            else:
                varscan_ex = varscan_indel
            
            if inclusion:
                varscan_in = bed_include(varscan_ex, inclusion, outdir + os.sep + 'indel.varscan.in.vcf')
                intermediate_files.add(varscan_in)
            else:
                varscan_in = varscan_ex
            
            indel_varscan_out = outdir + os.sep + 'indel.varscan.vcf'
            mod_varscan2.convert(varscan_in, indel_varscan_out)
            
            intermediate_files.add(indel_varscan_out)
            indel_intermediates.append(indel_varscan_out)
    
    if jsm:
        import vcfModifier.modify_JointSNVMix2 as mod_jsm

        if exclusion:
            jsm_ex = bed_exclude(jsm, exclusion, outdir + os.sep + 'snv.jsm.ex.vcf')
            intermediate_files.add(jsm_ex)
        else:
            jsm_ex = jsm
        
        if inclusion:
            jsm_in = bed_include(jsm_ex, inclusion, outdir + os.sep + 'snv.jsm.in.vcf')
            intermediate_files.add(jsm_in)
        else:
            jsm_in = jsm_ex
        
        jsm_out = outdir + os.sep + 'snv.jsm.vcf'
        mod_jsm.convert(jsm_in, jsm_out)
        
        intermediate_files.add(jsm_out)
        snv_intermediates.append(jsm_out)
        
    if sniper:
        import vcfModifier.modify_SomaticSniper as mod_sniper
        
        if exclusion:
            sniper_ex = bed_exclude(sniper, exclusion, outdir + os.sep + 'snv.somaticsniper.ex.vcf')
            intermediate_files.add(sniper_ex)
        else:
            sniper_ex = sniper
        
        if inclusion:
            sniper_in = bed_include(sniper_ex, inclusion, outdir + os.sep + 'snv.somaticsniper.in.vcf')
            intermediate_files.add(sniper_in)
        else:
            sniper_in = sniper_ex
        
        sniper_out = outdir + os.sep + 'snv.somaticsniper.vcf'
        mod_sniper.convert(sniper_in, sniper_out)
        
        intermediate_files.add(sniper_out)
        snv_intermediates.append(sniper_out)
        
    if vardict:
        import vcfModifier.modify_VarDict as mod_vardict

        if exclusion:
            vardict_ex = bed_exclude(vardict, exclusion, outdir + os.sep + 'vardict.ex.vcf')
            intermediate_files.add(vardict_ex)
        else:
            vardict_ex = vardict
        
        if inclusion:
            vardict_in = bed_include(vardict_ex, inclusion, outdir + os.sep + 'vardict.in.vcf')
            intermediate_files.add(vardict_in)
        else:
            vardict_in = vardict_ex
        
        snv_vardict_out   = outdir + os.sep + 'snv.vardict.vcf'
        indel_vardict_out = outdir + os.sep + 'indel.vardict.vcf'
        mod_vardict.convert(vardict_in, snv_vardict_out, indel_vardict_out)
        
        sorted_snv_vardict_out   = outdir + os.sep + 'snv.sort.vardict.vcf'
        sorted_indel_vardict_out = outdir + os.sep + 'indel.sort.vardict.vcf'
        
        vcfsorter(ref, snv_vardict_out,   sorted_snv_vardict_out)
        vcfsorter(ref, indel_vardict_out, sorted_indel_vardict_out)
        
        for file_i in snv_vardict_out, indel_vardict_out, sorted_snv_vardict_out, sorted_indel_vardict_out:
            intermediate_files.add(file_i)
        
        snv_intermediates.append(sorted_snv_vardict_out)
        indel_intermediates.append(sorted_indel_vardict_out)
        intermediate_vcfs['VarDict']['snv']   = sorted_snv_vardict_out
        intermediate_vcfs['VarDict']['indel'] = sorted_indel_vardict_out
        
    if muse:

        if exclusion:
            muse_ex = bed_exclude(muse, exclusion, outdir + os.sep + 'snv.muse.ex.vcf')
            intermediate_files.add(muse_ex)
        else:
            muse_ex = muse
        
        if inclusion:
            muse_in = bed_include(muse_ex, inclusion, outdir + os.sep + 'snv.muse.in.vcf')
            intermediate_files.add(muse_in)
        else:
            muse_in = muse_ex
        
        muse_out = outdir + os.sep + 'snv.muse.vcf'
        copy_TextFile.copy(muse_in, muse_out)
        
        intermediate_files.add(muse_out)
        snv_intermediates.append(muse_out)
        
    if lofreq_snv:
        
        if exclusion:
            lofreq_ex = bed_exclude(lofreq_snv, exclusion, outdir + os.sep + 'snv.lofreq.ex.vcf')
            intermediate_files.add(lofreq_ex)
        else:
            lofreq_ex = lofreq_snv
        
        if inclusion:
            lofreq_in = bed_include(lofreq_ex, inclusion, outdir + os.sep + 'snv.lofreq.in.vcf')
            intermediate_files.add(lofreq_in)
        else:
            lofreq_in = lofreq_ex

        snv_lofreq_out = outdir + os.sep + 'snv.lofreq.vcf'
        copy_TextFile.copy(lofreq_in, snv_lofreq_out)
        
        intermediate_files.add(snv_lofreq_out)
        snv_intermediates.append(snv_lofreq_out)

    if lofreq_indel:
        
        if exclusion:
            lofreq_ex = bed_exclude(lofreq_indel, exclusion, outdir + os.sep + 'indel.lofreq.ex.vcf')
            intermediate_files.add(lofreq_ex)
        else:
            lofreq_ex = lofreq_snv
        
        if inclusion:
            lofreq_in = bed_include(lofreq_ex, inclusion, outdir + os.sep + 'indel.lofreq.in.vcf')
            intermediate_files.add(lofreq_in)
        else:
            lofreq_in = lofreq_ex
        
        indel_lofreq_out = outdir + os.sep + 'indel.lofreq.vcf'
        copy_TextFile.copy(lofreq_in, indel_lofreq_out)
        
        intermediate_files.add(indel_lofreq_out)
        indel_intermediates.append(indel_lofreq_out)
        
    if scalpel:
        
        if exclusion:
            scalpel_ex = bed_exclude(scalpel, exclusion, outdir + os.sep + 'indel.scalpel.ex.vcf')
            intermediate_files.add(scalpel_ex)
        else:
            scalpel_ex = scalpel
        
        if inclusion:
            scalpel_in = bed_include(scalpel_ex, inclusion, outdir + os.sep + 'indel.scalpel.in.vcf')
            intermediate_files.add(scalpel_in)
        else:
            scalpel_in = scalpel_ex
        
        scalpel_out = outdir + os.sep + 'indel.scalpel.vcf'
        copy_TextFile.copy(scalpel_in, scalpel_out)
        
        intermediate_files.add(scalpel_out)
        indel_intermediates.append(scalpel_out)
    
    if strelka_snv or strelka_indel:
        
        import vcfModifier.modify_Strelka as mod_strelka

        if strelka_snv:
            
            if exclusion:
                strelka_ex = bed_exclude(strelka_snv, exclusion, outdir + os.sep + 'snv.strelka.ex.vcf')
                intermediate_files.add(strelka_ex)
            else:
                strelka_ex = strelka_snv
            
            if inclusion:
                strelka_in = bed_include(strelka_ex, inclusion, outdir + os.sep + 'snv.strelka.in.vcf')
                intermediate_files.add(strelka_in)
            else:
                strelka_in = strelka_ex
            
            snv_strelka_out = outdir + os.sep + 'snv.strelka.vcf'
            mod_strelka.convert(strelka_in, snv_strelka_out)
            
            intermediate_files.add(snv_strelka_out)
            snv_intermediates.append(snv_strelka_out)

        if strelka_indel:
            
            if exclusion:
                strelka_ex = bed_exclude(strelka_indel, exclusion, outdir + os.sep + 'indel.strelka.ex.vcf')
                intermediate_files.add(strelka_ex)
            else:
                strelka_ex = strelka_snv
            
            if inclusion:
                strelka_in = bed_include(strelka_ex, inclusion, outdir + os.sep + 'indel.strelka.in.vcf')
                intermediate_files.add(strelka_in)
            else:
                strelka_in = strelka_ex
                        
            indel_strelka_out = outdir + os.sep + 'indel.strelka.vcf'
            mod_strelka.convert(strelka_in, indel_strelka_out)
            
            intermediate_files.add(indel_strelka_out)
            indel_intermediates.append(indel_strelka_out)
            
    if tnscope:

        import vcfModifier.modify_MuTect2 as mod_mutect2
        
        if exclusion:
            tnscope_ex = bed_exclude(tnscope, exclusion, outdir + os.sep + 'tnscope.ex.vcf')
            intermediate_files.add(tnscope_ex)
        else:
            tnscope_ex = tnscope
        
        if inclusion:
            tnscope_in = bed_include(tnscope_ex, inclusion, outdir + os.sep + 'tnscope.in.vcf')
            intermediate_files.add(tnscope_in)
        else:
            tnscope_in = tnscope_ex
        
        snv_tnscope_out   = outdir + os.sep + 'snv.tnscope.vcf'
        indel_tnscope_out = outdir + os.sep + 'indel.tnscope.vcf'
        mod_mutect2.convert(tnscope_in, snv_tnscope_out, indel_tnscope_out, True)
        
        for file_i in snv_tnscope_out, indel_tnscope_out:
            intermediate_files.add( file_i )
        
        snv_intermediates.append(snv_tnscope_out)
        indel_intermediates.append(indel_tnscope_out)
        intermediate_vcfs['TNscope']['snv']   = snv_tnscope_out
        intermediate_vcfs['TNscope']['indel'] = indel_tnscope_out
    
    
    # Combine SNV/INDEL variant candidates
    snv_combined   = outdir + os.sep + 'unsorted.CombineVariants.snv.vcf'
    indel_combined = outdir + os.sep + 'unsorted.CombineVariants.indel.vcf'
    
    getUniqueVcfPositions.combine(snv_intermediates, snv_combined)
    getUniqueVcfPositions.combine(indel_intermediates, indel_combined)
    
    for file_i in snv_combined, indel_combined:
        intermediate_files.add( file_i )
    
    # Sort them:
    snv_combined_sorted = outdir + os.sep + 'CombineVariants.snv.vcf'
    indel_combined_sorted = outdir + os.sep + 'CombineVariants.indel.vcf'
    
    vcfsorter(ref, snv_combined,   snv_combined_sorted)
    vcfsorter(ref, indel_combined, indel_combined_sorted)
    
    if not keep_intermediates:
        for file_i in intermediate_files:
            subprocess.call( ('rm', '-v', file_i ) )
    
    return snv_combined_sorted, indel_combined_sorted, intermediate_vcfs, intermediate_files

