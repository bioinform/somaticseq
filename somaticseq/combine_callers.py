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
    
    intermediate_files  = []
    snv_intermediates   = []
    indel_intermediates = []

    intermediate_vcfs = {'MuTect2'  :{'snv': None, 'indel': None}, \
                         'VarScan2' :{'snv': None, 'indel': None}, \
                         'VarDict'  :{'snv': None, 'indel': None}, \
                         'LoFreq'   :{'snv': None, 'indel': None}, \
                         'Strelka'  :{'snv': None, 'indel': None}, }

    if mutect2:
        import vcfModifier.modify_ssMuTect2 as mod_mutect2
        snv_mutect_out   = outdir + os.sep + 'snv.mutect.vcf'
        indel_mutect_out = outdir + os.sep + 'indel.mutect.vcf'
        mod_mutect2.convert(mutect2, snv_mutect_out, indel_mutect_out)
        
        intermediate_files.extend( [snv_mutect_out, indel_mutect_out] )
        snv_intermediates.append(snv_mutect_out)
        indel_intermediates.append(indel_mutect_out)
        intermediate_vcfs['MuTect2']['snv']   = snv_mutect_out
        intermediate_vcfs['MuTect2']['indel'] = indel_mutect_out
        
    if varscan:
        import vcfModifier.modify_VarScan2 as mod_varscan2
        
        snv_temp          = outdir + os.sep + 'snv.varscan.temp.vcf'
        indel_temp        = outdir + os.sep + 'indel.varscan.temp.vcf'
        snv_varscan_out   = outdir + os.sep + 'snv.varscan.vcf'
        indel_varscan_out = outdir + os.sep + 'indel.varscan.vcf'

        splitVcf.split_into_snv_and_indel(varscan, snv_temp, indel_temp)
        mod_varscan2.convert(snv_temp,   snv_varscan_out)
        mod_varscan2.convert(indel_temp, indel_varscan_out)

        intermediate_files.extend( [snv_temp, indel_temp, snv_varscan_out, indel_varscan_out] )
        snv_intermediates.append(snv_varscan_out)
        indel_intermediates.append(indel_varscan_out)
        intermediate_vcfs['VarScan2']['snv']   = snv_varscan_out
        intermediate_vcfs['VarScan2']['indel'] = indel_varscan_out
        
    if vardict:
        import vcfModifier.modify_VarDict as mod_vardict
        
        snv_vardict_out   = outdir + os.sep + 'snv.vardict.vcf'
        indel_vardict_out = outdir + os.sep + 'indel.vardict.vcf'
        mod_vardict.convert(vardict, snv_vardict_out, indel_vardict_out)
        
        sorted_snv_vardict_out   = outdir + os.sep + 'snv.sort.vardict.vcf'
        sorted_indel_vardict_out = outdir + os.sep + 'indel.sort.vardict.vcf'
        
        vcfsorter(hg_dict, snv_vardict_out,   sorted_snv_vardict_out)
        vcfsorter(hg_dict, indel_vardict_out, sorted_indel_vardict_out)
                
        intermediate_files.extend( [snv_vardict_out, indel_vardict_out, sorted_snv_vardict_out, sorted_indel_vardict_out] )
        snv_intermediates.append(sorted_snv_vardict_out)
        indel_intermediates.append(sorted_indel_vardict_out)
        intermediate_vcfs['VarDict']['snv']   = sorted_snv_vardict_out
        intermediate_vcfs['VarDict']['indel'] = sorted_indel_vardict_out

    if lofreq:
        snv_lofreq_out   = outdir + os.sep + 'snv.lofreq.vcf'
        indel_lofreq_out = outdir + os.sep + 'indel.lofreq.vcf'

        splitVcf.split_into_snv_and_indel(lofreq, snv_lofreq_out, indel_lofreq_out)

        intermediate_files.extend( [snv_lofreq_out, indel_lofreq_out] )
        snv_intermediates.append(snv_lofreq_out)
        indel_intermediates.append(indel_lofreq_out)
        intermediate_vcfs['LoFreq']['snv']   = snv_varscan_out
        intermediate_vcfs['LoFreq']['indel'] = indel_varscan_out

    if scalpel:
        scalpel_out = outdir + os.sep + 'indel.scalpel.vcf'
        copy_TextFile.copy(scalpel, scalpel_out)
        intermediate_files.append(scalpel_out)
        indel_intermediates.append(scalpel_out)
        
    if strelka:
        import vcfModifier.modify_ssStrelka as mod_strelka
        
        snv_strelka_out   = outdir + os.sep + 'snv.strelka.vcf'
        indel_strelka_out = outdir + os.sep + 'indel.strelka.vcf'
        mod_strelka.convert(strelka, snv_strelka_out, indel_strelka_out)
        
        intermediate_files.extend( [snv_strelka_out, indel_strelka_out] )
        snv_intermediates.append(snv_strelka_out)
        indel_intermediates.append(indel_strelka_out)
        intermediate_vcfs['Strelka']['snv']   = snv_strelka_out
        intermediate_vcfs['Strelka']['indel'] = indel_strelka_out


    # Combine SNV/INDEL variant candidates
    snv_combined   = outdir + os.sep + 'unsorted.CombineVariants.snv.vcf'
    indel_combined = outdir + os.sep + 'unsorted.CombineVariants.indel.vcf'
    
    getUniqueVcfPositions.combine(snv_intermediates, snv_combined)
    getUniqueVcfPositions.combine(indel_intermediates, indel_combined)
    intermediate_files.extend( [snv_combined, indel_combined] )
    
    # Sort them:
    snv_combined_sorted = outdir + os.sep + 'CombineVariants.snv.vcf'
    indel_combined_sorted = outdir + os.sep + 'CombineVariants.indel.vcf'
    
    vcfsorter(hg_dict, snv_combined,   snv_combined_sorted)
    vcfsorter(hg_dict, indel_combined, indel_combined_sorted)
    
    lastSnv = snv_combined_sorted
    lastIndel = indel_combined_sorted
    
    if exclusion or inclusion:
        
        if exclusion:

            newSnv   = outdir + os.sep + 'CombineVariants.snv.ex.vcf'
            newIndel = outdir + os.sep + 'CombineVariants.indel.ex.vcf'
            
            os.system( 'intersectBed -header -a {} -b {} -v | uniq > {}'.format(lastSnv,   exclusion, newSnv) )
            os.system( 'intersectBed -header -a {} -b {} -v | uniq > {}'.format(lastIndel, exclusion, newIndel) )
            
            intermediate_files.extend( [lastSnv, lastIndel] )
            
            lastSnv = newSnv
            lastIndel = newIndel

        if inclusion:
            
            newSnv   = outdir + os.sep + 'CombineVariants.snv.in.vcf'
            newIndel = outdir + os.sep + 'CombineVariants.indel.in.vcf'            
            
            os.system( 'intersectBed -header -a {} -b {} | uniq > {}'.format(lastSnv,   inclusion, newSnv) )
            os.system( 'intersectBed -header -a {} -b {} | uniq > {}'.format(lastIndel, inclusion, newIndel) )
            
            intermediate_files.extend( [lastSnv, lastIndel] )
            
            lastSnv = newSnv
            lastIndel = newIndel


    if not keep_intermediates:
        for file_i in intermediate_files:
            subprocess.call( ('rm', '-v', file_i ) )
    
    return lastSnv, lastIndel, intermediate_vcfs, intermediate_files






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
    if mutect2:
        
        import vcfModifier.modify_MuTect2 as mod_mutect2
        
        mutect2_ex = bed_exclude(mutect2,    exclusion, outdir + os.sep + 'mutect.ex.vcf')
        mutect2_in = bed_include(mutect2_ex, inclusion, outdir + os.sep + 'mutect.in.vcf')
        
        snv_mutect_out   = outdir + os.sep + 'snv.mutect.vcf'
        indel_mutect_out = outdir + os.sep + 'indel.mutect.vcf'
        mod_mutect2.convert(mutect2_in, snv_mutect_out, indel_mutect_out, False)
        
        for file_i in mutect2_ex, mutect2_in, snv_mutect_out, indel_mutect_out
            intermediate_files.add( file_i )
        
        snv_intermediates.append(snv_mutect_out)
        indel_intermediates.append(indel_mutect_out)
        
        intermediate_vcfs['MuTect2']['snv']   = snv_mutect_out
        intermediate_vcfs['MuTect2']['indel'] = indel_mutect_out
    
    if varscan_snv or varscan_indel:
                
        import vcfModifier.modify_VarScan2 as mod_varscan2
        
        if varscan_snv:
            
            varscan_ex = bed_exclude(varscan_snv, exclusion, os.sep + 'snv.varscan.ex.vcf')
            varscan_in = bed_include(varscan_ex,  inclusion, os.sep + 'snv.varscan.in.vcf')
            
            snv_varscan_out = outdir + os.sep + 'snv.varscan.vcf'
            mod_varscan2.convert(varscan_in, snv_varscan_out)
            
            for file_i in varscan_ex, varscan_in, snv_varscan_out:
                intermediate_files.add(file_i)
            
            snv_intermediates.append(snv_varscan_out)
            
        if varscan_indel:

            varscan_ex = bed_exclude(varscan_indel, exclusion, os.sep + 'indel.varscan.ex.vcf')
            varscan_in = bed_include(varscan_ex,    inclusion, os.sep + 'indel.varscan.in.vcf')
            
            indel_varscan_out = outdir + os.sep + 'indel.varscan.vcf'
            mod_varscan2.convert(varscan_in, indel_varscan_out)
            
            for file_i in varscan_ex, varscan_in, indel_varscan_out:
                intermediate_files.add(file_i)
            
            indel_intermediates.append(indel_varscan_out)
    
    if jsm:
        import vcfModifier.modify_JointSNVMix2 as mod_jsm
        
        jsm_ex = bed_exclude(jsm,    exclusion, outdir + os.sep + 'snv.jsm.ex.vcf')
        jsm_in = bed_include(jsm_ex, inclusion, outdir + os.sep + 'snv.jsm.in.vcf')        
        
        jsm_out = outdir + os.sep + 'snv.jsm.vcf'
        mod_jsm.convert(jsm_in, jsm_out)
        
        for file_i in jsm_ex, jsm_in, jsm_out:
            intermediate_files.add(file_i)
            
        snv_intermediates.append(jsm_out)
        
    if sniper:
        import vcfModifier.modify_SomaticSniper as mod_sniper
        
        sniper_ex = bed_exclude(sniper,    exclusion, outdir + os.sep + 'snv.somaticsniper.ex.vcf')
        sniper_in = bed_include(sniper_ex, inclusion, outdir + os.sep + 'snv.somaticsniper.in.vcf')        
        
        sniper_out = outdir + os.sep + 'snv.somaticsniper.vcf'
        mod_sniper.convert(sniper_in, sniper_out)
        
        for file_i in sniper_ex, sniper_in, sniper_out:
            intermediate_files.add(file_i)
        
        snv_intermediates.append(sniper_out)
        
    if vardict:
        import vcfModifier.modify_VarDict as mod_vardict

        vardict_ex = bed_exclude(vardict,    exclusion, outdir + os.sep + 'snv.vardict.ex.vcf')
        vardict_in = bed_include(vardict_ex, inclusion, outdir + os.sep + 'snv.vardict.in.vcf')        
        
        snv_vardict_out   = outdir + os.sep + 'snv.vardict.vcf'
        indel_vardict_out = outdir + os.sep + 'indel.vardict.vcf'
        mod_vardict.convert(vardict_in, snv_vardict_out, indel_vardict_out)
        
        sorted_snv_vardict_out   = outdir + os.sep + 'snv.sort.vardict.vcf'
        sorted_indel_vardict_out = outdir + os.sep + 'indel.sort.vardict.vcf'
        
        vcfsorter(hg_dict, snv_vardict_out,   sorted_snv_vardict_out)
        vcfsorter(hg_dict, indel_vardict_out, sorted_indel_vardict_out)
        
        for file_i in vardict_ex, vardict_in, snv_vardict_out, indel_vardict_out, sorted_snv_vardict_out, sorted_indel_vardict_out
            intermediate_files.add(file_i)
        
        snv_intermediates.append(sorted_snv_vardict_out)
        indel_intermediates.append(sorted_indel_vardict_out)
        intermediate_vcfs['VarDict']['snv']   = sorted_snv_vardict_out
        intermediate_vcfs['VarDict']['indel'] = sorted_indel_vardict_out
        
    if muse:

        muse_ex = bed_exclude(muse,    exclusion, outdir + os.sep + 'snv.muse.ex.vcf')
        muse_in = bed_include(muse_ex, inclusion, outdir + os.sep + 'snv.muse.in.vcf')        
        
        muse_out = outdir + os.sep + 'snv.muse.vcf'
        copy_TextFile.copy(muse_in, muse_out)
        
        for file_i in muse_ex, muse_in, muse_out:
            intermediate_files.add(file_i)
        
        snv_intermediates.append(muse_out)
        
    if lofreq_snv:
        
        lofreq_ex = bed_exclude(lofreq_snv,    exclusion, outdir + os.sep + 'snv.lofreq.ex.vcf')
        lofreq_in = bed_include(lofreq_ex, inclusion, outdir + os.sep + 'snv.lofreq.in.vcf')        

        snv_lofreq_out = outdir + os.sep + 'snv.lofreq.vcf'
        copy_TextFile.copy(lofreq_in, snv_lofreq_out)
        
        for file_i in lofreq_ex, lofreq_in, snv_lofreq_out:
            intermediate_files.add(file_i)
        
        snv_intermediates.append(snv_lofreq_out)

    if lofreq_indel:
        
        lofreq_ex = bed_exclude(lofreq_indel, exclusion, outdir + os.sep + 'indel.lofreq.ex.vcf')
        lofreq_in = bed_include(lofreq_ex,    inclusion, outdir + os.sep + 'indel.lofreq.in.vcf')
        
        indel_lofreq_out = outdir + os.sep + 'indel.lofreq.vcf'
        copy_TextFile.copy(lofreq_in, indel_lofreq_out)
        
        for file_i in lofreq_ex, lofreq_in, indel_lofreq_out:
            intermediate_files.add(file_i)
        
        indel_intermediates.append(indel_lofreq_out)
        
    if scalpel:
        scalpel_out = outdir + os.sep + 'indel.scalpel.vcf'
        copy_TextFile.copy(scalpel, scalpel_out)
        intermediate_files.append(scalpel_out)
        indel_intermediates.append(scalpel_out)
    
    if strelka_snv or strelka_indel:
        
        import vcfModifier.modify_Strelka as mod_strelka

        if strelka_snv:
            snv_strelka_out = outdir + os.sep + 'snv.strelka.vcf'
            mod_strelka.convert(strelka_snv, snv_strelka_out)
            intermediate_files.append(snv_strelka_out)
            snv_intermediates.append(snv_strelka_out)

        if strelka_indel:
            indel_strelka_out = outdir + os.sep + 'indel.strelka.vcf'
            mod_strelka.convert(strelka_indel, indel_strelka_out)
            intermediate_files.append(indel_strelka_out)
            indel_intermediates.append(indel_strelka_out)
            
    if tnscope:

        import vcfModifier.modify_MuTect2 as mod_mutect2
        
        snv_tnscope_out   = outdir + os.sep + 'snv.tnscope.vcf'
        indel_tnscope_out = outdir + os.sep + 'indel.tnscope.vcf'
        mod_mutect2.convert(tnscope, snv_tnscope_out, indel_tnscope_out, True)
        
        intermediate_files.extend( [snv_tnscope_out, indel_tnscope_out] )
        snv_intermediates.append(snv_tnscope_out)
        indel_intermediates.append(indel_tnscope_out)
        intermediate_vcfs['TNscope']['snv']   = snv_tnscope_out
        intermediate_vcfs['TNscope']['indel'] = indel_tnscope_out
    
    
    # Combine SNV/INDEL variant candidates
    snv_combined   = outdir + os.sep + 'unsorted.CombineVariants.snv.vcf'
    indel_combined = outdir + os.sep + 'unsorted.CombineVariants.indel.vcf'
    
    getUniqueVcfPositions.combine(snv_intermediates, snv_combined)
    getUniqueVcfPositions.combine(indel_intermediates, indel_combined)
    intermediate_files.extend( [snv_combined, indel_combined] )
    
    # Sort them:
    snv_combined_sorted = outdir + os.sep + 'CombineVariants.snv.vcf'
    indel_combined_sorted = outdir + os.sep + 'CombineVariants.indel.vcf'
    
    vcfsorter(hg_dict, snv_combined,   snv_combined_sorted)
    vcfsorter(hg_dict, indel_combined, indel_combined_sorted)
    
    lastSnv = snv_combined_sorted
    lastIndel = indel_combined_sorted
    
    if exclusion or inclusion:
        
        if exclusion:

            newSnv   = outdir + os.sep + 'CombineVariants.snv.ex.vcf'
            newIndel = outdir + os.sep + 'CombineVariants.indel.ex.vcf'
            
            os.system( 'intersectBed -header -a {} -b {} -v | uniq > {}'.format(lastSnv,   exclusion, newSnv) )
            os.system( 'intersectBed -header -a {} -b {} -v | uniq > {}'.format(lastIndel, exclusion, newIndel) )
            
            intermediate_files.extend( [lastSnv, lastIndel] )
            
            lastSnv = newSnv
            lastIndel = newIndel

        if inclusion:
            
            newSnv   = outdir + os.sep + 'CombineVariants.snv.in.vcf'
            newIndel = outdir + os.sep + 'CombineVariants.indel.in.vcf'            
            
            os.system( 'intersectBed -header -a {} -b {} | uniq > {}'.format(lastSnv,   inclusion, newSnv) )
            os.system( 'intersectBed -header -a {} -b {} | uniq > {}'.format(lastIndel, inclusion, newIndel) )
            
            intermediate_files.extend( [lastSnv, lastIndel] )
            
            lastSnv = newSnv
            lastIndel = newIndel
        
    if not keep_intermediates:
        for file_i in intermediate_files:
            subprocess.call( ('rm', '-v', file_i ) )
    
    return lastSnv, lastIndel, intermediate_vcfs, intermediate_files

