import re

import somaticseq.genomicFileHandler.genomic_file_handlers as genome
from somaticseq.genomicFileHandler.read_info_extractor import (
    find_msi,
    find_msilen,
    find_shift3,
    mutect2_ecnt,
    mutect2_nlod,
    mutect2_str,
    mutect2_tlod,
)

nan = float("nan")

# Normal/Tumor index in the Merged VCF file, or any other VCF file that puts
# NORMAL first.
idxN, idxT = 0, 1

# Normal/Tumor index in VarDict VCF, or any other VCF file that puts TUMOR
# first.
vdT, vdN = 0, 1


"""
caller_variants is a dictionary, where the key is a tuple of ( (contig,
position), ref, alt ), and value is a genome.VcfLine object. 
"""


def countPASS(variant_id, generic_variants):
    # Most generic classification: 1 if PASS in FILTER, 0 everything else
    if variant_id in generic_variants:
        variant_i = generic_variants[variant_id]
        variant_classification = 1 if re.search(r"\bPASS\b", variant_i.filters) else 0
    else:
        variant_classification = 0

    return variant_classification


def countSOMATICPASS(variant_id, generic_variants):

    # 1 if PASS in FILTER and SOMATIC in INFO, 0 everything else
    if variant_id in generic_variants:
        variant_i = generic_variants[variant_id]
        variant_classification = (
            1
            if (
                re.search(r"\bPASS\b", variant_i.filters)
                and variant_i.get_info_value("SOMATIC")
            )
            else 0
        )
    else:
        variant_classification = 0

    return variant_classification


def MuTect(variant_id, mutect_variants):
    if variant_id in mutect_variants:
        mutect_variant_i = mutect_variants[variant_id]
        mutect_classification = (
            1
            if (
                mutect_variant_i.get_info_value("SOMATIC")
                or "PASS" in mutect_variant_i.filters
            )
            else 0
        )
        # MuTect2 has some useful information:
        nlod = mutect2_nlod(mutect_variant_i)
        tlod = mutect2_tlod(mutect_variant_i)
        tandem = mutect2_str(mutect_variant_i)
        ecnt = mutect2_ecnt(mutect_variant_i)
    else:
        # Not called by mutect
        mutect_classification = 0
        nlod = tlod = tandem = ecnt = nan

    return mutect_classification, nlod, tlod, tandem, ecnt


def ssMuTect(variant_id, mutect_variants):

    if variant_id in mutect_variants:
        mutect_variant_i = mutect_variants[variant_id]
        mutect_classification = 1 if mutect_variant_i.filters == "PASS" else 0
        tlod = mutect2_tlod(mutect_variant_i)
        ecnt = mutect2_ecnt(mutect_variant_i)
    else:
        # Not called by mutect
        mutect_classification = 0
        tlod = ecnt = nan

    return mutect_classification, tlod, ecnt


def VarScan(variant_id, varscan_variants):

    if variant_id in varscan_variants:
        varscan_variant_i = varscan_variants[variant_id]
        varscan_classification = 1 if varscan_variant_i.get_info_value("SOMATIC") else 0
    else:
        varscan_classification = 0

    return varscan_classification


def ssVarScan(variant_id, varscan_variants):
    if variant_id in varscan_variants:
        varscan_variant_i = varscan_variants[variant_id]
        varscan_classification = 1 if varscan_variant_i.filters == "PASS" else 0
        score_varscan2 = eval(varscan_variant_i.get_sample_value("PVAL"))
    else:
        varscan_classification = 0
        score_varscan2 = nan

    return varscan_classification, score_varscan2


def JSM(variant_id, jsm_variants):

    if variant_id in jsm_variants:
        jsm_variant_i = jsm_variants[variant_id]
        jointsnvmix2_classification = 1
        aaab = float(jsm_variant_i.get_info_value("AAAB"))
        aabb = float(jsm_variant_i.get_info_value("AABB"))
        jointsnvmix2_p = 1 - aaab - aabb
        score_jointsnvmix2 = genome.p2phred(jointsnvmix2_p, max_phred=50)

    else:
        jointsnvmix2_classification = 0
        score_jointsnvmix2 = nan

    return jointsnvmix2_classification, score_jointsnvmix2


def SomaticSniper(variant_id, sniper_variants):
    if variant_id in sniper_variants:
        sniper_variant_i = sniper_variants[variant_id]
        sniper_classification = (
            1 if sniper_variant_i.get_sample_value("SS", idxT) == "2" else 0
        )
        if sniper_classification == 1:
            score_somaticsniper = sniper_variant_i.get_sample_value("SSC", idxT)
            score_somaticsniper = (
                int(score_somaticsniper) if score_somaticsniper else nan
            )
        else:
            score_somaticsniper = nan

    else:
        sniper_classification = 0
        score_somaticsniper = nan

    return sniper_classification, score_somaticsniper


def VarDict(variant_id, vardict_variants):
    if variant_id in vardict_variants:
        vardict_variant_i = vardict_variants[variant_id]
        if (vardict_variant_i.filters == "PASS") and (
            "Somatic" in vardict_variant_i.info
        ):
            vardict_classification = 1
        elif "Somatic" in vardict_variant_i.info:
            vardict_filters = vardict_variant_i.filters.split(";")
            disqualifying_filters = (
                ("d7" in vardict_filters or "d5" in vardict_filters)
                or ("DIFF0.2" in vardict_filters)
                or ("LongAT" in vardict_filters)
                or ("MAF0.05" in vardict_filters)
                or ("MSI6" in vardict_filters)
                or ("NM4" in vardict_filters or "NM4.25" in vardict_filters)
                or ("pSTD" in vardict_filters)
                or ("SN1.5" in vardict_filters)
                or (
                    "P0.05" in vardict_filters
                    and float(vardict_variant_i.get_info_value("SSF")) >= 0.15
                )
                or (
                    ("v3" in vardict_filters or "v4" in vardict_filters)
                    and int(vardict_variant_i.get_sample_value("VD", 0)) < 3
                )
            )
            no_bad_filter = not disqualifying_filters
            filter_fail_times = len(vardict_filters)
            if no_bad_filter and filter_fail_times <= 2:
                vardict_classification = 0.5
            else:
                vardict_classification = 0
        else:
            vardict_classification = 0

        # Somatic Score:
        score_vardict = vardict_variant_i.get_info_value("SSF")
        if score_vardict:
            score_vardict = float(score_vardict)
            score_vardict = genome.p2phred(score_vardict, max_phred=100)
        else:
            score_vardict = nan

        # MSI, MSILEN, and SHIFT3:
        msi = find_msi(vardict_variant_i)
        msilen = find_msilen(vardict_variant_i)
        shift3 = find_shift3(vardict_variant_i)
    else:
        vardict_classification = 0
        msi = msilen = shift3 = score_vardict = nan

    return vardict_classification, msi, msilen, shift3, score_vardict


def ssVarDict(variant_id, vardict_variants):

    if variant_id in vardict_variants:
        vardict_variant_i = vardict_variants[variant_id]
        vardict_classification = 1 if vardict_variant_i.filters == "PASS" else 0

        # VarDict reported metrics:
        msi = vardict_variant_i.get_info_value("MSI")
        msi = msi if msi else nan
        msilen = vardict_variant_i.get_info_value("MSILEN")
        msilen = msilen if msilen else nan
        shift3 = vardict_variant_i.get_info_value("SHIFT3")
        shift3 = shift3 if shift3 else nan
        t_pmean = vardict_variant_i.get_info_value("PMEAN")
        t_pmean = t_pmean if t_pmean else nan
        t_pstd = vardict_variant_i.get_info_value("PSTD")
        t_pstd = t_pstd if t_pstd else nan
        t_qstd = vardict_variant_i.get_info_value("QSTD")
        t_qstd = t_qstd if t_qstd else nan
    else:
        # Not called by VarDict
        vardict_classification = 0
        msi = msilen = shift3 = t_pmean = t_pstd = t_qstd = nan

    return vardict_classification, msi, msilen, shift3, t_pmean, t_pstd, t_qstd


def MuSE(variant_id, muse_variants):
    if variant_id in muse_variants:
        muse_variant_i = muse_variants[variant_id]
        if muse_variant_i.filters == "PASS":
            muse_classification = 1
        elif muse_variant_i.filters == "Tier1":
            muse_classification = 0.9
        elif muse_variant_i.filters == "Tier2":
            muse_classification = 0.8
        elif muse_variant_i.filters == "Tier3":
            muse_classification = 0.7
        elif muse_variant_i.filters == "Tier4":
            muse_classification = 0.6
        elif muse_variant_i.filters == "Tier5":
            muse_classification = 0.5
        else:
            muse_classification = 0
    else:
        muse_classification = 0
    return muse_classification


def LoFreq(variant_id, lofreq_variants):

    if variant_id in lofreq_variants:

        lofreq_variant_i = lofreq_variants[variant_id]
        lofreq_classification = 1 if lofreq_variant_i.filters == "PASS" else 0

    else:
        lofreq_classification = 0

    return lofreq_classification


def ssLoFreq(variant_id, lofreq_variants):
    if variant_id in lofreq_variants:
        lofreq_variant_i = lofreq_variants[variant_id]
        lofreq_classification = 1 if lofreq_variant_i.filters == "PASS" else 0
    else:
        lofreq_classification = 0

    return lofreq_classification


def Scalpel(variant_id, scalpel_variants):
    if variant_id in scalpel_variants:
        scalpel_variant_i = scalpel_variants[variant_id]
        if scalpel_variant_i.get_info_value("SOMATIC"):
            if scalpel_variant_i.filters == "PASS":
                scalpel_classification = 1
            else:
                scalpel_classification = 0.5
        else:
            scalpel_classification = 0
    else:
        scalpel_classification = 0

    return scalpel_classification


def ssScalpel(variant_id, scalpel_variants):
    if variant_id in scalpel_variants:
        scalpel_variant_i = scalpel_variants[variant_id]
        scalpel_classification = 1 if scalpel_variant_i.filters == "PASS" else 0
    else:
        scalpel_classification = 0
    return scalpel_classification


def Strelka(variant_id, strelka_variants):
    if variant_id in strelka_variants:
        strelka_variant_i = strelka_variants[variant_id]
        strelka_classification = 1 if "PASS" in strelka_variant_i.filters else 0
        somatic_evs = strelka_variant_i.get_info_value("SomaticEVS")
        if somatic_evs == False:
            somatic_evs = nan
        qss = strelka_variant_i.get_info_value("QSS")
        tqss = strelka_variant_i.get_info_value("TQSS")
    else:
        strelka_classification = 0
        somatic_evs = qss = tqss = nan

    return strelka_classification, somatic_evs, qss, tqss


def ssStrelka(variant_id, strelka_variants):
    if variant_id in strelka_variants:
        strelka_variant_i = strelka_variants[variant_id]
        strelka_classification = 1 if "PASS" in strelka_variant_i.filters else 0
    else:
        strelka_classification = 0
    return strelka_classification


def TNscope(variant_id, tnscope_variants):
    if variant_id in tnscope_variants:
        tnscope_variant_i = tnscope_variants[variant_id]
        tnscope_classification = (
            1
            if (
                tnscope_variant_i.get_info_value("SOMATIC")
                or "PASS" in tnscope_variant_i.filters
            )
            else 0
        )
    else:
        # Not called by TNscope
        tnscope_classification = 0
    return tnscope_classification


def dbSNP(variant_id, dbsnp_variants):
    if variant_id in dbsnp_variants:
        dbsnp_variant_i = dbsnp_variants[variant_id]
        if_dbsnp = 1
        if_common = 1 if dbsnp_variant_i.get_info_value("COMMON") == "1" else 0
        rsID = dbsnp_variant_i.identifier.split(",")
    else:
        if_dbsnp = if_common = 0
        rsID = []
    return if_dbsnp, if_common, rsID


def anyInputVcf(variant_id, input_variants):
    if variant_id in input_variants:
        input_variant_i = input_variants[variant_id]
        if "REJECT" in input_variant_i.filters:
            classification = 0
        elif "LowQual" in input_variant_i.filters:
            classification = 0.5
        else:
            classification = 1
    else:
        classification = 0
    return classification


def COSMIC(variant_id, cosmic_variants):
    if variant_id in cosmic_variants:
        cosmic_variant_i = cosmic_variants[variant_id]
        # If designated as SNP, make it "non-cosmic" and make CNT=nan.
        if cosmic_variant_i.get_info_value("SNP"):
            if_cosmic = 0
            num_cases = nan
        else:
            if_cosmic = 1
            num_cases = cosmic_variant_i.get_info_value("CNT")
            num_cases = num_cases if num_cases else nan
        # COSMIC ID still intact:
        cosmicID = cosmic_variant_i.identifier.split(",")
    else:
        if_cosmic = num_cases = 0
        cosmicID = []
    return if_cosmic, num_cases, cosmicID
