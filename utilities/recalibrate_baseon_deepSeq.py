#!/usr/bin/env python3

# Recalibrate confidence-level of low VAF calls with high number of PASSES
# Rename StrongEvidence/WeakEvidence/NeutralEvidence/LikelyFalsePositive into HighConf/MedConf/LowConf/Unclassified

import sys, argparse, math, gzip, os
import regex as re
from copy import copy

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomic_file_handlers as genome

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-infile',   '--vcf-infile',   type=str, help='VCF in', required=True)
parser.add_argument('-outfile',  '--vcf-outfile',      type=str, help='VCF out', required=True)
parser.add_argument('-ref',      '--genome-reference',  type=str,   help='.fasta.fai file to get the contigs', required=True)

parser.add_argument('--bignova-bwa',    type=str, help='combined novaseq', )
parser.add_argument('--bignova-bowtie', type=str, help='combined novaseq', )
parser.add_argument('--bignova-novo',   type=str, help='combined novaseq', )

parser.add_argument('--spp-bwa',    type=str, help='300x', )
parser.add_argument('--spp-bowtie', type=str, help='300x', )
parser.add_argument('--spp-novo',   type=str, help='300x', )

args = parser.parse_args()

infile  = args.vcf_infile
outfile = args.vcf_outfile

nova_bwa    = args.bignova_bwa
nova_bowtie = args.bignova_bowtie
nova_novo   = args.bignova_novo
spp_bwa     = args.spp_bwa
spp_bowtie  = args.spp_bowtie
spp_novo    = args.spp_novo


fai_file  = args.genome_reference + '.fai'
chrom_seq = genome.faiordict2contigorder(fai_file, 'fai')


# def rename(line_i):
    
    # line_i = re.sub( 'StrongEvidence',      'HighConf',     line_i )
    # line_i = re.sub( 'WeakEvidence',        'MedConf',      line_i )
    # line_i = re.sub( 'NeutralEvidence',     'LowConf',      line_i )
    # line_i = re.sub( 'LikelyFalsePositive', 'Unclassified', line_i )

    # return line_i



def relabel(vcf_line, newLabel, additional_flag):
    
    vcf_i = genome.Vcf_line( vcf_line )
    item  = vcf_line.split('\t')
    
    filterColumn = re.sub(r'HighConf|MedConf|LowConf|Unclassified', newLabel, vcf_i.filters)
    item         = vcf_line.split('\t')
    item[6]      = filterColumn
    
    originalLabel = re.search(r'(HighConf|MedConf|LowConf|Unclassified)', vcf_i.filters).groups()[0]
    
    if 'FLAGS' in vcf_i.info:
        infoItems = vcf_i.info.split(';')
        for i, item_i in enumerate(infoItems):
            if item_i.startswith('FLAGS'):
                infoItems[i] = infoItems[i] + ',{}'.format(additional_flag)
        newInfo = ';'.join(infoItems)
    else:
        newInfo = vcf_i.info + ';FLAGS={}'.format(additional_flag)
    
    item[7] = newInfo
    line_i  = '\t'.join(item)
    
    return line_i



with genome.open_textfile(infile) as fin,  open(outfile, 'w') as fout:
    
    line_i = fin.readline().rstrip()

    # 450X NovaSeq data
    nova_bwa      = genome.open_textfile(nova_bwa)
    nova_bwa_line = nova_bwa.readline().rstrip()
    while nova_bwa_line.startswith('#'):
        nova_bwa_line = nova_bwa.readline().rstrip()
    
    nova_bowtie = genome.open_textfile(nova_bowtie)
    nova_bowtie_line = nova_bowtie.readline().rstrip()
    while nova_bowtie_line.startswith('#'):
        nova_bowtie_line = nova_bowtie.readline().rstrip()

    nova_novo   = genome.open_textfile(nova_novo)
    nova_novo_line = nova_novo.readline().rstrip()
    while nova_novo_line.startswith('#'):
        nova_novo_line = nova_novo.readline().rstrip()

    # 300X data sets from Genentech
    spp_bwa      = genome.open_textfile(spp_bwa)
    spp_bwa_line = spp_bwa.readline().rstrip()
    while spp_bwa_line.startswith('#'):
        spp_bwa_line = spp_bwa.readline().rstrip()

    spp_bowtie = genome.open_textfile(spp_bowtie)
    spp_bowtie_line = spp_bowtie.readline().rstrip()
    while spp_bowtie_line.startswith('#'):
        spp_bowtie_line = spp_bowtie.readline().rstrip()

    spp_novo = genome.open_textfile(spp_novo)
    spp_novo_line = spp_novo.readline().rstrip()
    while spp_novo_line.startswith('#'):
        spp_novo_line = spp_novo.readline().rstrip()


    # Copy the headlines, but add two additional lines
    while line_i.startswith('#'):

        #line_i = rename( line_i )
        fout.write( line_i + '\n' )
        line_i = fin.readline().rstrip()


    # First coordinate:
    coordinate_i = re.match( genome.pattern_chr_position, line_i )
    coordinate_i = coordinate_i.group() if coordinate_i else ''

    while line_i:
        
        #line_i = rename( line_i )
        my_vcf = genome.Vcf_line( line_i )
        
        my_coordinates = [(my_vcf.chromosome, my_vcf.position)]
        
        variants_at_my_coordinate = []
        variants_at_my_coordinate.append( my_vcf )

        # As long as the "coordinate" stays the same, it will keep reading until it's different.
        while my_coordinates[0] == (my_vcf.chromosome, my_vcf.position):

            line_i = fin.readline().rstrip()
            #line_i = rename( line_i )
            my_vcf = genome.Vcf_line( line_i )


            ###### VCF order checking ######
            coordinate_j = re.match( genome.pattern_chr_position, line_i )
            coordinate_j = coordinate_j.group() if coordinate_j else ''
            if genome.whoisbehind(coordinate_i, coordinate_j, chrom_seq) == 1:
                raise Exception( '{} does not seem to be properly sorted.'.format(infile) )
            coordinate_i = coordinate_j

            if my_coordinates[0] == (my_vcf.chromosome, my_vcf.position):
                variants_at_my_coordinate.append( my_vcf )     
            ###### VCF order checking ######


        ##### ##### ##### ##### ##### #####
        for my_coordinate in my_coordinates:

            ref_bases = []
            alt_bases = []
            
            for variant_i in variants_at_my_coordinate:

                ref_base  = variant_i.refbase
                first_alt = variant_i.altbase.split(',')[0]

                ref_bases.append( ref_base )
                alt_bases.append( first_alt )


            # deepSeq inputs
            got_nova_bwa,    nova_bwa_variants,    nova_bwa_line    = genome.find_vcf_at_coordinate(my_coordinate, nova_bwa_line,    nova_bwa,    chrom_seq)
            got_nova_bowtie, nova_bowtie_variants, nova_bowtie_line = genome.find_vcf_at_coordinate(my_coordinate, nova_bowtie_line, nova_bowtie, chrom_seq)
            got_nova_novo,   nova_novo_variants,   nova_novo_line   = genome.find_vcf_at_coordinate(my_coordinate, nova_novo_line,   nova_novo,   chrom_seq)

            got_spp_bwa,     spp_bwa_variants,     spp_bwa_line     = genome.find_vcf_at_coordinate(my_coordinate, spp_bwa_line,     spp_bwa,     chrom_seq)
            got_spp_bowtie,  spp_bowtie_variants,  spp_bowtie_line  = genome.find_vcf_at_coordinate(my_coordinate, spp_bowtie_line,  spp_bowtie,  chrom_seq)
            got_spp_novo,    spp_novo_variants,    spp_novo_line    = genome.find_vcf_at_coordinate(my_coordinate, spp_novo_line,    spp_novo,    chrom_seq)

            
            # Now, use pysam to look into the BAM file(s), variant by variant from the input:
            for ith_call, my_call in enumerate( variants_at_my_coordinate ):
                
                vcf_items = my_call.vcf_line.split('\t')
                
                tvaf      = float( my_call.get_info_value('TVAF') )
                nPASSES   = int( my_call.get_info_value('nPASSES') )
                nREJECTS  = int( my_call.get_info_value('nREJECTS') )

                bwaVDP,    bwaDP    = [ int(i) for i in my_call.get_info_value('bwaDP').split(',') ]
                bowtieVDP, bowtieDP = [ int(i) for i in my_call.get_info_value('bowtieDP').split(',') ]
                novoVDP,   novoDP   = [ int(i) for i in my_call.get_info_value('novoDP').split(',') ]
            
                VDP = bwaVDP + bowtieVDP + novoVDP
                DP  = bwaDP  + bowtieDP  + novoDP
                
                variant_id = ( (my_call.chromosome, my_call.position), my_call.refbase, my_call.altbase )
                ref_base   = ref_bases[ith_call]
                first_alt  = alt_bases[ith_call]


                # Combined NovaSeq BWA
                if variant_id in nova_bwa_variants:

                    nova_bwa_variant_i = nova_bwa_variants[variant_id]
                    nova_bwa_tvaf = float( nova_bwa_variant_i.get_sample_value('VAF', 1) )
                    
                    if 'PASS' in nova_bwa_variant_i.filters:
                        nova_bwa_PASS    = True
                        nova_bwa_REJECT  = False
                        nova_bwa_Missing = False
                    elif 'REJECT' in nova_bwa_variant_i.filters:
                        nova_bwa_PASS    = False
                        nova_bwa_REJECT  = True
                        nova_bwa_Missing = False
                    else:
                        nova_bwa_PASS    = False
                        nova_bwa_REJECT  = False
                        nova_bwa_Missing = False
                else:
                    nova_bwa_PASS    = False
                    nova_bwa_REJECT  = False
                    nova_bwa_Missing = True


                # Combined NovaSeq Bowtie
                if variant_id in nova_bowtie_variants:

                    nova_bowtie_variant_i = nova_bowtie_variants[variant_id]
                    nova_bowtie_tvaf = float( nova_bowtie_variant_i.get_sample_value('VAF', 1) )
                    
                    if 'PASS' in nova_bowtie_variant_i.filters:
                        nova_bowtie_PASS    = True
                        nova_bowtie_REJECT  = False
                        nova_bowtie_Missing = False
                    elif 'REJECT' in nova_bowtie_variant_i.filters:
                        nova_bowtie_PASS    = False
                        nova_bowtie_REJECT  = True
                        nova_bowtie_Missing = False
                    else:
                        nova_bowtie_PASS    = False
                        nova_bowtie_REJECT  = False
                        nova_bowtie_Missing = False
                else:
                    nova_bowtie_PASS    = False
                    nova_bowtie_REJECT  = False
                    nova_bowtie_Missing = True


                # Combined NovaSeq NovoAlign
                if variant_id in nova_novo_variants:

                    nova_novo_variant_i = nova_novo_variants[variant_id]
                    nova_novo_tvaf = float( nova_novo_variant_i.get_sample_value('VAF', 1) )
                    
                    if 'PASS' in nova_novo_variant_i.filters:
                        nova_novo_PASS    = True
                        nova_novo_REJECT  = False
                        nova_novo_Missing = False
                    elif 'REJECT' in nova_novo_variant_i.filters:
                        nova_novo_PASS    = False
                        nova_novo_REJECT  = True
                        nova_novo_Missing = False
                    else:
                        nova_novo_PASS    = False
                        nova_novo_REJECT  = False
                        nova_novo_Missing = False
                else:
                    nova_novo_PASS    = False
                    nova_novo_REJECT  = False
                    nova_novo_Missing = True



                # Combined SPP BWA
                if variant_id in spp_bwa_variants:

                    spp_bwa_variant_i = spp_bwa_variants[variant_id]
                    spp_bwa_tvaf = float( spp_bwa_variant_i.get_sample_value('VAF', 1) )

                    if 'PASS' in spp_bwa_variant_i.filters:
                        spp_bwa_PASS    = True
                        spp_bwa_REJECT  = False
                        spp_bwa_Missing = False
                    elif 'REJECT' in spp_bwa_variant_i.filters:
                        spp_bwa_PASS    = False
                        spp_bwa_REJECT  = True
                        spp_bwa_Missing = False
                    else:
                        spp_bwa_PASS    = False
                        spp_bwa_REJECT  = False
                        spp_bwa_Missing = False
                else:
                    spp_bwa_PASS    = False
                    spp_bwa_REJECT  = False
                    spp_bwa_Missing = True


                # Combined SPP Bowtie
                if variant_id in spp_bowtie_variants:

                    spp_bowtie_variant_i = spp_bowtie_variants[variant_id]
                    spp_bowtie_tvaf = float( spp_bowtie_variant_i.get_sample_value('VAF', 1) )

                    if 'PASS' in spp_bowtie_variant_i.filters:
                        spp_bowtie_PASS    = True
                        spp_bowtie_REJECT  = False
                        spp_bowtie_Missing = False
                    elif 'REJECT' in spp_bowtie_variant_i.filters:
                        spp_bowtie_PASS    = False
                        spp_bowtie_REJECT  = True
                        spp_bowtie_Missing = False
                    else:
                        spp_bowtie_PASS    = False
                        spp_bowtie_REJECT  = False
                        spp_bowtie_Missing = False
                else:
                    spp_bowtie_PASS    = False
                    spp_bowtie_REJECT  = False
                    spp_bowtie_Missing = True


                # Combined SPP Novo
                if variant_id in spp_novo_variants:

                    spp_novo_variant_i = spp_novo_variants[variant_id]
                    spp_novo_tvaf = float( spp_novo_variant_i.get_sample_value('VAF', 1) )

                    if 'PASS' in spp_novo_variant_i.filters:
                        spp_novo_PASS    = True
                        spp_novo_REJECT  = False
                        spp_novo_Missing = False
                    elif 'REJECT' in spp_novo_variant_i.filters:
                        spp_novo_PASS    = False
                        spp_novo_REJECT  = True
                        spp_novo_Missing = False
                    else:
                        spp_novo_PASS    = False
                        spp_novo_REJECT  = False
                        spp_novo_Missing = False
                else:
                    spp_novo_PASS    = False
                    spp_novo_REJECT  = False
                    spp_novo_Missing = True



                # By NovaSeq or SPP
                nova_hasPASS    = nova_bwa_PASS    or nova_bowtie_PASS    or nova_novo_PASS
                nova_hasREJECT  = nova_bwa_REJECT  or nova_bowtie_REJECT  or nova_novo_REJECT
                nova_hasMissing = nova_bwa_Missing or nova_bowtie_Missing or nova_novo_Missing

                spp_hasPASS     = spp_bwa_PASS    or  spp_bowtie_PASS
                spp_hasREJECT   = spp_bwa_REJECT  or  spp_bowtie_REJECT
                spp_hasMissing  = spp_bwa_Missing or  spp_bowtie_Missing

                # By aligner
                bwa_hasPASS     = nova_bwa_PASS    or spp_bwa_PASS
                bwa_REJECT      = nova_bwa_REJECT  or spp_bwa_REJECT
                bwa_Missing     = nova_bwa_Missing or spp_bwa_Missing

                bowtie_hasPASS  = nova_bowtie_PASS    or spp_bowtie_PASS
                bowtie_REJECT   = nova_bowtie_REJECT  or spp_bowtie_REJECT
                bowtie_Missing  = nova_bowtie_Missing or spp_bowtie_Missing

                novo_hasPASS    = nova_novo_PASS    or spp_novo_PASS
                novo_REJECT     = nova_novo_REJECT  or spp_novo_REJECT
                novo_Missing    = nova_novo_Missing or spp_novo_Missing


                
                if ( ('Unclassified' in my_call.filters) or ('LowConf' in my_call.filters) or ('MedConf' in my_call.filters) ) and \
                (tvaf <= 0.15 and nPASSES >= 15 and nREJECTS <= 10):

                    if ( ('LowConf' in my_call.filters) or ('MedConf' in my_call.filters) ) and nREJECTS == 0 and \
                       ( nova_bwa_PASS + spp_bwa_PASS + nova_bowtie_PASS + spp_bowtie_PASS + nova_novo_PASS + spp_novo_PASS == 6 ) :

                        confLabel_i = re.search(r'Unclassified|LowConf|MedConf', my_call.filters).group()
                        line_out    = relabel(my_call.vcf_line, 'HighConf', '{}_to_{}_by_300X'.format(confLabel_i, 'HighConf'))
                        
                    ##########
                    # Promote to "WeakEvidence"
                    elif ( nova_bwa_PASS + spp_bwa_PASS + nova_bowtie_PASS + spp_bowtie_PASS + nova_novo_PASS + spp_novo_PASS >=4 ) and \
                         ( nova_hasPASS and spp_hasPASS and bwa_hasPASS and bowtie_hasPASS and novo_hasPASS) and \
                         ( not (nova_hasREJECT or nova_hasMissing or spp_hasREJECT or spp_hasMissing) ):

                        confLabel_i = re.search(r'Unclassified|LowConf|MedConf', my_call.filters).group()
                        
                        if confLabel_i != 'MedConf':
                            line_out = relabel(my_call.vcf_line, 'MedConf', '{}_to_{}_by_300X'.format(confLabel_i, 'MedConf'))
                        else:
                            line_out = my_call.vcf_line

                    # Promote one ladder up if PASS by Burrows-Wheeler (bwa or bowtie) and NovoAlign
                    # "Missing" in 450X is worse than "missing" in 50X, so it's considered as "bad" as REJECT, i.e., two PASSES and one LowQual
                    elif ( nova_hasPASS and spp_hasPASS and bwa_hasPASS and bowtie_hasPASS and novo_hasPASS) and \
                         ( not (nova_hasREJECT or nova_hasMissing or spp_hasREJECT or spp_hasMissing) ):

                        confLabel_i = re.search(r'Unclassified|LowConf|MedConf', my_call.filters).group()
                        
                        if confLabel_i == 'Unclassified':
                            line_out = relabel(line_i, 'LowConf', '{}_to_{}_by_300X'.format(confLabel_i, 'LowConf'))
                        elif confLabel_i == 'LowConf':
                            line_out = relabel(my_call.vcf_line, 'MedConf', '{}_to_{}_by_300X'.format(confLabel_i, 'MedConf'))
                        else:
                            line_out = my_call.vcf_line

                    elif nova_hasPASS and spp_hasPASS and bwa_hasPASS and bowtie_hasPASS and novo_hasPASS:

                        confLabel_i = re.search(r'Unclassified|LowConf|MedConf', my_call.filters).group()

                        if confLabel_i == 'Unclassified':
                            line_out = relabel(my_call.vcf_line, 'LowConf', '{}_to_{}_by_300X'.format(confLabel_i, 'LowConf'))
                        else:
                            line_out = my_call.vcf_line

                    # Demotion
                    elif (not (nova_hasPASS or spp_hasPASS)) and (nova_hasREJECT or nova_hasMissing) and (spp_hasREJECT or spp_hasMissing):

                        confLabel_i = re.search(r'Unclassified|LowConf|MedConf', my_call.filters).group()

                        if confLabel_i == 'MedConf':
                            line_out = relabel(my_call.vcf_line, 'LowConf', '{}_to_{}_by_300X'.format(confLabel_i, 'LowConf'))
                        # elif confLabel_i == 'LowConf':
                            # line_out = relabel(line_i, 'Unclassified', '{}_to_{}_by_300X'.format(confLabel_i, 'Unclassified'))
                        else:
                            line_out = my_call.vcf_line
                            
                    else:
                        line_out = my_call.vcf_line


                elif 'LowConf' in my_call.filters  and tvaf >= 0.3:
                    
                    if (nova_hasREJECT or nova_hasMissing) and (spp_hasREJECT or spp_hasMissing) and (bwa_REJECT or bwa_Missing) and (bowtie_REJECT or bowtie_Missing) and (novo_REJECT or novo_Missing) and ( not (nova_hasPASS or spp_hasPASS) ):

                        confLabel_i = re.search(r'Unclassified|LowConf|MedConf', my_call.filters).group()

                        if confLabel_i == 'LowConf':
                            line_out = relabel(my_call.vcf_line, 'Unclassified', '{}_to_{}_by_300X'.format(confLabel_i, 'Unclassified'))
                        else:
                            line_out = my_call.vcf_line
                            
                    else:
                        line_out = my_call.vcf_line


                else:
                    line_out = my_call.vcf_line

                # Write
                fout.write( line_out + '\n' )


    opened_files = (nova_bwa, nova_bowtie, nova_novo)
    [ opened_file.close() for opened_file in opened_files if opened_file ]
