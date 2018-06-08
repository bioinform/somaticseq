#!/usr/bin/env python3

import sys, argparse, gzip, re, os

MY_DIR = os.path.dirname(os.path.realpath(__file__))
PRE_DIR = os.path.join(MY_DIR, os.pardir)
sys.path.append( PRE_DIR )

import genomic_file_handlers as genome
from read_info_extractor import * 

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-vcfin',    '--vcf-infile',  type=str, help='VCF in', required=True)
parser.add_argument('-filters',  '--vcf-filters', type=str, nargs='*', help='Does not actually work, hard coded actually', required=False, default=['StrongEvidence', 'WeakEvidence'])
parser.add_argument('-sample',    '--sample-out',  type=str, help='sample out', required=True)
parser.add_argument('-total',     '--total-out',   type=str, help='total out', required=True)


args = parser.parse_args()

vcf_in      = args.vcf_infile
vcf_filters = args.vcf_filters
sampleout   = args.sample_out
totalout    = args.total_out

increment = 100
nan = float('nan')

with genome.open_textfile(vcf_in) as vcfin:
    
    line_i = vcfin.readline().rstrip()
    
    while line_i.startswith('##'):
        line_i = vcfin.readline().rstrip()

    
    # This is the #CHROM line:
    header  = line_i.split('\t')
    samples = header[9::]
    
    bwa_tumors     = []
    bowtie_tumors  = []
    novo_tumors    = []
    
    for tSample_i in samples:
        if  tSample_i.endswith('.bwa'):
            bwa_tumors.append(tSample_i)
        elif tSample_i.endswith('.bowtie'):
            bowtie_tumors.append(tSample_i)
        elif tSample_i.endswith('.novo'):
            novo_tumors.append(tSample_i)
    
    tumor_indices = {}
    tumor_indices['bwa']     = [ samples.index(tSample_i) for tSample_i in bwa_tumors     ]
    tumor_indices['bowtie']  = [ samples.index(tSample_i) for tSample_i in bowtie_tumors  ]
    tumor_indices['novo']    = [ samples.index(tSample_i) for tSample_i in novo_tumors    ]
    

    line_i = vcfin.readline().rstrip()
    
    
    dp2vaf   = { 'bwa': {}, 'bowtie': {}, 'novo': {}, }
    dp2altdp = { 'bwa': {}, 'bowtie': {}, 'novo': {}, }
    
    totaldp2vaf   = { 'bwa': {}, 'bowtie': {}, 'novo': {}, }
    totaldp2altdp = { 'bwa': {}, 'bowtie': {}, 'novo': {}, }
    
    while line_i:
        
        my_vcf = genome.Vcf_line( line_i )
        
        
        if re.search(r'(Strong|Weak)Evidence', my_vcf.filters):
            
            for aligner_i in ('bwa', 'bowtie', 'novo'):
                
                total_altdp, total_dp = my_vcf.get_info_value( '{}DP'.format(aligner_i) ).split(',')
                total_altdp, total_dp = int(total_altdp), int(total_dp)
                
                overall_vaf = float(my_vcf.get_info_value( '{}TVAF'.format(aligner_i) ))
                
                if total_dp not in totaldp2altdp:
                    totaldp2altdp[ aligner_i ][ total_dp ] = set()
                    totaldp2vaf[ aligner_i ][ total_dp ] = set()
                    
                totaldp2altdp[ aligner_i ][ total_dp ].add(total_altdp)
                totaldp2vaf[ aligner_i ][ total_dp ].add(overall_vaf)
            
                for sample_i in tumor_indices[aligner_i]:
                    
                    somaticseqScore = my_vcf.get_sample_value('SCORE', sample_i)
                    
                    try:
                        somaticseqScore = float(somaticseqScore)
                        if somaticseqScore > genome.p2phred(0.3) and '1' in my_vcf.get_sample_value('GT', sample_i):
                            
                            try:
                                ref_1, ref_2, alt_1, alt_2 = my_vcf.get_sample_value('DP4', sample_i).split(',')
                                ref_1, ref_2, alt_1, alt_2 = int(ref_1), int(ref_2), int(alt_1), int(alt_2)
                                
                                dp = ref_1 + ref_2 + alt_1 + alt_2
                                alt_dp = alt_1 + alt_2
                                
                                vaf = float(my_vcf.get_sample_value('VAF', sample_i))
                                
                                if dp not in dp2vaf[aligner_i]:
                                    dp2vaf[aligner_i][dp] = set()
                                    dp2altdp[aligner_i][dp] = set()
                                
                                dp2vaf[aligner_i][dp].add( vaf )
                                dp2altdp[aligner_i][dp].add( alt_dp )
                                
                            except ValueError:
                                pass
                            
                    except ValueError:
                        pass
                    
        
        line_i = vcfin.readline().rstrip()



with open(sampleout, 'w') as sampleout:
    
    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format('DP', 'bwa.minAltDP', 'bowtie.minAltDP', 'novo.minAltDP', 'MinAltDP', 'bwa.minVAF', 'bowtie.minVAF', 'novo.minVAF', 'MinVAF'), file=sampleout )


    for dp_i in range(10,150):
        
        bwaMinDP = min( dp2altdp['bwa'][dp_i] )
        bowtieMinDP = min( dp2altdp['bowtie'][dp_i] )
        novoMinDP = min( dp2altdp['novo'][dp_i] )
        MinDP = min(bwaMinDP, bowtieMinDP, novoMinDP)
        
        bwaMinVAF = min( dp2vaf['bwa'][dp_i] )
        bowtieMinVAF = min( dp2vaf['bowtie'][dp_i] )
        novoMinVAF = min( dp2vaf['novo'][dp_i] )
        MinVAF = min( bwaMinVAF, bowtieMinVAF, novoMinVAF)
        
        print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(dp_i, bwaMinDP, bowtieMinDP, novoMinDP, MinDP, bwaMinVAF, bowtieMinVAF, novoMinVAF, MinVAF), file=sampleout )
        




with open(totalout, 'w') as totalout:
    
    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format('DP', 'bwa.minAltDP', 'bowtie.minAltDP', 'novo.minAltDP', 'MinAltDP', 'bwa.minVAF', 'bowtie.minVAF', 'novo.minVAF', 'MinVAF'), file=totalout )


    for dp_i in range(100,5001,increment):
        
        
        bwaMinDPs = set()
        bowtieMinDPs = set()
        novoMinDPs = set ()
        
        bwaMinVAFs = set()
        bowtieMinVAFs = set()
        novoMinVAFs = set ()
        
        
        for dp_j in range(dp_i, dp_i+increment):
        
            if dp_j in totaldp2altdp['bwa']:
                bwaMinDPs.add( min( totaldp2altdp['bwa'][dp_j] ) )
            
            if dp_j in totaldp2altdp['bowtie']:
                bowtieMinDPs.add( min( totaldp2altdp['bowtie'][dp_j] ) )
            
            if dp_j in totaldp2altdp['novo']:
                novoMinDPs.add( min( totaldp2altdp['novo'][dp_j] ) )
            
            
            if dp_j in totaldp2vaf['bwa']:
                bwaMinVAFs.add( min( totaldp2vaf['bwa'][dp_j] ) )
            
            if dp_j in totaldp2altdp['bowtie']:
                bowtieMinVAFs.add( min( totaldp2vaf['bowtie'][dp_j] ) )
            
            if dp_j in totaldp2altdp['novo']:
                novoMinVAFs.add( min( totaldp2vaf['novo'][dp_j] ) )
            
            
        
        bwaMinDP = min( bwaMinDPs ) if bwaMinDPs else nan
        bowtieMinDP = min( bowtieMinDPs ) if bowtieMinDPs else nan
        novoMinDP = min( novoMinDPs ) if novoMinDPs else nan
        
        MinDPs = set.union( bwaMinDPs, bowtieMinDPs , novoMinDPs)
        MinDP = min(MinDPs) if MinDPs else nan
        
        
        bwaMinVAF = min( bwaMinVAFs ) if bwaMinVAFs else nan
        bowtieMinVAF = min( bowtieMinVAFs ) if bowtieMinVAFs else nan
        novoMinVAF = min( novoMinVAFs ) if novoMinVAFs else nan
        
        MinVAFs = set.union( bwaMinVAFs, bowtieMinVAFs , novoMinVAFs)
        MinVAF = min( MinVAFs ) if MinVAFs else nan
        
        print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(dp_i, bwaMinDP, bowtieMinDP, novoMinDP, MinDP, bwaMinVAF, bowtieMinVAF, novoMinVAF, MinVAF), file=totalout )
        
        

