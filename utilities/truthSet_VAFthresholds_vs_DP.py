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


args = parser.parse_args()

vcf_in      = args.vcf_infile
vcf_filters = args.vcf_filters


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
    
    
    dp2vaf = { 'bwa': {}, 'bowtie': {}, 'novo': {}, }
    dp2altdp = { 'bwa': {}, 'bowtie': {}, 'novo': {}, }
    
    while line_i:
        
        my_vcf = genome.Vcf_line( line_i )
        
        
        if re.search(r'(Strong|Weak)Evidence', my_vcf.filters):
            
            
            for aligner_i in ('bwa', 'bowtie', 'novo'):
            
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


    print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format('DP', 'bwa.minAltDP', 'bowtie.minAltDP', 'novo.minAltDP', 'MinAltDP', 'bwa.minVAF', 'bowtie.minVAF', 'novo.minVAF', 'MinVAF') )


    for dp_i in range(10,150):
        
        bwaMinDP = min( dp2altdp['bwa'][dp_i] )
        bowtieMinDP = min( dp2altdp['bowtie'][dp_i] )
        novoMinDP = min( dp2altdp['novo'][dp_i] )
        MinDP = min(bwaMinDP, bowtieMinDP, novoMinDP)
        
        bwaMinVAF = min( dp2vaf['bwa'][dp_i] )
        bowtieMinVAF = min( dp2vaf['bowtie'][dp_i] )
        novoMinVAF = min( dp2vaf['novo'][dp_i] )
        MinVAF = min( bwaMinVAF, bowtieMinVAF, novoMinVAF)
        
        print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(dp_i, bwaMinDP, bowtieMinDP, novoMinDP, MinDP, bwaMinVAF, bowtieMinVAF, novoMinVAF, MinVAF) )
        
        
