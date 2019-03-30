#!/usr/bin/env python3

import sys, os, gzip, pysam

MY_DIR       = os.path.dirname(os.path.realpath(__file__))

# wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
refGene_file = MY_DIR + os.sep + 'data' + os.sep + 'hg38.refGene.txt.gz'

########################################################################  
"""
http://ucscbrowser.genenetwork.org/FAQ/FAQformat.html
The refGene table is an example of the genePredExt format.
Chromosome in refGene.txt uses 0-based indexing. 
    1)  string name;                "Name of gene (usually transcript_id from GTF)"
    2)  string chrom;               "Chromosome name"
    3)  char[1] strand;             "+ or - for strand"
    4)  uint txStart;               "Transcription start position"
    5)  uint txEnd;                 "Transcription end position"
    6)  uint cdsStart;              "Coding region start"
    7)  uint cdsEnd;                "Coding region end"
    8)  uint exonCount;             "Number of exons"
    9)  uint[exonCount] exonStarts; "Exon start positions"
    10) uint[exonCount] exonEnds;   "Exon end positions"
    11) uint id;                    "Unique identifier"
    12) string name2;               "Alternate name (e.g. gene_id from GTF)"
    13) string cdsStartStat;        "enum('none','unk','incmpl','cmpl')"
    14) string cdsEndStat;          "enum('none','unk','incmpl','cmpl')"
    15) lstring exonFrames;         "Exon frame offsets {0,1,2}"
"""

iStringName   =  1
iChrom        =  2
iStrand       =  3
iStart        =  4
iEnd          =  5
iCdsStart     =  6
iCdsEnd       =  7
iExonCount    =  8
iExonStarts   =  9
iExonEnds     = 10
iIdentifier   = 11
iGene         = 12
iCdsStartStat = 13
iCdsEndStat   = 14
iExonFrames   = 15


refGene_data = {}
with gzip.open(refGene_file, 'rt') as genedata:
    
    line_i = genedata.readline().rstrip('\n')
    
    while line_i:

        item = line_i.split('\t')
        
        gene_i = item[ iGene ]
        
        if gene_i not in refGene_data:
            refGene_data[ gene_i ] = []

        refGene_data[ gene_i ].append( item )
        
        line_i = genedata.readline().rstrip('\n')


class Gene:
    
    def __init__(self, refGene_gene):
        
        self.refGene_gene = refGene_gene
        
        self.transcripts = []
        for item in refGene_gene:
            
            start_position = int( item[iStart] )
            end_position   = int( item[iEnd] )
            cds_start      = int( item[iCdsStart] )
            cds_end        = int( item[iCdsEnd] )
            exon_count     = int( item[iExonCount] )
            exon_starts    = [ int(i) for i in item[iExonStarts].rstrip(',').split(',') ]
            exon_ends      = [ int(i) for i in item[iExonEnds].rstrip(',').split(',') ]
            exon_frames    = [ int(i) for i in item[iExonFrames].rstrip(',').split(',') ]
            
            transcript_i = {'name':           item[iStringName], \
                            'chrom':          item[iChrom], \
                            'strand':         item[iStrand], \
                            'start':          start_position, \
                            'end':            end_position, \
                            'cds_start':      cds_start, \
                            'cds_end':        cds_end, \
                            'exon_count':     exon_count, \
                            'exon_starts':    exon_starts, \
                            'exon_ends':      exon_ends, \
                            'identifier':     item[iIdentifier], \
                            'gene':           item[iGene], \
                            'cds_start_stat': item[iCdsStartStat], \
                            'cds_end_stat':   item[iCdsEndStat], \
                            'exon_frame':     exon_frames}
            
            self.transcripts.append( transcript_i )
