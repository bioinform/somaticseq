#!/usr/bin/env python3

import regex as re
import math
import sys, os, gzip
import genomic_file_handlers as genome

nan = float('nan')
inf = float('inf')


def seq(reads):
    
    '''A function to convert a string into a character generator.'''
    for read_i in reads:
        yield read_i


class Pileup_line:
    
    def __init__(self, pileup_line):
        
        '''Argument is a line in pileup file.'''
        self.pileup_line = pileup_line.rstrip('\n')
        
        try:
            self.chromosome, self.position, self.refbase, self.dp, self.reads, self.qualities = pileup_line.rstrip('\n').split('\t')
            self.position = int(self.position)
            self.dp = int(self.dp)
        
        except ValueError:
            
            try:
                self.chromosome, self.position, self.refbase, self.dp = pileup_line.rstrip('\n').split('\t')
                self.position = int(self.position)
                self.dp = int(self.dp)
                self.reads = self.qualities = ''
            
            except ValueError:
                self.chromosome = self.refbase = self.reads = self.qualities = ''
                self.position = None
                self.dp = 0
    
        
    def alt_read_count(self, pattern):
        alt_reads = re.findall(pattern, self.reads, re.I)
        total_count = len(alt_reads)
        return total_count
        
    def base_reads(self):
        
        # Make the base_calls into a generator:        
        base_calls = seq(self.reads)
        
        ref_forward_count = 0
        ref_reverse_count = 0
        alt_forward       = []
        alt_reverse       = []
        del_forward       = []
        del_reverse       = []
        ins_forward       = []
        ins_reverse       = []
        n_count           = 0
        N_count           = 0
        
        for base_i in base_calls:
            
            # Indicate beginning of a read:
            if base_i == '^':
                
                # Skip the MQ ascii:
                next(base_calls)
            
            # Reference forward reads
            elif base_i == '.':
                ref_forward_count = ref_forward_count + 1
            
            # Reference reverse reads
            elif base_i == ',':
                ref_reverse_count = ref_reverse_count + 1
            
            # Random alt base
            elif base_i == 'n':
                n_count = n_count + 1
            elif base_i == 'N':
                N_count = N_count + 1
            
            # SNP
            elif re.match(r'[GCTAU]', base_i):
                alt_forward.append(base_i)
            elif re.match(r'[gctau]', base_i):
                alt_reverse.append(base_i)
            
            # Deletion
            elif base_i == '-':
                
                num = ''
                base_i = next(base_calls)
                
                while base_i.isnumeric():
                    num = num + base_i
                    base_i = next(base_calls)
                
                num = int(num)
                deleted_seq = base_i
                for i in range(num-1):
                    deleted_seq = deleted_seq + next(base_calls)
                    
                if deleted_seq.isupper():
                    del_forward.append(deleted_seq)
                elif deleted_seq.islower():
                    del_reverse.append(deleted_seq)
            
            # Insertion
            elif base_i == '+':
                num = ''
                base_i = next(base_calls)
                
                while base_i.isnumeric():
                    num = num + base_i
                    base_i = next(base_calls)
                
                num = int(num)
                inserted_seq = base_i
                for i in range(num-1):
                    inserted_seq = inserted_seq + next(base_calls)
                    
                if inserted_seq.isupper():
                    ins_forward.append(inserted_seq)
                elif inserted_seq.islower():
                    ins_reverse.append(inserted_seq)
                    
        return ref_forward_count, ref_reverse_count, alt_forward, alt_reverse, del_forward, del_reverse, ins_forward, ins_reverse, n_count, N_count
                
    
    def total_insertion_calls(self):
        ins = re.findall(r'\+[0-9]+[ACGTNacgtn]+', self.reads)
        total_count = len(ins)
        return total_count
    
    def total_deletion_calls(self):
        dels = re.findall(r'-[0-9]+[ACGTNacgtn]+', self.reads)
        total_count = len(dels)
        return total_count
    
    def indel_fraction(self):
        ins = self.total_insertion_calls()
        dels = self.total_deletion_calls()
        
        try:
            fraction = (ins + dels ) / self.dp
        except ZeroDivisionError:
            fraction = 0
        return fraction
