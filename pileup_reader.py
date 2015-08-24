#!/usr/bin/env python3

import regex as re
import math
import sys, os, gzip

nan = float('nan')
inf = float('inf')


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
                self.position = nan
                self.dp = 0
    
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
        
    def alt_read_count(self, pattern):
        alt_reads = re.findall(pattern, self.reads, re.I)
        total_count = len(alt_reads)
        return total_count
        
    
