#!/usr/bin/env python3

import sys, re


for line_i in sys.stdin:
    
    if line_i.startswith('#'):
        print(line_i, end='')
        
    else:
        
        item     = line_i.rstrip().split('\t')
        item[3]  = re.sub(r'[^gctanGCTAN,0-9]', 'N', item[3])
        item[4]  = re.sub(r'[^gctanGCTAN,0-9]', 'N', item[4])
        line_out = '\t'.join(item)
        
        print(line_out)
