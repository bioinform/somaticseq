#!/usr/bin/env python3

import os


def leftAlign(infile, outfile, ref, gatk3):

    assert infile != outfile

    exit_code = os.system( '''java -jar {} -T LeftAlignAndTrimVariants -R {} --variant {} | egrep -v '^[0-9]+ variants|^INFO' > {}'''.format(gatk3, ref, infile, outfile) )
    
    assert exit_code
    
    return outfile
