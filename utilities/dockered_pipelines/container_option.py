import sys, argparse, os, re
import uuid
import utilities.split_Bed_into_equal_regions as split_bed
import somaticseq._version.__version__ as VERSION
from pathlib import Path

def container_header( container, tech='docker', file_paths=[], extra_args='' ):
    
    
    
    
    if tech == 'docker':
        pass
        
    elif tech == 'singularities':
        pass
        
    pass
