import sys, argparse, os, re
import uuid
from datetime import datetime
from shutil import move

import utilities.split_Bed_into_equal_regions as split_bed
import somaticseq._version

VERSION = somaticseq._version.__version__

ts = re.sub(r'[:-]', '.', datetime.now().isoformat() )


def tumor_normal(input_parameters):
    
    pass
