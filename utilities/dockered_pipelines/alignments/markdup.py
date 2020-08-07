#!/usr/bin/env python3

import sys, argparse, os, re
import subprocess
import uuid
from pathlib import Path
from datetime import datetime
import utilities.dockered_pipelines.container_option as container
from somaticseq._version import __version__ as VERSION

ts = re.sub(r'[:-]', '.', datetime.now().isoformat() )

DEFAULT_PARAMS = {'bwa_image'               : 'lethalfang/picard:2.22.7',
                  'MEM'                     : 16,
                  'output_directory'        : os.curdir,
                  'out_bam'                 : 'aligned.markdup.bam',
                  'action'                  : 'echo',
                  'extra_docker_options'    : '',
                  'extra_picard_arguments'  : '',
                  'threads'                 : 1,
                  'script'                  : 'markdup.{}.cmd'.format(ts),
                  }


def picard( input_parameters, tech='docker' ):


    # "Run" the script that was generated
    command_item = (input_parameters['action'], outfile)
    returnCode   = subprocess.call( command_item )

    return outfile
