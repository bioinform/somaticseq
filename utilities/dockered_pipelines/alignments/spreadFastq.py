#!/usr/bin/env python3

import sys, argparse, os, re
import subprocess
import uuid
from pathlib import Path
from datetime import datetime
import utilities.dockered_pipelines.container_option as container
import genomicFileHandler.concat as concat

from somaticseq._version import __version__ as VERSION

ts = re.sub(r'[:-]', '.', datetime.now().isoformat(sep='.', timespec='milliseconds') )


DEFAULT_PARAMS = {'tabix_image'             : 'lethalfang/tabix:1.10',
                  'MEM'                     : 4,
                  'output_directory'        : os.curdir,
                  'action'                  : 'echo',
                  'extra_docker_options'    : '',
                  'script'                  : 'mergeFastqs.{}.cmd'.format(ts),
                  'action'                  : 'echo',
                  'threads'                 : 1,
                  }

