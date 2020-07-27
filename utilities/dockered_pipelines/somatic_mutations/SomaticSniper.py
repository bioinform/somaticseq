import sys, argparse, os, re
import subprocess
from datetime import datetime
import utilities.dockered_pipelines.container_option as container
from somaticseq._version import __version__ as VERSION

ts = re.sub(r'[:-]', '.', datetime.now().isoformat() )


DEFAULT_PARAMS = {'somaticsniper_image'     : 'lethalfang/somaticsniper:1.0.5.0-2',
                  'MEM'                     : '4G',
                  'threads'                 : 1,
                  'normal_bam'              : None,
                  'tumor_bam'               : None,
                  'genome_reference'        : None,
                  'inclusion_region'        : None,
                  'output_directory'        : os.curdir,
                  'outfile'                 : 'VarScan2.vcf',
                  'action'                  : 'echo',
                  'extra_arguments'         : '',
                  'extra_docker_options'    : '',
                  'script'                  : 'varscan2.{}.cmd'.format(ts),
                  'min_MQ'                  : 1,
                  'min_BQ'                  : 20,
                  }

