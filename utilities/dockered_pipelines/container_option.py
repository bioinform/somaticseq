import sys, argparse, os, re
import uuid
import utilities.split_Bed_into_equal_regions as split_bed
from somaticseq._version import __version__ as VERSION
from pathlib import Path




def container_params( container_image, tech='docker', files=[], extra_args='', singularity_image_loc='docker://' ):
    
    file_Paths    = [ Path(i) for i in files ]
    file_names    = [ i.name for i in file_Paths ]
    file_dirs     = [ i.parent for i in file_Paths ]
    file_abs_dirs = [ i.absolute().parent for i in file_Paths ]
    random_dirs   = [ '/'+uuid.uuid4().hex for i in files ]

    fileDict = {}
    
    for file_i, path_i, filename_i, dir_i, abs_dir_i, random_dir_i in zip(files, file_Paths, file_names, file_dirs, file_abs_dirs, random_dirs):
        fileDict[ file_i ] = {'filepath': path_i, 'filename': filename_i, 'dir': dir_i, 'abs_dir': abs_dir_i, 'mount_dir': random_dir_i, 'mount_path': os.path.join(random_dir_i, filename_i) }
    
    
    if tech == 'docker':
        
        MOUNT_STRING = ''
        for file_i in fileDict:
            sys_dir = fileDict[ file_i ][ 'abs_dir' ]
            container_dir = fileDict[ file_i ][ 'mount_dir' ]
            MOUNT_STRING = MOUNT_STRING + f' -v {sys_dir}:{container_dir}'
        
        container_string = f'docker run {MOUNT_STRING} -u $UID --rm {extra_args} {container_image}'
    
    
    elif tech == 'singularity':
        
        MOUNT_STRING = ''
        for file_i in fileDict:
            sys_dir = fileDict[ file_i ][ 'abs_dir' ]
            container_dir = fileDict[ file_i ][ 'mount_dir' ]
            MOUNT_STRING = MOUNT_STRING + f' --bind {sys_dir}:{container_dir}'
        
        container_string = f'singularity exec --cleanenv {MOUNT_STRING} {extra_args} {singularity_image_loc}{container_image}'
    
    
    return container_string, fileDict
