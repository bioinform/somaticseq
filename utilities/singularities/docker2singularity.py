#!/usr/bin/env python3

import re, os, sys, subprocess

def docker2singularity(line_i):
    
    docker_cmd = re.search(r'^(.*)(docker\s+run)\s+(.+)\s+([\w/]+:[^\s]+)\s*(.*)', line_i)
    
    if docker_cmd:
        
        content_prior  = docker_cmd.groups()[0]
        docker_run     = docker_cmd.groups()[1]
        docker_options = docker_cmd.groups()[2]
        docker_image   = docker_cmd.groups()[3]
        docker_command = docker_cmd.groups()[4]
        
        mounts  = re.findall(r'(-v|--volume)[\s=]+([^\s:]+:[^\s:]+)', docker_options)
        mount_i = ''
        for i in mounts:
            mount_i = mount_i + '--bind ' + i[1] + ' '
        
        workdir_i = re.search(r'(-w|--workdir)[\s=]+([^\s:]+)', docker_options)
        if workdir_i:
            workdir_i = '--pwd ' + workdir_i.groups()[1]
        else:
            workdir_i = ''
        
        line_j = '{PRIOR}singularity exec {MOUNT} {WORKDIR} docker://{IMAGE} {CMD}\n'.format(PRIOR=content_prior, MOUNT=mount_i, WORKDIR=workdir_i, IMAGE=docker_image, CMD=docker_command)
        
        #line_k = line_j.encode('unicode_escape').decode()
        
        return line_j
        
    else:
        return line_i



sh_files = subprocess.getoutput( 'find ../dockered_pipelines -name "*.sh"' )
sh_files = sh_files.split('\n')

for file_i in sh_files:
    
    file_j = re.sub(r'../dockered_pipelines', '.', file_i)
    directory_j = os.sep.join( file_j.split(os.sep)[:-1] )
    os.system( 'mkdir -p {}'.format( directory_j ) )
    
    with open(file_i) as dockerFile, open(file_j, 'w') as singularityFile:
        
        line_i = dockerFile.readline()
        
        while line_i:
            line_j = docker2singularity(line_i)
            singularityFile.write( line_j )
        
            line_i = dockerFile.readline() 
