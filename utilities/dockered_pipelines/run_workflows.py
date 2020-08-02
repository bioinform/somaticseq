#!/usr/bin/env python3

import sys, argparse, os
import logging
import subprocess
from datetime import datetime
from multiprocessing import Pool

def seconds_to_readable( second_i ):

    days    = second_i // 86400
    hours   = second_i // 3600 % 24
    minutes = second_i // 60 % 60
    seconds = second_i % 60

    out_time = {'days': int(days), 'hours': int(hours), 'minutes': int(minutes), 'seconds': round(seconds, 3)}

    out_line = ''
    for item_i in ('days', 'hours', 'minutes', 'seconds'):
        if out_time[ item_i ] > 0:
            out_line = out_line + ' {} {}'.format(out_time[ item_i ], str(item_i))
            
    return out_line.lstrip(' ')



def run_script( bash_script ):
    
    logger = logging.getLogger(run_script.__name__)
    
    cmd_items = ('bash', bash_script)
    logger.info(' '.join(cmd_items) )
    
    start_time = datetime.now()
    returnCode = subprocess.call( cmd_items )
    end_time   = datetime.now()
    run_time   = end_time - start_time
    
    time_line  = seconds_to_readable( run_time.total_seconds() )
    
    logger.info('FINISHED RUNNING {} in {} with an exit code of {}.'.format(bash_script, time_line, returnCode) )

    return returnCode




def run_workflows( list_of_ListOfTasks, threads = 1 ):
    '''
    The input should be a list of "list of tasks."
    Each task shall be executed with the syntax of "bash TASK.sh".
    The order of the "list of tasks" imply precedence, i.e., all tasks in the 1st list must be complete before any task in the 2nd list may proceed.
    '''
    pool = Pool(processes = threads)
    
    for list_of_tasks in list_of_ListOfTasks:
        bash_async = pool.map_async(run_script, list_of_tasks)
        outputs    = bash_async.get()

    pool.close()

    return True
