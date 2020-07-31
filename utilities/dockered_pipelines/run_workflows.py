#!/usr/bin/env python3

import sys, argparse, os
import logging
import subprocess
from multiprocessing import Pool

def run_workflows( list_of_ListOfTasks, threads = 1 ):
    '''
    The input should be a list of "list of tasks."
    Each task shall be executed with the syntax of "bash TASK.sh".
    The order of the "list of tasks" imply precedence, i.e., all tasks in the 1st list must be complete before any task in the 2nd list may proceed.
    '''
    
    return True
