#!/usr/bin/env python3

import argparse
import logging
import subprocess
from datetime import datetime
from functools import partial
from multiprocessing import Pool


def seconds_to_readable(second_i):

    days = second_i // 86400
    hours = second_i // 3600 % 24
    minutes = second_i // 60 % 60
    seconds = second_i % 60

    out_time = {
        "days": int(days),
        "hours": int(hours),
        "minutes": int(minutes),
        "seconds": round(seconds, 3),
    }

    out_line = ""
    for item_i in ("days", "hours", "minutes", "seconds"):
        if out_time[item_i] > 0:
            out_line = out_line + " {} {}".format(out_time[item_i], str(item_i))

    return out_line.lstrip(" ")


def run_script(bash_script, shell="bash"):

    logger = logging.getLogger(run_script.__name__)

    cmd_line = "{} {}".format(shell, bash_script)
    logger.info(cmd_line)

    start_time = datetime.now()
    returnCode = subprocess.call(cmd_line, shell=True, bufsize=0)
    end_time = datetime.now()
    run_time = end_time - start_time

    time_line = seconds_to_readable(run_time.total_seconds())

    logger.info(
        "FINISHED RUNNING {} in {} with an exit code of {}.".format(
            bash_script, time_line, returnCode
        )
    )

    return returnCode


def run_workflows(list_of_ListOfTasks, threads=1, shell="bash"):
    """
    The input should be a list of "list of tasks."
    Each task shall be executed with the syntax of "bash TASK.sh".
    The order of the "list of tasks" imply precedence, i.e., all tasks in the 1st list must be complete before any task in the 2nd list may proceed.
    """

    pool = Pool(processes=threads)

    run_script_i = partial(run_script, shell=shell)

    for list_of_tasks in list_of_ListOfTasks:
        bash_async = pool.map_async(run_script_i, list_of_tasks)
        outputs = bash_async.get()

    pool.close()

    return True


def cumsum(num_in_each_list):
    """Calculate cumulative sum of a list like np.cumsum"""

    summed_list = []
    n = len(num_in_each_list)
    summed_list = [sum(num_in_each_list[0:i]) for i in range(0, n + 1)]

    return summed_list[1:]


def partition_list_to_lists(super_list, num_in_each_list):
    """
    Given a list, will partition it into multiple lists according to num_in_each_list, i.e.,
    Ex.
    super_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    num_in_each_list = [3, 6, 1], which means first 3 items go to 1st list, next 6 items to go 2nd list, and the next 1 item to go the 3rd list
    Return [ [1, 2, 3], [4, 5, 6, 7, 8, 9], [10] ]
    """
    assert sum(num_in_each_list) == len(super_list)

    num_start = cumsum(num_in_each_list)
    num_start.insert(0, 0)
    num_start = num_start[:-1]

    list_of_list = []
    for i, j in zip(num_start, num_in_each_list):
        list_i = super_list[i : i + j]
        list_of_list.append(list_i)

    return list_of_list


def run():

    parser = argparse.ArgumentParser(
        description="This module has two use cases:\n1) the first use case is simple, e.g., run the 4 scripts using 2 threads in parallel: run_workflows.py --scripts 1.sh 2.sh 3.sh 4.sh -nt 2.\n2) the second use case is more complex, e.g., with -nt 3 threads, the following 10 scripts have to be completed in the following orders: the first 3 need to complete first, and then the next 4 need to complete, and then finally the next 3 will need to complete: run_workflows.py -scripts 1.sh 2.sh 3.sh 4.sh 5.sh 6.sh 7.sh 8.sh 9.sh 10.sh -parts 3 4 3 -nt 3",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # INPUT FILES and Global Options
    parser.add_argument("-scripts", "--list-of-scripts", nargs="*", type=str)
    parser.add_argument("-nt", "--threads", default=1, type=int)
    parser.add_argument("-parts", "--partition-numbering", nargs="*", type=int)
    parser.add_argument("-sh", "--shell", default="bash", type=str)

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    FORMAT = "%(levelname)s %(asctime)-15s %(name)-20s %(message)s"
    logger = logging.getLogger("run_workflow.py")

    logger.setLevel(logging.DEBUG)
    logging.basicConfig(level=logging.INFO, format=FORMAT)

    args = run()

    if args.partition_numbering:
        list_of_workflows = partition_list_to_lists(
            args.list_of_scripts, args.partition_numbering
        )
        run_workflows(list_of_workflows, args.threads, args.shell)

    else:
        run_workflows(
            [
                args.list_of_scripts,
            ],
            args.threads,
            args.shell,
        )
