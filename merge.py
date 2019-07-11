# This is the replacement for the following:
# merge.sh, mergeSubs.sh, mergeAntiSub.sh

import configparser
import concurrent.futures
import threading
import os
import logging
import re
import collections
from file_manager import FileManager
from collections import defaultdict

logging.basicConfig(format='%(asctime)s - %(levelname)s: %(message)s',
                    filename='sanefalcon.log', filemode='w', level=logging.DEBUG)

logger = logging.getLogger("merge")

MAX_THREAD_NUMBER = 10
MAX_JOB_NUMBER = 10
chromosomes = range(1, 23)


def launch_multithreads(runs, name='merge'):
    logger.info('Starting multithreaded {} '.format(name))
    if name == 'merge':
        fn = _merge
    elif name == 'merge_subs':
        fn = _merge_subs
    elif name == 'merge_anti_subs':
        fn = _merge_anti_subs
    else:
        exit('Unknown function ({}) for multithread'.format(name))

    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_THREAD_NUMBER) as executor:
        jobs = {}
        runs_left = len(runs)
        runs_iter = iter(runs)

        while runs_left:
            for run in runs_iter:
                job = executor.submit(fn, *run)
                jobs[job] = run
                if len(jobs) > MAX_JOB_NUMBER:
                    break

            for job in concurrent.futures.as_completed(jobs):
                runs_left -= 1
                _ = job.result()
                run = jobs[job]
                logger.debug('Ended job {}'.format(run))
                del jobs[job]
                break

    logger.info('Multithreaded {} ended'.format(name))


def sort_and_write(data, outfile):
    # data = list(set(data))  # remove duplicates and sort
    data.sort()
    to_be_written = map(lambda x: str(x) + "\n", data)
    with open(outfile, 'w') as out:
        out.writelines(to_be_written)


def read_all_files(file_list):
    data = []
    for f in file_list:
        data.extend([int(line.strip()) for line in open(f, 'r')])
    return data


def pick_manip_pairs(manips, list_):
    res = {}
    for subdir, manip_list in manips.items():
        for manip in manip_list:
            files = []
            manip_name = os.path.basename(manip)
            for path in list_:
                base = os.path.basename(path)
                if os.path.commonprefix([base, manip_name]) == manip_name:
                    files.append(path)
            res[manip] = files
    return res


def _merge(files, subdir, chrom):
    # print("{} launched on {} (chrom {})".format(threading.current_thread(), subdir, chrom))
    logger.debug("merging chrom {} start files in {}".format(chrom, subdir))
    outfile = os.path.join(subdir, "merge.{}".format(chrom))
    if os.path.isfile(outfile):
        logger.debug("merge. file {} already there. Skipping...".format(outfile))
        return
    data = read_all_files(files)
    logger.debug("merge into {}".format(outfile))
    sort_and_write(data, outfile)
    # print("{} completed on {} (chrom {})".format(threading.current_thread(), subdir, chrom))


def _prepare_jobs_arguments(files_to_merge):
    runs = []
    for subdir, dic in files_to_merge.items():
        for chrom in range(1, 23):
            fwd = dic['fwd'][chrom]
            rev = dic['rev'][chrom]
            files = fwd + rev
            runs.append((files, subdir, chrom))
    return runs


def merge(files_to_merge):
    runs = _prepare_jobs_arguments(files_to_merge)
    logger.debug('Submitting {} runs to merge'.format(len(runs)))
    launch_multithreads(runs, name='merge')


def _merge_subs(chrom, files, trainfolder):
    # print("{} launched on {} ".format(threading.current_thread(), chrom))
    logger.debug("merge subs chrom {}".format(chrom))
    outfile = os.path.join(trainfolder, "merge.{}".format(chrom))
    if os.path.isfile(outfile):
        logger.debug("merge_subs. file {} already there. Skipping...".format(outfile))
        return
    data = read_all_files(files)
    sort_and_write(data, outfile)
    # print("{} completed on {}".format(threading.current_thread(), chrom))


def merge_subs(merge_subs_files, trainfolder):
    tmp = collections.defaultdict(list)
    for f in merge_subs_files:
        tmp[f.split('.')[1]].append(f)  # grep .chrom_number
    runs = [(k, v, trainfolder) for k, v in tmp.items()]
    logger.debug('Submitting {} runs to merge_subs'.format(len(runs)))
    launch_multithreads(runs, name='merge_subs')


def _merge_anti_subs(folder, chrom, files):
    # print("{} launched on {} ".format(threading.current_thread(), chrom))
    logger.debug("merge anti subs chrom: {}, folder: {}".format(chrom, folder))
    outfile = os.path.join(folder, "anti.{}".format(chrom))
    if os.path.isfile(outfile):
        logger.debug("merge_anti_subs. file {} already there. Skipping...".format(outfile))
        return
    data = read_all_files(files)
    sort_and_write(data, outfile)
    # print("{} completed on {}".format(threading.current_thread(), chrom))


def merge_anti_subs(merge_subs_files, trainfolder):
    """
    input: [sanefalcontrain/b/merge.chr1, sanefalcontrain/c/merge.chr1, ...]
    output: sanefalcontrain/a/anti.chr1
    :param trainfolder:
    :return:
    """
    subfolders = set()
    tmp = collections.defaultdict(list)
    for f in merge_subs_files:
        subfolders.add(os.path.join(trainfolder, os.path.basename(os.path.dirname(f))))
        tmp[f.split('.')[1]].append(f)

    runs = []
    for folder in list(subfolders):
        pattern = re.compile('(?!{})'.format(folder))  # matching "not folder"
        reduced_tmp = {k: list(filter(lambda x: re.match(pattern, x), v)) for k, v in tmp.items()}
        run = [(folder, k, v) for k, v in reduced_tmp.items()]
        runs.extend(run)

    logger.debug('Submitting {} runs to merge_anti_subs'.format(len(runs)))
    launch_multithreads(runs, name='merge_anti_subs')


def merge_all(fm):
    logger.info("starting merge all")
    merge_file_list = fm.get_start_files_per_subdir()
    merge(merge_file_list)
    logger.debug("merge done")
    merge_subs_files = fm.find_merge_files_in_subdirectories()
    merge_subs(merge_subs_files, fm.trainfolder)
    logger.debug("merge_subs done")
    merge_anti_subs(merge_subs_files, fm.trainfolder)
    logger.debug("merge_anti_subs done")
    logger.info("merging completed")


if __name__ == "__main__":
    conf_file = 'tests/data/test.conf'
    config = configparser.ConfigParser()
    config.read(conf_file)

    samtools = config['default']['samtools']
    datafolder = config['default']['datafolder']
    trainfolder = config['default']['trainfolder']
    rspfolder = config['default']['rspfolder']
    fm = FileManager(config)
    merge_all(fm)
