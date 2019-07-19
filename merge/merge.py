# This is the replacement for the following:
# merge.sh, mergeSubs.sh, mergeAntiSub.sh

import configparser
import concurrent.futures
import os
import logging
import re
import collections
from file_manager import FileManager
from multi.multithread import launch_multithreads as lt
from .utils import Utils

logger = logging.getLogger("merge")

chromosomes = range(1, 23)


def _merge(files, subdir, chrom):
    logger.debug("merging chrom {} start files in {}".format(chrom, subdir))
    outfile = os.path.join(subdir, "merge.{}".format(chrom))
    if os.path.isfile(outfile):
        logger.debug("merge. file {} already there. Skipping...".format(outfile))
        return
    data = Utils.read_all_files(files)
    logger.debug("merge into {}".format(outfile))
    Utils.sort_and_write(data, outfile)


def merge(files_to_merge):
    runs = Utils.prepare_jobs_arguments(files_to_merge)
    logger.debug('Submitting {} runs to merge'.format(len(runs)))
    lt(runs, _merge)


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
    lt(runs, _merge_subs)


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
    lt(runs, _merge_anti_subs)


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
