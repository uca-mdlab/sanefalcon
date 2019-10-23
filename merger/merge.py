import os
import logging
from collections import defaultdict
from multi.multithread import launch_multithreads as lt
from .utils import Utils

logger = logging.getLogger(__name__)

chromosomes = range(1, 23)


def _merge(files, subdir, chrom):
    outfile = os.path.join(subdir, "merge.{}".format(chrom))
    if os.path.isfile(outfile):
        logger.debug("_merge: {} already there. Skipping...".format(outfile))
    else:
        logger.debug("_merge: merging chrom {} start files in {} -> {}".format(chrom, subdir, outfile))
        data = Utils.read_all_files(files)
        Utils.sort_and_write(data, outfile)
    return subdir, outfile


def merge(files_to_merge):
    runs = Utils.prepare_jobs_arguments(files_to_merge)
    logger.debug('Submitting {} runs to merge'.format(len(runs)))
    list_of_tup = lt(runs, _merge)
    merge_files = defaultdict(list)
    for subdir, mergefile in list_of_tup:
        merge_files[subdir].append(mergefile)
    return merge_files


def _merge_subs(outfile, files):
    if os.path.isfile(outfile):
        logger.debug("merge_subs. file {} already there. Skipping...".format(outfile))
    else:
        logger.debug("merge subs file {}".format(outfile))
        data = Utils.read_all_files(files)
        Utils.sort_and_write(data, outfile)
    return outfile


def merge_subs(merged, trainfolder):
    runs = [(os.path.join(trainfolder, fname), list_) for fname, list_ in merged.items()]
    logger.debug('Submitting {} runs to merge_subs'.format(len(runs)))
    list_of_merge_files = lt(runs, _merge_subs)
    return list_of_merge_files


def _merge_anti_subs(folder, chrom, files):
    # print("{} launched on {} ".format(threading.current_thread(), chrom))
    outfile = os.path.join(folder, "anti.{}".format(chrom))
    if os.path.isfile(outfile):
        logger.debug("merge_anti_subs. file {} already there. Skipping...".format(outfile))
    else:
        logger.debug("merge anti subs chrom: {}, folder: {} -> {}".format(chrom, folder, outfile))
        data = Utils.read_all_files(files)
        Utils.sort_and_write(data, outfile)
    return outfile


def merge_anti_subs(antisubs):
    runs = [(path, chrom, files) for path, chrom, files in antisubs]
    logger.debug('Submitting {} runs to merge_anti_subs'.format(len(runs)))
    antisubfiles = lt(runs, _merge_anti_subs)

    # logger.info('Submitting {} runs to merge_anti_subs (sequential)'.format(len(runs)))
    # antisubfiles = []
    # lrun = len(runs)
    # for r in runs:
    #     antisub = _merge_anti_subs(*r)
    #     antisubfiles.append(antisub)
    #     logger.debug('Antisub done: {}. ({} to go...)'.format(antisub, lrun - len(antisubfiles)))
    return antisubfiles

