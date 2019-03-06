# This is the replacement for the following:
# merge.sh, mergeSubs.sh, mergeAntiSub.sh

import configparser
import concurrent.futures
import threading
import os
import logging
import re
import collections

logging.basicConfig(format='%(asctime)s - %(levelname)s: %(message)s',
                    filename='sanefalcon.log', filemode='w', level=logging.DEBUG)

logger = logging.getLogger("merge")

MAX_THREAD_NUMBER = 15
chromosomes = range(1, 23)


def sort_and_write(data, outfile):
    data = list(set(data))  # remove duplicates and sort
    data.sort()
    to_be_written = map(lambda x: str(x) + "\n", data)
    with open(outfile, 'w') as out:
        out.writelines(to_be_written)


def find_all_manips(trainfolder):
    manips = {}
    subfolders = [f.path for f in os.scandir(trainfolder) if f.is_dir()]
    for subdir in subfolders:
        manips[subdir] = [os.path.basename(f.path) for f in os.scandir(subdir) if f.is_dir()]
    return manips


def find_merge_files_in_subdirectories(trainfolder):
    merge_files = []
    pattern_name = re.compile("merge.\d{1,2}")
    pattern_subdir = re.compile("/[a-z]/")
    for root, subdirs, files in os.walk(trainfolder):
        for fname in files:
            filename = os.path.join(root, fname)
            if re.match(pattern_name, fname) and re.search(pattern_subdir, filename):
                merge_files.append(filename)
    return merge_files

#
# def search_manip_name(manips, fname):
#     base = os.path.basename(fname)
#     logger.debug("Searching for manip {}".format(base))
#     for subdir, manip_names in manips.items():
#         logger.debug("{} - {}".format(subdir, manip_names))
#         basenames = [os.path.basename(f) for f in manip_names]
#         logger.warning("any: {}".format(any([os.path.commonprefix([bn, base]) == bn for bn in basenames])))
#         if any([os.path.commonprefix([bn, base]) == bn for bn in basenames]):
#             return subdir
#         else:
#             return None


def prepare_file_lists(trainfolder):
    """

    :param trainfolder:
    :return: {'chr1': {'a': [chr1.start.fwd, chr1.start.rev], 'b': [...]}}
    """
    res = {}
    manips = find_all_manips(trainfolder)
    files_dic = dict.fromkeys(list(map(str, chromosomes)), list())
    all_start_files = []
    for root, subdirs, files in os.walk(trainfolder):
        for fname in files:
            if fname.endswith('.fwd') or fname.endswith('.rev'):
                all_start_files.append(os.path.join(root, fname))

    for chrom in chromosomes:
        string_pattern = "\.bam\.{}\.start".format(chrom)
        pattern = re.compile(string_pattern)
        files_dic[str(chrom)] = list(filter(lambda x: re.search(pattern, x), all_start_files))

    for chrom, chrom_start_files in files_dic.items():
        logger.debug('Preparing start_files from chrom {}'.format(chrom))
        files_per_subdir = {}
        for subdir, manip_names in manips.items():
            start_per_manip = []
            for n in manip_names:
                start_per_manip.extend(list(filter(lambda x: re.search(n, x), chrom_start_files)))
            if subdir in files_per_subdir:
                files_per_subdir[subdir].extend(start_per_manip)
            else:
                files_per_subdir[subdir] = start_per_manip
        res[chrom] = files_per_subdir

    return res


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
    data = read_all_files(files)
    outfile = os.path.join(subdir, "merge.{}".format(chrom))
    logger.debug("merge into {}".format(outfile))
    sort_and_write(data, outfile)
    # print("{} completed on {} (chrom {})".format(threading.current_thread(), subdir, chrom))


def _prepare_jobs_arguments(files_to_merge):
    runs = []
    for chrom, dic in files_to_merge.items():
        for subdir, files in dic.items():
            runs.append((files, subdir, chrom))
    return runs


def merge(files_to_merge):
    runs = _prepare_jobs_arguments(files_to_merge)
    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_THREAD_NUMBER) as executor:
        jobs = {executor.submit(_merge, *run): run for run in runs}
        for job in concurrent.futures.as_completed(jobs):
            try:
                _ = job.result()
            except Exception as ex:
                logger.error('Future Exception {}'.format(ex.__cause__))


def _merge_subs(chrom, files, trainfolder):
    # print("{} launched on {} ".format(threading.current_thread(), chrom))
    logger.debug("merge subs chrom {}".format(chrom))
    data = read_all_files(files)
    outfile = os.path.join(trainfolder, "merge.{}".format(chrom))
    sort_and_write(data, outfile)
    # print("{} completed on {}".format(threading.current_thread(), chrom))


def merge_subs(trainfolder):
    merge_files = find_merge_files_in_subdirectories(trainfolder)
    tmp = collections.defaultdict(list)
    for f in merge_files:
        tmp[f.split('.')[1]].append(f)

    runs = [(k, v, trainfolder) for k, v in tmp.items()]
    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_THREAD_NUMBER) as executor:
        jobs = {executor.submit(_merge_subs, *run): run for run in runs}
        for job in concurrent.futures.as_completed(jobs):
            try:
                _ = job.result()
            except Exception as ex:
                logger.error('Future Exception {}'.format(ex.__cause__))
    return merge_files


def _merge_anti_subs(folder, chrom, files):
    # print("{} launched on {} ".format(threading.current_thread(), chrom))
    logger.debug("merge anti subs chrom: {}, folder: {}".format(chrom, folder))
    data = read_all_files(files)
    outfile = os.path.join(folder, "anti.{}".format(chrom))
    sort_and_write(data, outfile)
    # print("{} completed on {}".format(threading.current_thread(), chrom))


def merge_anti_subs(merge_files, trainfolder):
    """
    input: [sanefalcontrain/b/merge.chr1, sanefalcontrain/c/merge.chr1, ...]
    output: sanefalcontrain/a/anti.chr1
    :param trainfolder:
    :return:
    """
    tmp = collections.defaultdict(list)
    for f in merge_files:
        tmp[f.split('.')[1]].append(f)
    # print(tmp)
    subfolders = [f.path for f in os.scandir(trainfolder) if f.is_dir()]

    runs = []
    for folder in subfolders:
        pattern = re.compile('(?!{})'.format(folder))  # matching "not folder"
        reduced_tmp = {k: list(filter(lambda x: re.match(pattern, x), v)) for k, v in tmp.items()}
        run = [(folder, k, v) for k, v in reduced_tmp.items()]
        runs.extend(run)

    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_THREAD_NUMBER) as executor:
        jobs = {executor.submit(_merge_anti_subs, *run): run for run in runs}
        for job in concurrent.futures.as_completed(jobs):
            try:
                _ = job.result()
            except Exception as ex:
                logger.error('Future Exception {}'.format(ex.__cause__))


def merge_all(trainfolder):
    logger.info("starting merge all")
    dic = prepare_file_lists(trainfolder)
    merge(dic)
    logger.debug("merge done")
    merge_files = merge_subs(trainfolder)
    logger.debug("merge_subs done")
    merge_anti_subs(merge_files, trainfolder)
    logger.debug("merge_anti_subs done")
    logger.info("merging completed")


if __name__ == "__main__":
    conf_file = 'sanefalcon.conf'
    config = configparser.ConfigParser()
    config.read(conf_file)

    samtools = config['default']['samtools']
    datafolder = config['default']['datafolder']
    trainfolder = config['default']['trainfolder']

    # files_to_merge = prepare_file_lists(trainfolder)
    # for k, v in files_to_merge['22'].items():
    #     print(k)
    #     for n in sorted(v):
    #         print(n)
    # merge(files_to_merge)
    # merge_files = merge_subs(trainfolder)
    # merge_files = find_merge_files_in_subdirectories(trainfolder)
    # merge_anti_subs(merge_files, trainfolder)
    merge_all(trainfolder)
    # merge(dic)
