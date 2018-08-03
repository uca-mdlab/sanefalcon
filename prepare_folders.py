import os
import argparse
import string
import re

import time
import datetime

import random
import logging
logging.basicConfig(format='%(asctime)s - %(levelname)s: %(message)s',
                    filename='sanefalcon.log', filemode='w', level=logging.DEBUG)

logger = logging.getLogger(__name__)

letters = list(string.ascii_lowercase)


def list_files_to_use(bamdir):
    """

    :param bamdir: the root of all manips
    :return: list of bam files to use for training
    :return: list of manip names
    """
    first_date = datetime.date(2017, 1, 13)
    timestamp = time.mktime(first_date.timetuple())

    exclude = set([d for d in os.listdir(bamdir) if re.search('_tn_', d) or not re.search('Auto_user_', d) or
                   not re.search('DPNI', d) or os.path.getmtime(os.path.join(bamdir, d)) <= timestamp])

    files_to_link = []
    for root, subdir, files in os.walk(bamdir):
        subdir[:] = set(subdir) - exclude
        for f in files:
            if f.endswith(".bam") or f.endswith(".bai"):
                files_to_link.append(os.path.join(root, f))
            if root.count(os.sep) >= 1:  # stop at first level
                del subdir[:]

    manip_list = list(set([f.split('/')[-2] for f in files_to_link]))  # list of manip with .bam files
    logger.info("Found {} manip with {} bam files".format(len(manip_list), len(files_to_link)))
    return files_to_link, manip_list


def prepare_batches(manip_list):
    """Yield successive 5-sized chunks from manip_list."""
    for i in range(0, len(manip_list), 5):
        yield manip_list[i:i + 5]


def main(bamdir, traindir):
    files_to_link, manip_list = list_files_to_use(bamdir)

    batches = {}
    for num_batch, batch in enumerate(prepare_batches(manip_list)):
        batches[letters[num_batch]] = batch
        logger.debug("Batch {}: {}".format(letters[num_batch], batch))

    for batch_name, batch_list in batches.items():
        workingdir = os.path.join(traindir, batch_name)
        try:
            os.mkdir(workingdir)
            logger.debug("Created folder: {}".format(workingdir))
        except FileExistsError:
            logger.warning("Folder {} exists, skipping...".format(workingdir))

        for manip in batch_list:
            manip_regex = re.compile(manip)
            files = [f for f in files_to_link if re.search(manip_regex, f)]
            logger.debug("Batch {}: {}".format(batch_name, files))
            runpath = ''
            for fname in files:
                run = fname.split(bamdir)[1].split('/')[1]
                runpath = os.path.join(workingdir, run)
                if not os.path.isdir(runpath):
                    os.mkdir(runpath)
                os.symlink(fname, os.path.join(runpath, fname.split('/')[-1]))
        logger.info("Batches created with symlinks")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Prepare the environment and start the plugin',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('bamdir', type=str, default="/results/analysis/output/Home", help='path of the manipulation')
    parser.add_argument('traindir', type=str, default="/tmp/sanefalcontrain", help='path of the train subtree')

    args = parser.parse_args()
    bamdir = args.bamdir
    traindir = args.traindir

    # UNCOMMENT for complete run
    # main(bamdir, traindir)
    logger.warning("Skipping prepare_folder.py")
    print(traindir)

