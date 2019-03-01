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

def readfile(bamlist):
    with open (bamlist) as f:
        myfile=[x.strip() for x in f.readlines()]
    f.close()

    return myfile

def list_files_to_use(bamlist,bamdir):
    list_files=[x.split('/')[-1] for x in bamlist]
    files_to_link=[]
    for root,subdir,files in os.walk(bamdir):
        for f in files:
            if f in list_files:
                files_to_link.append(os.path.join(root,f))

        # files_to_link=[os.path.join(root,f) for f in files if f in list_files]


    # files_to_link=["/".join((bamdir,x)) for x in bamlist if os.path.isfile("/".join((bamdir,x)))]

    manip_list = list(set([f.split('/')[-2] for f in files_to_link if os.path.isfile(f)]))  # list of manip with .bam files
    logger.info("Found {} manip with {} bam files".format(len(manip_list), len(files_to_link)))

    return files_to_link, manip_list


def prepare_batches(manip_list):
    # Yield successive 5-sized chunks from manip_list.
    for i in range(0, len(manip_list), 5):
        yield manip_list[i:i + 5]


def prepare_train_folder(bamlist, bamdir, traindir):
# def prepare_train_folder(bamlist, traindir):
    letters = list(string.ascii_lowercase)
    files_to_link, manip_list = list_files_to_use(readfile(bamlist),bamdir)

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
            for fname in files:
                run = fname.split(bamdir)[1].split('/')[1]
                runpath = os.path.join(workingdir, run)
                if not os.path.isdir(runpath):
                    os.mkdir(runpath)
                os.symlink(fname, os.path.join(runpath, fname.split('/')[-1]))
                fname=fname+".bai"
                os.symlink(fname, os.path.join(runpath, fname.split('/')[-1]))
        logger.info("Batches created with symlinks")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Prepare the environment and start the plugin',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('bamdir', type=str, default="/results/analysis/output/Home", help='path of the manipulation')
    parser.add_argument('bamlist', type=str, help='list of bamfiles to use, with run names and sample names')
    parser.add_argument('traindir', type=str, default="/tmp/sanefalcontrain", help='path of the train subtree')

    args = parser.parse_args()
    bamdir = args.bamdir
    bamlist = args.bamlist
    traindir = args.traindir
    #
    prepare_train_folder(bamlist,bamdir, traindir)
