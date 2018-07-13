import os
import argparse
import string
import re

import time
import datetime

import random

letters = list(string.ascii_lowercase)

parser = argparse.ArgumentParser(description='Prepare the environment and start the plugin',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('bamdir', type=str, default="/results/analysis/output/Home", help='path of the manipulation')
parser.add_argument('traindir', type=str, default="/tmp/sanefalcontrain", help='path of the train subtree')

args = parser.parse_args()
bamdir = args.bamdir
traindir = args.traindir



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
            if f.endswith(".bam"):
                files_to_link.append(os.path.join(root, f))
            if root.count(os.sep) >= 1:  # stop at first level
                del subdir[:]

    # print(len(files_to_link))
    # for i in files_to_link:
    #     print(i)

    manip_list = list(set([f.split('/')[-2] for f in files_to_link]))  # list of manip with .bam files
    return files_to_link, manip_list


def prepare_batch_folders(manip_list):
    """Yield successive 5-sized chunks from manip_list."""
    for i in range(0, len(manip_list), 5):
        yield manip_list[i:i + 5]


files_to_link, manip_list = list_files_to_use(bamdir)

print("{} bam files in {} manips".format(len(files_to_link), len(manip_list)))

for batch in prepare_batch_folders(manip_list):
    print('batch')
    print(batch)

try:
    lst = sorted(os.listdir(traindir))[-1]
except IndexError:
    lst = None

exit(0)

if not lst:
    workingdir = os.path.join(traindir, 'a')
else:
    workingdir = os.path.join(traindir, letters[letters.index(lst) + 1 % len(letters)])

os.mkdir(workingdir)
# print('Training subfolder created')

runpath = ''
for fname in files_to_link:
    run = fname.split(bamdir)[1].split('/')[1]
    runpath = os.path.join(workingdir, run)
    if not os.path.isdir(runpath):
        os.mkdir(runpath)
    os.symlink(fname, os.path.join(runpath, fname.split('/')[-1]))

print(runpath)
