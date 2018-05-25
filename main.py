import os
import argparse
import string

import random

letters = list(string.ascii_lowercase)
parser = argparse.ArgumentParser(description='Prepare the environment and start the plugin',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('bamdir', type=str, default="/results/analysis/output/Home", help='path of the manipulation')
parser.add_argument('traindir', type=str, default="/tmp", help='path of the train subtree')

args = parser.parse_args()

bamdir = args.bamdir
traindir = args.traindir

files_to_link = []
for root, subdir, files in os.walk(bamdir):
    for f in files:
        files_to_link.append(os.path.join(root, f))

try:
    lst = sorted(os.listdir(traindir))[-1]
except IndexError:
    lst = None

if not lst:
    workingdir = os.path.join(traindir, 'a')
else:
    workingdir = os.path.join(traindir, letters[letters.index(lst) + 1 % len(letters)])

os.mkdir(workingdir)
# print('Training subfolder created')

for fname in files_to_link:
    run = fname.split(bamdir)[1].split('/')[1]
    runpath = os.path.join(workingdir, run)
    if not os.path.isdir(runpath):
        os.mkdir(runpath)
    os.symlink(fname, os.path.join(runpath, fname.split('/')[-1]))

print(runpath)