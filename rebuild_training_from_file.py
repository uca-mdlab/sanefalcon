import os
import sys
from collections import defaultdict
import re


data_dir = '/home/mdlab/storage/data_nice/dedup'
trainingdir = '/home/mdlab/storage/sanefalcon/training'
batches = defaultdict(list)
with open(sys.argv[1], 'r') as infile:
    lines = [x.strip() for x in infile.readlines()]

for line in lines:
    if not re.match(r'[a-z]/.*\.bam', line):
        exit(f'Wrong file format. Got: {line}')
    name, f = line.split('/')
    batches[name].append(os.path.join(data_dir, f))

for k, v in batches.items():
    current_dir = os.path.join(trainingdir, k)
    try:
        os.makedirs(current_dir)
    except FileExistsError:
        pass

    for fname in v:
        os.symlink(fname, os.path.join(current_dir, os.path.basename(fname)))
        os.symlink(fname + '.bai', os.path.join(current_dir, os.path.basename(fname) + '.bai'))
