import os
import re
import csv
from collections import Counter, defaultdict
import logging

logging.basicConfig(format='%(asctime)s - %(levelname)s: %(message)s',
                    filename='sanefalcon.log', filemode='w', level=logging.DEBUG)

logger = logging.getLogger(__name__)


def sum_profile_file(proffile):
    s = 0
    with open(proffile, 'r') as in_:
        reader = csv.reader(in_, delimiter=',')
        for row in reader:
            s += sum(list(map(float, row)))
    return s


def get_profile_files_per_chrom(fm, stream='up'):
    fnames = [os.path.join(fm.profilefolder, fname) for fname in os.listdir(fm.profilefolder)]
    for chrom in range(1, 23):
        pattern = re.compile(r"\.{}\.".format(chrom))
        chrom_files = list(filter(lambda x: re.search(pattern, x), fnames))
        if stream == 'up':
            fwd = [f for f in chrom_files if f.endswith('.fwd')]
            rev = [f for f in chrom_files if f.endswith('.rev')]
            yield chrom, fwd, rev
        else:
            ifwd = [f for f in chrom_files if f.endswith('.ifwd')]
            irev = [f for f in chrom_files if f.endswith('.irev')]
            yield chrom, ifwd, irev


def combine_profile_files(arr):
    d = {}
    pattern = re.compile(r"\w*\.bam")
    for fname in arr:
        name = re.match(pattern, os.path.basename(fname)).group()
        d[name] = sum_profile_file(fname)
    return d


def create_upstream(fm):
    res = defaultdict(dict)
    for chrom, fwd, rev in get_profile_files_per_chrom(fm, stream='up'):
        d1 = combine_profile_files(fwd)
        d2 = combine_profile_files(rev)
        c = Counter(d1)
        c.update(Counter(d2))
        for k, v in dict(c).items():
            res[k].update({chrom: v})
    return res


def create_downstream(fm):
    res = defaultdict(dict)
    for chrom, ifwd, irev in get_profile_files_per_chrom(fm, stream='down'):
        d1 = combine_profile_files(ifwd)
        d2 = combine_profile_files(irev)
        c = Counter(d1)
        c.update(Counter(d2))
        for k, v in dict(c).items():
            res[k].update({chrom: v})
    return res


def create_streams(fm):
    streams = {}
    upstream = create_upstream(fm)
    downstream = create_downstream(fm)
    for k, v in upstream.items():
        nums = [x for _, x in v.items()]
        nums.reverse()
        logger.debug('Stream up (rev): {}. length {}. [0]: {}, [-1]: {}'.format(k, len(nums), nums[0], nums[-1]))
        streams[k] = nums[:-1]
    for k, v in downstream.items():
        nums = [x for _, x in v.items()]
        logger.debug('Stream down: {}. length {}. [0]: {}, [-1]: {}'.format(k, len(nums), nums[0], nums[-1]))
        streams[k].extend(nums)
    return streams


if __name__ == '__main__':
    # fname = '/tmp/testsuite/profiles/p2.bam.15.start.fwd.15.fwd'
    # s = sum_profile_file(fname)
    # print(s)
    # exit()
    import configparser
    from file_manager import FileManager

    config = configparser.ConfigParser()
    config.read('tests/data/test.conf')
    f = FileManager(config)
    streams = create_streams(f)
    for k, v in streams.items():
        print(k, v)

