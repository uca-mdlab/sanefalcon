import os
import re
import csv
from collections import Counter, defaultdict
import logging

logging.basicConfig(format='%(asctime)s - %(levelname)s: %(message)s',
                    filename='sanefalcon.log', filemode='w', level=logging.DEBUG)

logger = logging.getLogger(__name__)


def get_profile_file(proffile):
    l = []
    with open(proffile, 'r') as in_:
        reader = list(csv.reader(in_, delimiter=','))
        assert len(reader) == 1
        for row in reader:
            l = list(map(float, row))
    return l


def get_profile_files_per_chrom(fnames, stream='up'):
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
        d[name] = get_profile_file(fname)
    return d


def create_upstream(fnames):
    upstream = defaultdict(list)
    tmp = defaultdict(list)
    for chrom, fwd, rev in get_profile_files_per_chrom(fnames, stream='up'):
        d1 = combine_profile_files(fwd)
        d2 = combine_profile_files(rev)
        for k, v in d1.items():
            tmp[k].append(v)
        for k, v in d2.items():
            tmp[k].append(v)
    for k, v in tmp.items():
        logger.debug('Upstream {}: {} files with {} entries each.'.format(k, len(v), set([len(x) for x in v])))
        upstream[k] = [sum(x) for x in zip(*v)]
    return upstream


def create_downstream(fnames):
    downstream = defaultdict(list)
    tmp = defaultdict(list)
    for chrom, ifwd, irev in get_profile_files_per_chrom(fnames, stream='down'):
        d1 = combine_profile_files(ifwd)
        d2 = combine_profile_files(irev)
        for k, v in d1.items():
            tmp[k].append(v)
        for k, v in d2.items():
            tmp[k].append(v)
    for k, v in tmp.items():
        logger.debug('Downstream {}: {} files with {} entries each.'.format(k, len(v), set([len(x) for x in v])))
        downstream[k] = [sum(x) for x in zip(*v)]
    return downstream


def create_streams(fm):
    fnames = [os.path.join(fm.profilefolder, fname) for fname in os.listdir(fm.profilefolder)]
    streams = {}
    upstream = create_upstream(fnames)
    downstream = create_downstream(fnames)
    # print(upstream['s48_np_090.bam'][0], downstream['s48_np_090.bam'][0])
    # print(upstream['s48_np_090.bam'][-1], downstream['s48_np_090.bam'][-1])
    import matplotlib.pyplot as plt
    x = range(294)
    for k in upstream.keys():
        up = upstream[k]
        down = downstream[k]
        up.reverse()
        y = up + down
        plt.plot(x, y)
    plt.show()
    return streams


if __name__ == '__main__':
    import configparser
    from file_manager import FileManager

    config = configparser.ConfigParser()
    config.read('local.conf')
    f = FileManager(config)
    streams = create_streams(f)
    # for k, v in streams.items():
    #     print(k, v)

