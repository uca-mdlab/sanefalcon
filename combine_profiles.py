import os
import re
import csv
from collections import Counter, defaultdict


def sum_profile_file(proffile):
    s = 0
    with open(proffile, 'r') as in_:
        reader = csv.reader(in_, delimiter=',')
        for row in reader:
            s += sum(list(map(float, row)))
    return s


def get_profile_files_per_chrom(fm):
    fnames = [os.path.join(fm.profilefolder, fname) for fname in os.listdir(fm.profilefolder)]
    for chrom in range(1, 23):
        pattern = re.compile(r"\.{}\.".format(chrom))
        chrom_files = list(filter(lambda x: re.search(pattern, x), fnames))
        fwd = [f for f in chrom_files if f.endswith('.fwd')]
        rev = [f for f in chrom_files if f.endswith('.rev')]
        ifwd = [f for f in chrom_files if f.endswith('.ifwd')]
        irev = [f for f in chrom_files if f.endswith('.irev')]
        yield chrom, fwd, rev, ifwd, irev


def do_it(fwd):
    d = {}
    pattern = re.compile(r"\w*\.bam")
    for fname in fwd:
        name = re.match(pattern, os.path.basename(fname)).group()
        d[name] = sum_profile_file(fname)
    return d


def create_upstream(fm):
    res = defaultdict(dict)
    for chrom, fwd, rev, _, _ in get_profile_files_per_chrom(fm):
        d1 = do_it(fwd)
        d2 = do_it(rev)
        c = Counter(d1)
        c.update(Counter(d2))
        for k, v in dict(c).items():
            res[k].update({chrom: v})
    return res


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
    res = create_upstream(f)
    for k, v in res.items():
        print(k, v)
