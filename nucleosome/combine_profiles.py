import os
import re
import csv
from collections import Counter, defaultdict
import matplotlib.pyplot as plt

from log_setup import setup_logger
logger = setup_logger(__name__, 'logs/nucleosome.log')


def map_to_int(arr):
    return list(map(int, arr))


def get_profile_file(proffile):
    logger.debug(f'Get profile file {proffile}')
    l = []
    with open(proffile, 'r') as in_:
        reader = list(csv.reader(in_, delimiter=','))
        assert len(reader) == 1
        for row in reader:
            l = list(map(float, row))
            l = map_to_int(l)
    return l


def get_profile_files_per_chrom(fnames, stream='up'):
    for chrom in range(1, 23):
        pattern = re.compile(r"\.{}\.".format(chrom))
        chrom_files = list(filter(lambda x: re.search(pattern, x), fnames))
        logger.debug(f'Got {len(chrom_files)} on chrom {chrom}')
        if stream == 'up':
            fwd = [f for f in chrom_files if f.endswith('.fwd')]
            rev = [f for f in chrom_files if f.endswith('.rev')]
            logger.debug(f'Yielding {len(fwd)} fwds and {len(rev)} revs')
            yield chrom, fwd, rev
        else:
            ifwd = [f for f in chrom_files if f.endswith('.ifwd')]
            irev = [f for f in chrom_files if f.endswith('.irev')]
            logger.debug(f'Yielding {len(ifwd)} ifwds and {len(irev)} irevs')
            yield chrom, ifwd, irev


def combine_profile_files(arr):
    d = {}
    name = ''
    pattern = re.compile(r"\w*\.bam")      # FIXME have to move pattern outside this
    pattern_mar = re.compile(r'P201.-[0-9]*_dedup\.bam')
    for fname in arr:
        try:
            name = re.match(pattern, os.path.basename(fname)).group()
        except AttributeError:
            name = re.match(pattern_mar, os.path.basename(fname)).group()
        except:
            logger.error(f'regexp problem for: {fname}')
            exit(f'AttributeError on regexp {fname}: check logs.')
        if name:
            d[name] = get_profile_file(fname)
        else:
            exit(f'AttributeError on regexp {fname}: check logs.')

    return d


def create_upstream(fnames):
    upstream = defaultdict(list)
    tmp = defaultdict(list)
    for chrom, fwd, rev in get_profile_files_per_chrom(fnames, stream='up'):
        if len(fwd) == len(rev) == 0:
            logger.error(f'No file found in {fnames} for chrom {chrom}')
            continue
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


def create_streams(fm, training=True):
    if training:
        fnames = [os.path.join(fm.profilefolder, fname) for fname in os.listdir(fm.profilefolder)]
    else:
        fnames = [os.path.join(fm.testprofilefolder, fname) for fname in os.listdir(fm.testprofilefolder)]
    streams = {}
    upstream = create_upstream(fnames)
    downstream = create_downstream(fnames)

    for sample in upstream.keys():
        up = upstream[sample]
        up.reverse()
        down = downstream[sample]
        streams[sample] = up[:-1] + down

    return streams


def save_streams_to_file(fm, outfile, training=True):
    streams = create_streams(fm, training)
    with open(outfile, 'w') as out:
        for sample, stream in streams.items():
            out.write('{},{}\n'.format(sample, ','.join(map(str, stream))))
    img_profile = outfile + '.png'
    plot_streams(streams, img_profile)
    return outfile, img_profile


def plot_streams(streams, figname):
    x = range(293)  # 147 * 2 - 1  # we got rid of center value
    for sample, stream in streams.items():
        y = stream
        plt.plot(x, y, label=sample)
    plt.legend(loc='upper center')
    plt.savefig(figname)
    # plt.show()


if __name__ == "__main__":
    logger = setup_logger('main', None)

    fnames = ['/data/tempff/comparison/new/p0.22.rev', '/data/tempff/comparison/new/p0.22.fwd',
              '/data/tempff/comparison/new/p0.22.irev', '/data/tempff/comparison/new/p0.22.ifwd']

    outfile ='/data/tempff/comparison/new/p0.22.profile'

    streams = {}
    upstream = create_upstream(fnames)
    downstream = create_downstream(fnames)

    for sample in upstream.keys():
        up = upstream[sample]
        up.reverse()
        down = downstream[sample]
        streams[sample] = up[:-1] + down

    with open(outfile, 'w') as out:
        for sample, stream in streams.items():
            out.write('{},{}\n'.format(sample, ','.join(map(str, stream))))

