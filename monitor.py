import configparser
import logging
import os
import re

from collections import defaultdict

from manager.file_manager import FileManager
from manager.utils import Utils
from manager.rsp import RspBuilder
from merger.merge import merge, merge_subs

logging.basicConfig(format='%(asctime)s - %(name)-21s - %(levelname)-8s: %(message)s',
                    filename='sanefalcon.log', filemode='w', level=logging.DEBUG)

logger = logging.getLogger('monitor')


def get_rsp_batches_mapping(batches, rspfiles):
    mapping = defaultdict(dict)
    for subdir, list_of_fnames in batches.items():
        bam_names = [os.path.basename(f) for f in list_of_fnames]
        for bam in bam_names:
            rsp = list(filter(lambda x: re.search(bam, x), rspfiles))
            mapping[subdir][bam] = rsp

    logger.info('Mapping complete.')
    return mapping


def get_fwd_rev_files(dic):
    fwd = defaultdict(list)
    rev = defaultdict(list)
    fwds = []
    revs = []
    for samplename, list_ in dic.items():
        fwds.extend([fname for fname in list_ if fname.endswith('.fwd')])
        revs.extend([fname for fname in list_ if fname.endswith('.rev')])
        for chrom in range(1, 23):
            pattern = re.compile(r'\.{}\.'.format(chrom))
            fwd[chrom].extend(list(filter(lambda x: re.search(pattern, x), fwds)))
            rev[chrom].extend(list(filter(lambda x: re.search(pattern, x), revs)))
    return fwd, rev


def launch_merge(mapping):
    merge_file_list = defaultdict(dict)
    for subdir, dic in mapping.items():
        fwd, rev = get_fwd_rev_files(dic)
        merge_file_list[subdir] = {'fwd': fwd, 'rev': rev}
    res = merge(merge_file_list)
    return res


def launch_merge_subs(mergeddic, trainfolder):
    merged = defaultdict(list)
    for subdir, files in mergeddic.items():
        for f in files:
            merged[os.path.basename(f)].append(f)

    merge_subs_files = merge_subs(merged, trainfolder)
    return merge_subs_files


def launch_merge_anti_subs(mergeddic):
    for subdir, files in mergeddic.items():
        print(subdir, files)




if __name__ == '__main__':
    config = configparser.ConfigParser()
    config.read('local.conf')
    logger.info('Starting...')
    f = FileManager(config)
    f.check_paths()
    rs = RspBuilder(config)
    # print(f.get_data())
    # print(f.get_data(mask='.bam$'))
    bamlist = Utils.readfile(config['training']['bamlist'])
    batchsize = int(config['training']['batchsize'])
    batches = f.prepare_train_folder(bamlist, batchsize)
    rspfiles = rs.prepare_samples(f.datafolder, f.rspfolder)
    mapping = get_rsp_batches_mapping(batches, rspfiles)
    merged = launch_merge(mapping)
    merged_subs = launch_merge_subs(merged, f.trainfolder)
    merged_anti = launch_merge_anti_subs(merged)
