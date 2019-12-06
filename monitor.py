import os
import re

from collections import defaultdict

from manager.file_manager import FileManager
from manager.utils import Utils
from manager.rsp import RspBuilder
from merger.merge import merge, merge_subs, merge_anti_subs
from nucleosome.tracker import Tracker
from nucleosome.profiler import Profiler

from log_setup import setup_logger

logger = setup_logger(__name__, 'logs/sanefalcon.log')


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
    for samplename, list_ in dic.items():
        fwds = [fname for fname in list_ if fname.endswith('.fwd')]
        revs = [fname for fname in list_ if fname.endswith('.rev')]
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
    anti = defaultdict(list)
    anti.update((k, []) for k in mergeddic.keys())
    for path, pathfiles in mergeddic.items():
        for antipath in anti.keys():
            if antipath is not path:
                anti[antipath].extend(pathfiles)

    antisubs = []
    for path, files in anti.items():
        for chrom in range(1, 23):
            pattern = re.compile(r'\.{}$'.format(chrom))
            per_chrom = [f for f in files if re.search(pattern, f)]
            antisubs.append((path, chrom, per_chrom))

    anti_sub_files = merge_anti_subs(antisubs)
    return anti_sub_files


def prepare_and_merge(fm, rsb, config):
    bamlist = Utils.readfile(config['training']['bamlist'])
    batchsize = int(config['training']['batchsize'])
    batches = fm.prepare_train_folder(bamlist, batchsize)
    logger.info('Batches prepared')
    rspfiles = rsb.prepare_samples(fm.datafolder, fm.rspfolder)
    logger.info('Samples prepared')
    mapping = get_rsp_batches_mapping(batches, rspfiles)
    merged = launch_merge(mapping)
    # merged_subs = launch_merge_subs(merged, fm.trainfolder)
    merged_anti = launch_merge_anti_subs(merged)
    logger.info('Merge terminated')
    return mapping, merged, merged_anti


def create_nucleosome_profiles(fm, mapping, training=True):
    t = Tracker(fm)
    t.create_tracks(training)
    p = Profiler(fm, t)
    p.compute_profiles(mapping, training)
    return p.combine(training)


def training(config):
    logger.info('Starting...')
    f = FileManager(config)
    f.check_paths()
    rs = RspBuilder(config)
    mapping, merged, anti = prepare_and_merge(f, rs, config)
    nucleosome_file, img_file = create_nucleosome_profiles(f, mapping)
    logger.info('Nucleosome file created: {}'.format(nucleosome_file))
    return nucleosome_file


# Testing phase
def filter_out_rsp_files(fm):
    train_names = [os.path.basename(f) for f in Utils.readfile(fm.config['training']['bamlist'])]
    logger.debug('Found {} names in training'.format(len(train_names)))
    pattern = re.compile(r'\.\d{1,2}\.')
    testing_rsp = defaultdict(list)
    rsps = []
    for f in fm.get_read_start_positions_data():
        basename = os.path.basename(f)
        try:
            token = re.search(pattern, basename).group()
            name = basename.split(token)[0]
            if name not in train_names:
                testing_rsp[name].append(f)
                rsps.append(f)
        except AttributeError:
            pass

    batch = fm.create_fake_batch_for_testing(list(testing_rsp.keys()))
    mapping = get_rsp_batches_mapping(batch, rsps)
    fwd, rev = get_fwd_rev_files(testing_rsp)
    merge_file_list = defaultdict(dict)
    merge_file_list[fm.testfolder] = {'fwd': fwd, 'rev': rev}

    return mapping, merge_file_list


def testing(config):
    fm = FileManager(config)
    mapping, rsp_to_merge = filter_out_rsp_files(fm)
    merged = merge(rsp_to_merge)
    logger.info('Merged testing: {}'.format(merged))
    nucleosome_file, img_file = create_nucleosome_profiles(fm, mapping, training=False)
    logger.info('Nucleosome file created: {}'.format(nucleosome_file))
    return nucleosome_file




