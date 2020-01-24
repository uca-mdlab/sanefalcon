from collections import defaultdict
import os
import re

from log_setup import setup_logger

logger = setup_logger(__name__, 'logs/sanefalcon.log')

chromosomes = range(1, 23)


def merge_streams(profile_dir, samplename, ext, rev=False):
    logger.debug(f'Merge stream on {profile_dir} for {samplename} (ext = {ext}')
    dic = defaultdict(float)
    files = [os.path.join(profile_dir, f) for f in os.listdir(profile_dir) if
             re.match(samplename, f) and (f.endswith(ext[0]) or f.endswith(ext[1]))]
    logger.debug(f'Merging data from {len(files)} files')
    logger.debug(f'{files}')
    for f in files:
        with open(f, 'r') as infile:
            arr = infile.readlines()
            assert len(arr) == 1
            chrom_profile = list(map(float, arr[0].strip().split(',')))
            for i, v in enumerate(chrom_profile):
                dic[i] += v
            # logger.debug(f'len(dic) = {len(dic)}')

    merged = list(dic.values())
    if rev:
        merged.reverse()
    return merged


def combine_profiles(samplename, profile_dir):
    ext1 = ['.rev', '.fwd']
    ext2 = ['.irev', '.ifwd']
    profile = merge_streams(profile_dir, samplename, ext1, rev=True) + merge_streams(profile_dir, samplename, ext2)
    normalized = [x / sum(profile) for x in profile]
    return normalized


def parse_model_file(model_file):
    with open(model_file, 'r') as infile:
        arr = infile.readlines()
    assert len(arr) == 2
    corr_string, scal_string = arr
    correlations = [float(x) for x in corr_string.strip().split()]
    scalars = [float(x) for x in scal_string.strip().split()]
    return correlations, scalars


def predict(model_file, sample, profile_dir):
    nucleosome_profile = combine_profiles(sample, profile_dir)
    correlations, scalars = parse_model_file(model_file)

    tot = sum(val * corr for val, corr in zip(nucleosome_profile, correlations))
    ff = scalars[0] * tot + scalars[1]
    return ff, nucleosome_profile


def compute_ff(modelname, profile_dir):
    result = {}
    samplenames = list({f.split('.bam')[0] + '.bam' for f in os.listdir(profile_dir)})
    logger.debug(f'Found {len(samplenames)} sample names for testing')
    for samplename in samplenames:
        ff, nuclprofile = predict(modelname, samplename, profile_dir)
        with open(os.path.join(profile_dir, samplename + '.ff'), 'w') as out:
            out.write('Fetal Fraction: {}\n'.format(ff))
            out.write('Nucleosome Profile: {}'.format(nuclprofile))
        result[samplename] = ff
        logger.debug(f"ff({samplename}) = {ff}")
    return result


