# This is the replacement for the following:
# merge.sh, mergeSubs.sh, mergeAntiSub.sh

import configparser
import concurrent.futures
import os
import logging
import re


chromosomes = range(1, 23)


def sort_and_write(data, outfile):
    data = list(set(data))  # remove duplicates and sort
    data.sort()
    to_be_written = map(lambda x: str(x) + "\n", data)
    with open(outfile, 'w') as out:
        out.writelines(to_be_written)


def find_merge_files_in_subdirectories(trainfolder):
    merge_files = []
    pattern_name = re.compile("merge.\d{1,2}")
    pattern_subdir = re.compile("/[a-z]/")
    for root, subdirs, files in os.walk(trainfolder):
        for fname in files:
            filename = os.path.join(root, fname)
            if re.match(pattern_name, fname) and re.search(pattern_subdir, filename):
                merge_files.append(filename)
    return merge_files


def prepare_file_lists(trainfolder):
    """

    :param trainfolder: sanefalcontrain
    :return: {chr1: {a: [sanefalcontrain/a/chr1.fwd,sanefalcontrain/a/chr1.rev], b: [...]},
              chr2: {a: [sanefalcontrain/a/chr2.fwd,sanefalcontrain/a/chr2.rev], b: [...]},
              ...
             }
    """

    files_dic = dict.fromkeys(chromosomes, dict())

    pattern = re.compile("\.bam\.\d{1,2}\.start")
    for root, subdirs, files in os.walk(trainfolder):
        for fname in files:
            if fname.endswith('.fwd') or fname.endswith('.rev'):
                filename = os.path.join(root, fname)
                match = re.search(pattern, fname)
                chrom = int(match.group(0).split(".")[2])
                subdir = os.path.dirname(filename)
                if subdir in files_dic[chrom]:
                    files_dic[chrom][subdir].append(filename)
                else:
                    files_dic[chrom] = {subdir: [filename]}

    return files_dic


def merge(files_dic):
    """
    input: [sanefalcontrain/a/sample1.start.fwd, sanefalcontrain/a/sample1.start.rev]
    output: sanefalcontrain/a/merge.chr1
    :param trainfolder:
    :return:
    """
    for chrom, dic in files_dic.items():
        for dir, files in dic.items():
            data = []
            for f in files:
                data.extend([int(line.strip()) for line in open(f, 'r')])
            outfile = os.path.join(dir, "merge.{}".format(chrom))
            sort_and_write(data, outfile)


def merge_subs(trainfolder, files_dic):
    """
    input: [sanefalcontrain/a/merge.chr1, sanefalcontrain/b/merge.chr1, sanefalcontrain/c/merge.chr1, ...]
    output: sanefalcontrain/merge.chr1
    :param trainfolder:
    :return:
    """
    for chrom in chromosomes:
        files_to_merge = []
        for dir, _ in files_dic[chrom].items():
            subdir_merge_file = os.path.join(dir, "merge.{}".format(chrom))
            if not os.path.isfile(subdir_merge_file):
                exit("{} not found".format(subdir_merge_file))
            files_to_merge.append(subdir_merge_file)

        data = []
        for mergefile in files_to_merge:
            data.extend([int(line.strip()) for line in open(mergefile, 'r')])
        outfile = os.path.join(trainfolder, "merge.{}".format(chrom))
        sort_and_write(data, outfile)


def merge_anti_subs(trainfolder):
    """
    input: [sanefalcontrain/b/merge.chr1, sanefalcontrain/c/merge.chr1, ...]
    output: sanefalcontrain/a/anti.chr1
    :param trainfolder:
    :return:
    """
    all_merge_files = find_merge_files_in_subdirectories(trainfolder)
    subdirs = [subdir.path for subdir in os.scandir(trainfolder) if subdir.is_dir()]
    done = dict.fromkeys(subdirs, False)
    for subdir in subdirs:
        if not done[subdir]:
            current = subdir
            merge_files = [x for x in all_merge_files if not re.match(current, x)]
            for chrom in chromosomes:
                data = []
                files_to_merge = [f for f in merge_files if f.endswith(".{}".format(chrom))]
                for mergefile in files_to_merge:
                    data.extend([int(line.strip()) for line in open(mergefile, 'r')])
                if len(data) > 0:
                    outfile = os.path.join(subdir, "anti.{}".format(chrom))
                    sort_and_write(data, outfile)
                done[subdir] = True

if __name__ == "__main__":
    conf_file = 'sanefalcon.conf'
    config = configparser.ConfigParser()
    config.read(conf_file)

    samtools = config['default']['samtools']
    datafolder = config['default']['datafolder']
    trainfolder = config['default']['trainfolder']
    nucleosomefolder = config['default']['nucleosomefolder']

    dic = prepare_file_lists(trainfolder)
    merge(dic)
    print("merge done")
    merge_subs(trainfolder, dic)
    print("merge_subs done")
    merge_anti_subs(trainfolder)
    print("merge_anti_subs done")
