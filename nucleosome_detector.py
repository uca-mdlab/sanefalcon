import configparser
import os
import logging
import operator
import re
import collections
import concurrent.futures
import threading

MAX_THREAD_NUMBER = 15

logger = logging.getLogger('nucleosome_detector')
logging.basicConfig(format='%(asctime)s - %(levelname)s:%(message)s',
                    level=logging.DEBUG,
                    filename='nucleosome.log',
                    filemode='w')

chromosomes = range(1, 23)

areaSize = 73  # 60
sideSize = 20  # 25
padding = areaSize + sideSize
outerLen = 2 * sideSize  # len(left)+len(right)
innerLen = areaSize * 2 + 1  # len(window)


def find_merge_anti_files(trainfolder):
    merge_files = []
    anti_files = []
    root_merge_files = []
    pattern_name = re.compile("merge.\d{1,2}")
    anti_pattern_name = re.compile("anti.\d{1,2}")
    pattern_subdir = re.compile("/[a-z]/")
    for root, subdirs, files in os.walk(trainfolder):
        for fname in files:
            filename = os.path.join(root, fname)
            if re.match(pattern_name, fname) and re.search(pattern_subdir, filename):
                merge_files.append(filename)
            elif re.match(pattern_name, fname) and not re.search(pattern_subdir, filename):
                root_merge_files.append(filename)
            if re.match(anti_pattern_name, fname) and re.search(pattern_subdir, filename):
                anti_files.append(filename)
    return merge_files, anti_files, root_merge_files


def flush(area, endPoint, allNucl):
    # Sliding window settings, area is half the nucleosome size
    bins = [0] * padding + area + [0] * padding
    areaLen = len(area)
    startPoint = endPoint - areaLen + 1  # TODO check +1?

    # Let's move out
    bpScores = []
    extra = []
    peaks = []

    # Slide  a window of 20\147/20 over the region and score positions
    for i in range(padding, areaLen + padding):
        window = bins[i - areaSize: i + areaSize + 1]

        leftVal = sum(bins[i - areaSize - sideSize: i - areaSize])  # calcul de 0 à 20
        rightVal = sum(bins[i + areaSize + 1: i + areaSize + sideSize + 1])  # calcul de 167 à 187
        innerVal = sum(window)  # calcul de 20 à 167
        outerVal = leftVal + rightVal  # min(leftVal,rightVal)
        # if leftVal > 25 and rightVal > 25:
        # Solve zero devision and score this pos twice as high
        #	innerVal=max(innerVal,0.5)
        score = ((outerVal + 1) / float(outerLen)) / ((innerVal + 1) / float(innerLen))
        bpScores.append(score)
        extra.append([leftVal, innerVal, rightVal])
        # else:
        #	bpScores.append(0)

    # Pick best scoring positions in the region until no positive values exist
    def findCenters(start, end):
        if end - start < 1:
            return []
        maxIndex, bpMax = max(enumerate(bpScores[int(start):int(end)]), key=operator.itemgetter(1))
        if bpMax > 1:

            # If the max is part of a flat line take it's center position instead
            tmpIndex = maxIndex
            while tmpIndex < (end - start - 1) and bpScores[tmpIndex + 1] == bpMax:
                tmpIndex += 1
            maxIndex = (maxIndex + tmpIndex) / 2

            left = start + maxIndex - innerLen
            right = start + maxIndex + innerLen + 1
            leftList = findCenters(start, left)
            rightList = findCenters(right, end)
            thisList = [start + maxIndex, bpMax]
            return leftList + [thisList] + rightList
        return []

    newCenters = findCenters(0, areaLen)
    if newCenters != [] and newCenters[-1] != [] and newCenters[-1][0] > len(extra):
        print(newCenters[-1], len(extra), len(bpScores), bpScores[-5:])
    allNucl.extend([[x[0] + startPoint] + x[1:] + extra[int(x[0])] for x in newCenters])

    # if len(sys.argv) > 3:
    #     plotRegion(startPoint, endPoint, area, bpScores, newCenters)
    return allNucl


def assemble_runs(trainfolder, files, file_template, subdirs=None):
    runs = []
    tmp = collections.defaultdict(list)
    print(tmp,"tmp1")
    for f in files:
        tmp[f.split('.')[1]].append(f)
    print(tmp,"tmp2")
    print(subdirs,'SUBDIRS')

    if subdirs:
        for folder in subdirs:
            if len(folder) == 1:
                print(folder,'FOLDERS')
                pattern = re.compile('({})'.format(folder))  # matching "folder"
                reduced_tmp = {k: list(filter(lambda x: re.match(pattern, x), v)) for k, v in tmp.items()}
                fname_stub = os.path.join(folder, file_template)
                print(reduced_tmp)
                run = [(folder, k, v[0], fname_stub) for k, v in reduced_tmp.items()]
                runs.extend(run)
            else:
                logger.error("there is a problem with folder {} it's longer than 1".format(folder))
    else:
        fname_stub = os.path.join(trainfolder, file_template)
        run = [(trainfolder, k, v[0], fname_stub) for k, v in tmp.items()]
        runs.extend(run)
    return runs


def _create_nucleosome_file(folder, chrom, mergefile, fname):
    print('starting {} on file {}'.format(os.getpid(), mergefile))
    logger.info("Creating nucleosome file {} for {}".format(fname, mergefile))
    curArea = [0]
    lastPos = 0
    maxDist = 190  # A little over our sliding window size
    allNucl = []
    with open(mergefile, 'r') as infile:
        logger.debug('working on {}'.format(infile.name))
        positions = list(map(int, infile.readlines()))
        length = len(positions)
    s = time.time()
    count = 0
    for position in positions:
        count += 1
        # position = int(line.strip())
        distance = position - lastPos
        if distance > maxDist:
            allNucl = flush(curArea, lastPos, allNucl)
            curArea = [1]
        # Add read to current region
        else:
            curArea += [0 for x in range(distance)]
            # print(curArea)
            curArea[-1] += 1
        lastPos = position
        if count % 100000 == 0:
            logger.debug("[p:{} f:{}], line {} ({}). time = {}".format(os.getpid(), mergefile, count, '{:.1%}'
                                                                       .format(count/length), time.time() - s))
            s = time.time()

    allNucl = flush(curArea, lastPos, allNucl)

    # Dump results to a file
    with open(os.path.join(folder, fname + ".{}".format(chrom)), "w") as output_file:
        for nucl in allNucl:
            output_file.write("\t".join([str(x) for x in nucl]) + "\n")
    logger.info('Nucleosome saved for chrom {} in {}'.format(chrom, output_file.name))


def create_nucleosome_files(trainfolder, nucl_file_template,anti_file_template):
    subdirs = [f.path for f in os.scandir(trainfolder) if f.is_dir()]
    merge_files, anti_files, root_merge_files = find_merge_anti_files(trainfolder)
    runs = assemble_runs(trainfolder, merge_files, nucl_file_template, subdirs)
    runs.extend(assemble_runs(trainfolder, anti_files, anti_file_template, subdirs))
    runs.extend(assemble_runs(trainfolder, root_merge_files, nucl_file_template))
    with concurrent.futures.ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        jobs = {executor.submit(_create_nucleosome_file, *run): run for run in runs}
        for job in concurrent.futures.as_completed(jobs):
            try:
                _ = job.result()
            except Exception as ex:
                logger.error('Future Exception {}'.format(ex.__cause__))

if __name__ == "__main__":
    import time
    conf_file = 'sanefalcon.conf'
    config = configparser.ConfigParser()
    config.read(conf_file)

    trainfolder = config['default']['trainfolder']
    nucl_file_template = config['default']['nucltemplate']
    anti_file_template = nucl_file_template + '_anti'

    subdirs = [f.path for f in os.scandir(trainfolder) if f.is_dir()]
    merge_files, anti_files, root_merge_files = find_merge_anti_files(trainfolder)
    runs = []
    #runs = assemble_runs(trainfolder, merge_files, nucl_file_template, subdirs)
    # runs.extend(assemble_runs(trainfolder, anti_files, anti_file_template, subdirs))
    runs.extend(assemble_runs(trainfolder, root_merge_files, nucl_file_template))
    # print(len(runs))
    # exit()
    # s = time.time()
    # create_nucleosome_file(*tup)
    # print('done:', time.time() - s)
    # exit()
    s = time.time()
    with concurrent.futures.ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        jobs = {executor.submit(_create_nucleosome_file, *run): run for run in runs}
        for job in concurrent.futures.as_completed(jobs):
            try:
                _ = job.result()
            except Exception as ex:
                logger.error('Future Exception {}'.format(ex.__cause__))
    print('End: ', time.time() - s)
    # create_nucl_files(trainfolder, nucl_file_template, 'merge')
    # anti_template = nucl_file_template + '_anti'
    # create_nucl_files(trainfolder, anti_template, 'anti')
