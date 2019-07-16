import os
import logging
import operator
import re
import collections
import concurrent.futures
import time

MAX_THREAD_NUMBER = 10

logger = logging.getLogger('nucleosome_detector')
logging.basicConfig(format='%(asctime)s - %(levelname)s:%(message)s',
                    level=logging.DEBUG,
                    filename='sanefalcon.log',
                    filemode='w')

chromosomes = range(1, 23)

areaSize = 73  # 60
sideSize = 20  # 25
padding = areaSize + sideSize
outerLen = 2 * sideSize  # len(left)+len(right)
innerLen = areaSize * 2 + 1  # len(window)


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
    for f in files:
        tmp[f.split('.')[1]].append(f)

    if subdirs:
        for folder in subdirs:
            logger.debug('Assembling runs on {}'.format(folder))
            if len(os.path.split(folder)[1]) == 1:
                pattern = re.compile('({})'.format(folder))  # matching "folder"
                reduced_tmp = {k: list(filter(lambda x: re.match(pattern, x), v)) for k, v in tmp.items()}
                fname_stub = os.path.join(folder, file_template)
                run = [(folder, k, v[0], fname_stub) for k, v in reduced_tmp.items()]
                runs.extend(run)
            else:
                logger.error("there is a problem with folder {} it's longer than 1".format(folder))  # useless FIXME
    else:
        fname_stub = os.path.join(trainfolder, file_template)
        run = [(trainfolder, k, v[0], fname_stub) for k, v in tmp.items()]
        runs.extend(run)
    return runs


def _create_nucleosome_file(folder, chrom, mergefile, fname):
    outfile = os.path.join(folder, fname + ".{}".format(chrom))
    if os.path.isfile(outfile):
        logger.info("Nucleosome file for {} = {} already there. Skipping..".format(mergefile, outfile))
        return
    logger.info("p:{} creating nucleosome file for {} = {}".format(os.getpid(), mergefile, outfile))
    curArea = [0]
    lastPos = 0
    maxDist = 190  # A little over our sliding window size
    allNucl = []
    with open(mergefile, 'r') as infile:
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


def create_nucleosome_files(fm, training=True):
    subdirs = list(fm.find_all_manips_per_subfolder().keys())

    merge_files, anti_files, root_merge_files = fm.find_merge_anti_files()

    runs = assemble_runs(fm.trainfolder, anti_files, fm.anti_file_template, subdirs)
    if not training:
        logger.info('Training set to FALSE. Extending nucleosome runs')
        runs.extend(assemble_runs(fm.trainfolder, merge_files, fm.nucl_file_template, subdirs))
        runs.extend(assemble_runs(fm.trainfolder, root_merge_files, fm.nucl_file_template))

    tot = len(runs)
    logger.info('Creating nucleosomes: {} runs to be launched'.format(tot))
    with concurrent.futures.ThreadPoolExecutor(max_workers=os.cpu_count() * 2) as executor:
        jobs = {executor.submit(_create_nucleosome_file, *run): run for run in runs}
        for job in concurrent.futures.as_completed(jobs):
            try:
                _ = job.result()
            except Exception as ex:
                logger.error('Future Exception {}'.format(ex.__cause__))


if __name__ == "__main__":
    import configparser
    from file_manager import FileManager
    conf_file = 'tests/data/test.conf'
    config = configparser.ConfigParser()
    config.read(conf_file)

    fm = FileManager(config)
    create_nucleosome_files(fm)
