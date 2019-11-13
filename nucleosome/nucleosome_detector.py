import os
import logging
import operator
import time

logger = logging.getLogger(__name__)

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


def create_nucleosome_file(chrom, mergefile, outfile):
    if os.path.isfile(outfile):
        logger.info("Nucleosome file for {} = {} already there. Skipping..".format(mergefile, outfile))
        return outfile
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
        if count % 1000000 == 0:
            logger.debug("[p:{} f:{}], line {} ({}). time = {}".format(os.getpid(), mergefile, count, '{:.1%}'
                                                                       .format(count/length), time.time() - s))
            s = time.time()

    allNucl = flush(curArea, lastPos, allNucl)

    # Dump results to a file
    with open(outfile, "w") as out_:
        for nucl in allNucl:
            out_.write("\t".join([str(x) for x in nucl]) + "\n")
    logger.info('Nucleosome saved for chrom {} in {}'.format(chrom, outfile))
    return outfile


