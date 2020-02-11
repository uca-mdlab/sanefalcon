import pickle
from collections import defaultdict, Counter
import concurrent.futures
from itertools import islice
import os
import time
import operator
from log_setup import setup_logger


logger = setup_logger('nucl.detect', 'logs/nucldetect.log')

count = None
scores = defaultdict(list)


def check_pos(readstart):
    global count
    global scores
    low = readstart - 93
    high = readstart + 73

    left = {k: count[k] for k in range(low, low + 20 + 1)}
    right = {k: count[k] for k in range(high, high + 20 + 1)}
    centers = {k: count[k] for k in range(low + 20, high + 1)}

    linker = sum(left.values()) + sum(right.values())
    center = sum(centers.values())
    if center > 0:
        scores[readstart] = linker / center
    else:
        scores[readstart] = 0
    return


def work(chunk):
    for rs in chunk:
        check_pos(rs)


def prepare_chunk(lis, size):
    it = iter(lis)
    return iter(lambda: tuple(islice(it, size)), ())


def extract_scores(mergefile):
    global count
    global scores
    s0 = time.time()
    # scores = defaultdict(list)
    scorefile = os.path.basename(mergefile) + '.scores'
    if os.path.isfile(scorefile):
        scores = pickle.load(open(scorefile, 'rb'))
    else:
        positions = pickle.load(open(mergefile, 'rb'))
        length = len(positions)
        count = Counter(positions)
        l = len(count)
        pool_size = 50
        logger.debug(f'{l} positions to check: ({int(l/pool_size)}) rounds in {mergefile}')
        i = 0
        s = s0
        for chunk in prepare_chunk(list(count.keys()), pool_size):
            with concurrent.futures.ThreadPoolExecutor(max_workers=pool_size) as executor:
                jobs = {}
                sub = executor.submit(work, chunk)
                jobs[sub] = chunk

                for job in concurrent.futures.as_completed(jobs):
                    _ = job.result()
                    # rs, score = job.result()
                    # scores[rs] = score
            if i % 10000 == 0:
                elaps = time.time() - s
                logger.debug("chunk {} (pos:{}) time: {}".format(i, chunk[0], elaps))
                s = time.time()
            i += 1

        if scores:
            pickle.dump(scores, open(scorefile, 'wb'))

    logger.info(f'{mergefile} scores extracted in {round(time.time() - s0)} sec. ')
    return scores


def apply_mask(readstart):
    global scores
    low = readstart - 147
    high = readstart + 147
    for k in range(low, high + 1):
        if k in scores:
            scores.pop(k)


def find_nucleosomes(mergefile):
    global scores
    scores = extract_scores(mergefile)
    nucleosomes = defaultdict(float)

    s0 = time.time()
    s = s0
    i = 0
    sort_list = [(k, v) for k, v in sorted(scores.items(), key=operator.itemgetter(1), reverse=True) if v >= 1]
    while sort_list:
        bp, max_score = sort_list.pop(0)
        nucleosomes[bp] = max_score
        to_remove = list(filter(lambda x: bp - 147 <= x[0] <= bp + 147, sort_list))
        indexes = sorted([sort_list.index(x) for x in to_remove], reverse=True)
        for index in indexes:
            del sort_list[index]
        if i % 10000 == 0:
            logger.debug("max score: {}; len(scores): {}; time: {}".format(max_score, len(sort_list), round(time.time() - s)))
            s = time.time()
        i += 1
    """
    max_score = 2
    while max_score >= 1:
        bp, max_score = max(scores.items(), key=operator.itemgetter(1))
        # bp, max_score = sorted(scores.items(), key=lambda item: item[1])[-1]  # last element
        nucleosomes[bp] = max_score
        if i % 1000 == 0:
            print(max_score, len(scores), time.time() - s)
            s = time.time()
        #scores = mask(bp, scores)
        apply_mask(bp)
        i += 1
    """
    logger.info(f'{mergefile} nucleosome computed in {round(time.time() - s0)} sec. ')
    return nucleosomes


def create_nucleosome_file(chrom, mergefile, outfile):
    if not os.path.isfile(outfile):
        nucleosomes = find_nucleosomes(mergefile)
        pickle.dump(nucleosomes, open(outfile, 'wb'))
    else:
        print(f'{outfile} already there. Skipping')
    return outfile


if __name__ == "__main__":
    mergefile = '/data/tempff/anti.22'
    s = time.time()
    nucleosomes = find_nucleosomes(mergefile)
    pickle.dump(nucleosomes, open('nucl.anti.22', 'wb'))
    print(time.time() - s)