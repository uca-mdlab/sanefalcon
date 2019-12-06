import os
import concurrent.futures
import threading
from collections import defaultdict

from log_setup import setup_logger
logger = setup_logger(__name__, 'logs/nucleosome.log')

# parameters and functions taken from the original script
minSide = 25
minCenter = 100

areaSize = 73.  # 60
sideSize = 20.  # 25

filtOut = 0


def nuclFilt(prop):
    global filtOut
    if min(prop[0], prop[2])/sideSize > prop[1]/147.:  # and prop[1] > minCenter: #ex9
        return True
    filtOut += 1
    return False


def loadNucl(nuclLine):
    return nuclLine[0:2]


def cast_line_to_numbers(list_of_strings):
    def cast_to_number(string_literal):
        try:
            return int(string_literal)
        except ValueError:
            return float(string_literal)

    return list(map(cast_to_number, list_of_strings))


def load_data(nucl_ex_file, fwd_rev_file):
    thread_name = threading.current_thread().name
    logger.debug("(t-{}): Loading data from nucl_file: {}; fwd_rev_file: {}".format(thread_name, nucl_ex_file, fwd_rev_file))
    with open(nucl_ex_file, 'r') as infile:
        lines = list(map(cast_line_to_numbers, [x.split() for x in infile.readlines()]))
        peaks = []
        for l in lines:
            if nuclFilt(l[2:]):
                l[0] = int(l[0])
                peaks += [loadNucl(l)]

    with open(fwd_rev_file, 'r') as infile:
        reads = [int(line) for line in infile]
    logger.debug("{}: {} peaks, {} filtOut".format(os.path.basename(fwd_rev_file), len(peaks), filtOut))
    return lines, peaks, reads


def process_forward(peaks, reads, outfile):
    thread_name = threading.current_thread().name
    logger.debug('process_forward for {} (t-{})'.format(outfile, thread_name))
    logger.debug('process_forward (t-{}): peaks ([0]: {}, len: {})'.format(thread_name, peaks[0], len(peaks)))
    logger.debug('process_forward (t-{}): reads ([0]: {}, len: {})'.format(thread_name, reads[0], len(reads)))

    # before = len(peaks)  # just to be sure
    peaks.append([-1, 0])
    # logger.debug('process_forward ({}) (fwd append) id: {} (len += {}) [0] = {}; [-1] = {}'.format(thread_name, id(peaks), len(peaks) - before, peaks[0], peaks[-1]))

    maxDist = 147
    sumPeak = [0. for _ in range(maxDist)]
    read = reads[0]  # +shift
    j = 0
    nuclHit = []
    logger.debug('process_forward (t-{}): starting cycle on peaks'.format(thread_name))
    for i, peakPair in enumerate(peaks):
        # if i == 0 or i == len(peaks) - 1:
        #    logger.debug('process_forward {} (fwd enumerate id: {}) i: {} peakPair {}'.format(thread_name, id(peaks), i, peakPair))
        peak = peakPair[0]
        peakWeight = peakPair[1]
        thisPeak = [0. for x in range(maxDist)]
        while read <= peak:
            if read > peak - maxDist:
                thisPeak[peak - read] += 1
            j += 1
            if j >= len(reads):
                break
            read = reads[j]  # +shift

        thisSum = float(sum(thisPeak))
        if thisSum == 0:
            continue

        nuclHit.append(peak)
        # thisPeakNorm=[x/thisSum for x in thisPeak]
        # thisPeak=[x*peakWeight for x in thisPeak]
        sumPeak = [sumPeak[x] + thisPeak[x] for x in range(len(thisPeak))]

    logger.debug('process_forward (t-{}): completed cycle on peaks'.format(thread_name))
    logger.debug('process_forward (t-{}): len (nuclHit: {}, sumPeak: {})'.format(thread_name, len(nuclHit), len(sumPeak)))
    with open(outfile, 'w') as outPeaks:
        outPeaks.write(",".join([str(x) for x in sumPeak]))
    return sumPeak


def process_reverse(peaks, reads, outfile):
    thread_name = threading.current_thread().name
    if peaks[-1] == [-1, 0]:  # it happens to have already the [-1, 0] pair at the end. This is a hotfix. FIXME
        peaks = peaks[:-1]

    peaks.reverse()
    reads.reverse()
    before = len(peaks)  # just to be sure
    peaks.append([-1, 0])
    logger.debug('process_reverse {} (rev append) id: {} (len += {}) [0] = {}; [-1] = {}'.format(thread_name, id(peaks), len(peaks) - before, peaks[0], peaks[-1]))

    maxDist = 147
    sumPeak = [0. for _ in range(maxDist)]
    read = reads[0]  # +shift
    j = 0
    nuclHit = []

    for i, peakPair in enumerate(peaks):
        if i == 0 or i == len(peaks):
            logger.debug('process_reverse {} (rev enumerate id: {}) i: {} peakPair {}'.format(thread_name, id(peaks), i, peakPair))
        peak = peakPair[0]
        peakWeight = peakPair[1]
        thisPeak = [0. for x in range(maxDist)]
        while read >= peak:
            if read < peak + maxDist:
                thisPeak[read - peak] += 1
            j += 1
            if j >= len(reads):
                break
            read = reads[j]  # +shift
        thisSum = float(sum(thisPeak))
        if thisSum == 0:
            continue

        nuclHit.append(peak)
        # thisPeakNorm=[x/thisSum for x in thisPeak]
        # thisPeak=[x*peakWeight for x in thisPeak]
        sumPeak = [sumPeak[x] + thisPeak[x] for x in range(len(thisPeak))]
    nuclHit.reverse()

    with open(outfile, 'w') as outPeaks:
        outPeaks.write(",".join([str(x) for x in sumPeak]))
    return sumPeak


def run_forward(chrom, outdir, fwd_file, nucl_ex_file):
    output_stubname = os.path.join(outdir, os.path.basename(fwd_file) + '.{}'.format(chrom))
    sumPeakFwd = []
    sumPeakRev = []
    fwd_out = output_stubname + ".fwd"
    ifwd_out = output_stubname + ".ifwd"

    if not (os.path.isfile(fwd_out) and os.path.isfile(ifwd_out)):
        logger.debug("run_forward chrom {}. outdir={}, fwd_file={}, nucl_ex_file={}, "
                     "output_stubname={}".format(chrom, outdir, fwd_file, nucl_ex_file, output_stubname))
        lines, peaks, reads = load_data(nucl_ex_file, fwd_file)
        logger.debug(
            'run_forward chrom {}-{}, peaks id {}'.format(nucl_ex_file.split('.')[-1], os.path.basename(fwd_file),
                                                          id(peaks)))
        try:
            sumPeakFwd = process_forward(peaks, reads, fwd_out)  # -> .fwd
            logger.debug(
                "run_forward - chrom {} id {}: process_forward completed: sumPeaks [{}..{}]".format(chrom, id(peaks),
                                                                                                    sumPeakFwd[0],
                                                                                                    sumPeakFwd[-1]))
        except IndexError:
            logger.error("run_forward IndexError: .fwd: peaks {}, reads {}".format(len(peaks), len(reads)))

        try:
            sumPeakRev = process_reverse(peaks, reads, ifwd_out)  # -> .ifwd
            logger.debug("run_forward - chrom {} id {}: process_reverse completed: sumPeaks [{}..{}]".format(chrom, id(peaks), sumPeakRev[0], sumPeakRev[-1]))
        except IndexError:
            logger.error("IndexError: .ifwd: peaks {}, reads {}".format(len(peaks), len(reads)))
    else:
        logger.debug('Profile already stored: {}'.format(fwd_out))
        logger.debug('Profile already stored: {}'.format(ifwd_out))
    return (fwd_out, ifwd_out), len(sumPeakFwd), len(sumPeakRev)


def run_reverse(chrom, outdir, rev_file, nucl_ex_file):
    output_stubname = os.path.join(outdir, os.path.basename(rev_file) + '.{}'.format(chrom))
    irev_out = output_stubname + ".irev"
    rev_out = output_stubname + ".rev"
    sumPeakFwd = []
    sumPeakRev = []

    if not (os.path.isfile(irev_out) and os.path.isfile(rev_out)):
        lines, peaks, reads = load_data(nucl_ex_file, rev_file)
        logger.debug(
            'run_reverse chrom {}-{}. peaks id {}'.format(nucl_ex_file.split('.')[-1], os.path.basename(rev_file),
                                                          id(peaks)))

        try:
            sumPeakFwd = process_forward(peaks, reads, irev_out)  # -> .irev
            logger.debug("run_reverse - chrom {} id {}: process_forward completed: sumPeaks [{}..{}]".format(chrom, id(peaks), sumPeakFwd[0], sumPeakFwd[-1]))
        except IndexError:
            logger.error("IndexError: .irev: peaks {}, reads {}".format(len(peaks), len(reads)))

        try:
            sumPeakRev = process_reverse(peaks, reads, rev_out)  # -> .rev
            logger.debug("run_reverse - chrom {} id {}: process_reverse completed: sumPeaks [{}..{}]".format(chrom, id(peaks), sumPeakRev[0], sumPeakRev[-1]))
        except IndexError:
            logger.error("IndexError: .rev: peaks {}, reads {}".format(len(peaks), len(reads)))
    else:
        logger.debug('Profile already stored: {}'.format(irev_out))
        logger.debug('Profile already stored: {}'.format(rev_out))
    return (irev_out, rev_out), len(sumPeakFwd), len(sumPeakRev)


def process(chrom, fwdrevdic, outdir):
    MAX_FILES_PROCESSED = 15  # limit parallelism not to overload memory
    d = defaultdict(list)
    for subdir, list_ in fwdrevdic.items():
        for dic in list_:
            nucl_ex_file = dic['nucl_file']
            fwd_files = dic['fwd']
            fwd_files_left = len(fwd_files)
            fwd_files_iter = iter(fwd_files)

            rev_files = dic['rev']
            rev_files_left = len(rev_files)
            rev_files_iter = iter(rev_files)
            logger.debug('Starting multithread on {} with {}'.format(subdir, nucl_ex_file))
            with concurrent.futures.ThreadPoolExecutor(max_workers=5, thread_name_prefix='fwd') as executor:
                jobs = {}
                while fwd_files_left:
                    for this_fwd_file in fwd_files_iter:
                        job = executor.submit(run_forward, chrom, outdir, this_fwd_file, nucl_ex_file)
                        jobs[job] = this_fwd_file
                        if len(jobs) > MAX_FILES_PROCESSED:
                            break
                    for job in concurrent.futures.as_completed(jobs):
                        fwd_files_left -= 1
                        try:
                            l_tup = job.result()
                        except Exception as ex:
                            logger.error('[fwd] Exception {}'.format(ex.__cause__))
                            l_tup = (None,)
                        del jobs[job]
                        d[subdir].append(l_tup[0])
            logger.info('End of forward concurrent phase for chrom {}'.format(chrom))

            with concurrent.futures.ThreadPoolExecutor(max_workers=5, thread_name_prefix='rev') as executor:
                jobs = {}
                while rev_files_left:
                    for this_rev_file in rev_files_iter:
                        job = executor.submit(run_reverse, chrom, outdir, this_rev_file, nucl_ex_file)
                        jobs[job] = this_rev_file
                        if len(jobs) > MAX_FILES_PROCESSED:
                            break
                    for job in concurrent.futures.as_completed(jobs):
                        rev_files_left -= 1
                        try:
                            l_tup = job.result()
                        except Exception as ex:
                            logger.error('[rev] Exception {}'.format(ex.__cause__))
                            l_tup = (None,)
                        del jobs[job]
                        d[subdir].append(l_tup[0])
            logger.info('End of reverse concurrent phase for chrom {}'.format(chrom))

    logger.info('Finished processing chrom {}'.format(chrom))
    return d


def submit_process(tup):
    return process(*tup)

