import os
import sys
import logging
import re
import concurrent.futures
import threading
import argparse

logger = logging.getLogger('getProfileParallel')
logging.basicConfig(format='%(asctime)s - %(levelname)s:%(message)s',
                    level=logging.DEBUG,
                    filename='profileparallel.log',
                    filemode='w')

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
    # int(nuclLine[0])


def cast_line_to_numbers(list_of_strings):
    def cast_to_number(string_literal):
        try:
            return int(string_literal)
        except ValueError:
            return float(string_literal)

    return list(map(cast_to_number, list_of_strings))


def load_data(nucl_ex_file, fwd_rev_file):
    with open(nucl_ex_file, 'r') as infile:
        lines = list(map(cast_line_to_numbers, [x.split() for x in infile.readlines()]))

        # peaks = [loadNucl(l) for l in lines if nuclFilt(l[2:])]
        peaks=[]
        for l in lines:
            if nuclFilt(l[2:]):
                l[0]=int(l[0])
                peaks+=[loadNucl(l)]
        # print(peaks,"TESTPEAKS")

    with open(fwd_rev_file, 'r') as infile:
        reads = [int(line) for line in infile]
    # reads = [int(line) for line in open(fwd_rev_file)]
    logger.debug("{}: {} peaks, {} filtOut".format(os.path.basename(fwd_rev_file), len(peaks), filtOut))
    return lines, peaks, reads


def process_forward(peaks, reads, outfile):
    thread_name = threading.current_thread().name
    before = len(peaks)  # just to be sure
    peaks.append([-1, 0])
    logger.debug('process_forward {} (fwd append) id: {} (len += {}) [0] = {}; [-1] = {}'.format(thread_name, id(peaks), len(peaks) - before, peaks[0], peaks[-1]))

    maxDist = 147
    sumPeak = [0. for _ in range(maxDist)]
    read = reads[0]  # +shift
    j = 0
    nuclHit = []
    for i, peakPair in enumerate(peaks):
        if i == 0 or i == len(peaks) - 1:
            logger.debug('process_forward {} (fwd enumerate id: {}) i: {} peakPair {}'.format(thread_name, id(peaks), i, peakPair))
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


def get_data(train_folder, outfolder, nucl_stub):
    """
    Avoid calling multiple times the getProfileSubmit.sh script, as it loads all the data in memory.

    :param train_folder: sanefalcontrain/a
    :param outfolder: /tmp/...
    :return: a dictionary with all the data packed and organized for processing
    """
    # nucl_files=[]
    # nucl_files = [os.path.join(train_folder, f) for f in os.listdir(train_folder) if os.path.isfile(os.path.join(train_folder, f)) and f.startswith(nucl_stub)]
    # nucl_files = [os.path.join(train_folder, f) for root,sub,f in os.walk(train_folder) if os.path.isfile(os.path.join(train_folder, f)) and f.startswith(nucl_stub)]
    for root, subdir, files in os.walk(train_folder):
        for f in files:
            if os.path.isfile(os.path.join(root, f)) and f.startswith(nucl_stub):
                nucl_files.append(os.path.join(train_folder, f))
    print(nucl_files,"NUCL_FILE")
    print(train_folder,"trainfolder")
    print(outfolder,"outfolder")
    print(nucl_stub,"nuclstub")

    logger.debug('Found {} nucl_files'.format(len(nucl_files)))
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)
    fwd_files = []
    rev_files = []
    for root, subdir, files in os.walk(train_folder):
        for f in files:
            if f.endswith('fwd'):
                fwd_files.append(os.path.join(root, f))
            elif f.endswith('rev'):
                rev_files.append(os.path.join(root, f))
            else:
                continue
    logger.debug('Found {} .start.fwd and {} .start.rev files'.format(len(fwd_files), len(rev_files)))
    chromosomes = range(1, 23) ## 1, 23
    d = dict.fromkeys(chromosomes)
    for c in chromosomes:
        nucl_file = [f for f in nucl_files if int(f.split('.')[-1]) == c][0]  # get the nucl file for the chromosome
        regexp = re.compile(".bam.{}.start".format(c))
        c_fwd_files = [f for f in fwd_files if re.search(regexp, f)]          # add fwd files
        c_rev_files = [f for f in rev_files if re.search(regexp, f)]          # add rev files
        logger.debug("Consistency c_fwd_files, c_rev_files = {} - {}".format(all([f.endswith('.fwd') for
                                                                                  f in c_fwd_files]),
                                                                             all([f.endswith('.rev')
                                                                                  for f in c_rev_files])))
        d[c] = {'nucl_file': nucl_file, 'fwd': c_fwd_files, 'rev': c_rev_files}  # pack together in a dictionary
        logger.debug('Chrom {}: {} fwd files, {} rev files'.format(c, len(c_fwd_files), len(c_rev_files)))
    return d


def run_forward(chrom, outdir, fwd_file, nucl_ex_file):
    output_stubname = os.path.join(outdir, os.path.basename(fwd_file) + '.{}'.format(chrom))
    logger.debug("run_forward chrom {}. outdir={}, fwd_file={}, nucl_ex_file={}, "
                 "output_stubname={}".format(chrom, outdir, fwd_file,nucl_ex_file, output_stubname))
    lines, peaks, reads = load_data(nucl_ex_file, fwd_file)
    logger.debug('run_forward chrom {}-{}, peaks id {}'.format(nucl_ex_file.split('.')[-1], os.path.basename(fwd_file), id(peaks)))
    sumPeakFwd = []
    sumPeakRev = []
    fwd_out = output_stubname + ".fwd"
    try:
        sumPeakFwd = process_forward(peaks, reads, fwd_out)  # -> .fwd
        logger.debug(
            "run_forward - chrom {} id {}: process_forward completed: sumPeaks [{}..{}]".format(chrom, id(peaks),
                                                                                                sumPeakFwd[0],
                                                                                                sumPeakFwd[-1]))
    except IndexError:
        logger.error("run_forward IndexError: .fwd: peaks {}, reads {}".format(len(peaks), len(reads)))
    ifwd_out = output_stubname + ".ifwd"
    try:
        sumPeakRev = process_reverse(peaks, reads, ifwd_out)  # -> .ifwd
        logger.debug("run_forward - chrom {} id {}: process_reverse completed: sumPeaks [{}..{}]".format(chrom, id(peaks), sumPeakRev[0], sumPeakRev[-1]))
    except IndexError:
        logger.error("IndexError: .ifwd: peaks {}, reads {}".format(len(peaks), len(reads)))

    return len(sumPeakFwd), len(sumPeakRev)


def run_reverse(chrom, outdir, rev_file, nucl_ex_file):
    output_stubname = os.path.join(outdir, os.path.basename(rev_file) + '.{}'.format(chrom))
    lines, peaks, reads = load_data(nucl_ex_file, rev_file)
    logger.debug('run_reverse chrom {}-{}. peaks id {}'.format(nucl_ex_file.split('.')[-1], os.path.basename(rev_file), id(peaks)))
    sumPeakFwd = []
    sumPeakRev = []
    irev_out = output_stubname + ".irev"
    try:
        sumPeakFwd = process_forward(peaks, reads, irev_out)  # -> .irev
        logger.debug("run_reverse - chrom {} id {}: process_forward completed: sumPeaks [{}..{}]".format(chrom, id(peaks), sumPeakFwd[0], sumPeakFwd[-1]))
    except IndexError:
        logger.error("IndexError: .irev: peaks {}, reads {}".format(len(peaks), len(reads)))
    rev_out = output_stubname + ".rev"
    try:
        sumPeakRev = process_reverse(peaks, reads, rev_out)  # -> .rev
        logger.debug("run_reverse - chrom {} id {}: process_reverse completed: sumPeaks [{}..{}]".format(chrom, id(peaks), sumPeakRev[0], sumPeakRev[-1]))
    except IndexError:
        logger.error("IndexError: .rev: peaks {}, reads {}".format(len(peaks), len(reads)))
    return len(sumPeakFwd), len(sumPeakRev)


def process(chrom, dic, outdir):
    MAX_FILES_PROCESSED = 10  # limit parallelism not to overload memory

    nucl_ex_file = dic['nucl_file']
    fwd_files = dic['fwd']
    fwd_files_left = len(fwd_files)
    fwd_files_iter = iter(fwd_files)

    rev_files = dic['rev']
    rev_files_left = len(rev_files)
    rev_files_iter = iter(rev_files)

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
                this_fwd_file = jobs[job]
                del jobs[job]
                logger.info('[fwd] Got future for'.format(this_fwd_file))

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
                this_rev_file = jobs[job]
                del jobs[job]
                logger.info('[rev] Got future for'.format(this_rev_file))

    logger.info('End of reverse concurrent phase for chrom {}'.format(chrom))
    logger.info('Finished processing chrom {}'.format(chrom))
    return chrom


def submit_process(tup):
    return process(*tup)


if __name__ == "__main__":
    import multiprocessing as mp

    parser = argparse.ArgumentParser(description='Get Nucleosome Profile for Sanefalcon')
    parser.add_argument('train_folder', type=str, help='Folder containing data for the Nucleosome profile')
    parser.add_argument('outfolder', type=str, help='Output data folder')
    parser.add_argument('-n', '--nucltrack', dest='nucl_stub', type=str, default="nucl_ex4",
                        help='Give a nucltrack name without the extension')
    args = parser.parse_args()

    # train_folder = sys.argv[1]
    train_folder = args.train_folder
    # outfolder = sys.argv[2]
    outfolder = args.outfolder
    logger.info('Starting on {}'.format(train_folder))

    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)
        logger.info('Created out folder {}'.format(outfolder))

    nucl_stub = args.nucl_stub
    data = get_data(train_folder, outfolder, nucl_stub)  # all the available data

    # for chrom, dic in data.items():
    #     process(chrom, dic, outfolder)

    input_list = [(chrom, dic, outfolder) for chrom, dic in data.items()]

    num_cores = mp.cpu_count()
    #num_cores = 2
    with mp.Pool(num_cores) as pool:
        finished = pool.map(submit_process, input_list)
        print(len(finished) == 22)
