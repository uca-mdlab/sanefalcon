# old: prepSamples.sh + retro.py

import configparser
import os
import subprocess
import sys
import logging
from threading import Thread
import concurrent.futures

logging.basicConfig(format='%(asctime)s - %(levelname)s: %(message)s',
                    filename='sanefalcon.log', filemode='w', level=logging.DEBUG)

logger = logging.getLogger("prepare_samples")


def prepare_fwd(samtools, bamfile, outfile, chrom):
    one = [samtools, "view", bamfile, "chr{}".format(chrom), "-F", "20", "-q", "1"]
    two = [sys.executable, "retro.py", "-fwd"]
    p1 = subprocess.Popen(one, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(two, stdin=p1.stdout, stdout=open(outfile, 'w'))
    output = p2.communicate()[0]


def prepare_rev(samtools, bamfile, revoutfile, chrom):
    one = [samtools, "view", bamfile, "chr{}".format(chrom), "-f", "16", "-F", "4", "-q", "1"]
    two = [sys.executable, "retro.py", "-rev"]
    p1 = subprocess.Popen(one, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(two, stdin=p1.stdout, stdout=open(revoutfile, 'w'))
    output = p2.communicate()[0]


def prepare_samples(datafolder, rspfolder, samtools):
    if not os.path.isdir(rspfolder):
        os.makedirs(rspfolder)
    chromosomes = range(1, 23)
    fwd_jobs = []
    rev_jobs = []

    for root, subdir, files in os.walk(datafolder):
        for fname in files:
            if fname.endswith('.bam'):
                bamfile = os.path.join(root, fname)
                outfile_name = fname.split('.')[0] + '.sort.bam'
                for chrom in chromosomes:
                    outfile = os.path.join(rspfolder, outfile_name + '.{}.start.fwd'.format(chrom))
                    revoutfile = os.path.join(rspfolder, outfile_name + '.{}.start.rev'.format(chrom))
                    if not os.path.isfile(outfile):
                        run = (samtools, bamfile, outfile, chrom)
                        fwd_jobs.append(run)

                    if not os.path.isfile(revoutfile):
                        run = (samtools, bamfile, revoutfile, chrom)
                        rev_jobs.append(run)

                logger.debug("prep_samples for {}".format(fname))

    logger.info('Spawning {} threads'.format(len(fwd_jobs) + len(rev_jobs)))
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        jobs = {}
        for job in fwd_jobs:
            sub = executor.submit(prepare_fwd, *job)
            jobs[sub] = job

        for job in rev_jobs:
            sub = executor.submit(prepare_rev, *job)
            jobs[sub] = job

        for job in concurrent.futures.as_completed(jobs):
            _ = job.result()

    logger.info('End multithreading. Samples prepared.')


if __name__ == "__main__":
    conf_file = 'tests/data/test.conf'
    config = configparser.ConfigParser()
    config.read(conf_file)

    samtools = config['default']['samtools']
    datafolder = config['default']['datafolder']
    trainfolder = config['default']['trainfolder']
    rspfolder = config['default']['rspfolder']

    prepare_samples(datafolder, rspfolder, samtools)
