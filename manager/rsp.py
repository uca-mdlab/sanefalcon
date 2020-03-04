import os
import sys
import subprocess
import concurrent.futures

from log_setup import setup_logger
import pysam
import pickle
from retro import in_mem_retro

logger = setup_logger(__name__, 'logs/manager.log')


class RspBuilder:
    def __init__(self, config):
        self.samtools = config['default']['samtools']

    def prepare_fwd_new(self, bamfile, outfile, chrom):
        reads = []
        infile = pysam.AlignmentFile(bamfile, 'rb')
        for read in infile.fetch(f'chr{chrom}'):
            if not read.is_reverse and read.mapping_quality >= 1:
                start = read.reference_start + 1  # pysam is 0-indexed. samtools 1-indexed
                length = read.reference_length
                reads.append((start, length))
        return reads

    def prepare_fwd(self, bamfile, outfile, chrom):
        reads = self.prepare_fwd_new(bamfile, outfile, chrom)
        in_mem_retro(reads, outfile, flag=True)

    def prepare_rev_new(self, bamfile, outfile, chrom):
        reads = []
        infile = pysam.AlignmentFile(bamfile, 'rb')
        for read in infile.fetch(f'chr{chrom}'):
            if read.is_reverse and read.mapping_quality >= 1:
                start = read.reference_start
                length = read.reference_length + 1
                reads.append((start, length))
        return reads

    def prepare_rev_tmp(self, bamfile, outfile, chrom):
        reads = self.prepare_rev_new(bamfile, outfile, chrom)
        in_mem_retro(reads, outfile, False)

    # def prepare_fwd(self, bamfile, outfile, chrom):
    #     one = [self.samtools, "view", bamfile, "chr{}".format(chrom), "-F", "20", "-q", "1"]
    #     two = [sys.executable, "retro.py", "-fwd"]
    #     p1 = subprocess.Popen(one, stdout=subprocess.PIPE)
    #     with open(outfile, 'w') as out_:
    #         p2 = subprocess.Popen(two, stdin=p1.stdout, stdout=out_)
    #         output = p2.communicate()[0]
    #     return outfile

    def prepare_rev(self, bamfile, revoutfile, chrom):
        one = [self.samtools, "view", bamfile, "chr{}".format(chrom), "-f", "16", "-F", "4", "-q", "1"]
        two = [sys.executable, "retro.py", "-rev"]
        p1 = subprocess.Popen(one, stdout=subprocess.PIPE)
        with open(revoutfile, 'w') as out_:
            p2 = subprocess.Popen(two, stdin=p1.stdout, stdout=out_)
            output = p2.communicate()[0]
        return revoutfile

    def prepare_samples(self, datafolder, rspfolder):
        chromosomes = range(1, 23)
        fwd_jobs = []
        rev_jobs = []
        outfiles = []
        for root, subdir, files in os.walk(datafolder):
            for fname in files:
                if fname.endswith('.bam'):
                    bamfile = os.path.join(root, fname)
                    for chrom in chromosomes:
                        outfile = os.path.join(rspfolder, fname + '.{}.start.fwd'.format(chrom))
                        revoutfile = os.path.join(rspfolder, fname + '.{}.start.rev'.format(chrom))
                        if not os.path.isfile(outfile):
                            run = (bamfile, outfile, chrom)
                            fwd_jobs.append(run)
                        if not os.path.isfile(revoutfile):
                            run = (bamfile, revoutfile, chrom)
                            rev_jobs.append(run)

                    logger.debug("prep_samples for {}".format(fname))

        tot = len(fwd_jobs) + len(rev_jobs)
        if tot > 0:
            logger.info('Spawning {} threads'.format(tot))
            with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
                jobs = {}
                for job in fwd_jobs:
                    sub = executor.submit(self.prepare_fwd, *job)
                    jobs[sub] = job

                for job in rev_jobs:
                    sub = executor.submit(self.prepare_rev, *job)
                    jobs[sub] = job

                for job in concurrent.futures.as_completed(jobs):
                    outfile = job.result()
                    outfiles.append(outfile)

            logger.info('End multithreading. Samples prepared.')
        else:
            logger.info('Samples already prepared.')

        if len(outfiles) == 0:
            return [os.path.join(rspfolder, fname) for fname in os.listdir(rspfolder)]
        else:
            return outfiles


if __name__ == '__main__':
    config = {'default': {'samtools': '/user/mmilanes/home/sw/samtools-1.9/samtools'}}
    infile = '/data/tempff/bam/p1_065_dedup.bam'
    chrom = 1
    rspbuilder = RspBuilder(config)
    #rspbuilder.prepare_fwd_tmp(infile, '/data/tempff/1.tmp.fwd', chrom)
    rspbuilder.prepare_rev_tmp(infile, '/data/tempff/1.tmp.rev', chrom)
    #fwdfile = rspbuilder.prepare_fwd(infile, '/data/tempff/1.fwd', chrom)
    #revfile = rspbuilder.prepare_rev(infile, '/data/tempff/1.rev', chrom)
