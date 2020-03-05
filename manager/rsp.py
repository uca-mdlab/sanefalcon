import os
import sys
import subprocess
import concurrent.futures

from log_setup import setup_logger
import pysam
from retro import in_mem_retro

logger = setup_logger(__name__, 'logs/manager.log')


class RspBuilder:
    def __init__(self, config):
        pass
        # self.samtools = config['default']['samtools']

    def prepare_fwd(self, bamfile, outfile, chrom):
        reads = []
        infile = pysam.AlignmentFile(bamfile, 'rb')
        for read in infile.fetch(f'chr{chrom}'):
            if not read.is_reverse and read.mapping_quality >= 1:
                start = read.reference_start + 1  # pysam is 0-indexed. samtools 1-indexed
                length = read.reference_length
                reads.append((start, length))

        in_mem_retro(reads, outfile, flag=True)
        return outfile

    def prepare_rev(self, bamfile, outfile, chrom):
        reads = []
        infile = pysam.AlignmentFile(bamfile, 'rb')
        for read in infile.fetch(f'chr{chrom}'):
            if read.is_reverse and read.mapping_quality >= 1:
                start = read.reference_start + 1
                length = read.query_length
                reads.append((start, length))
        in_mem_retro(reads, outfile, False)
        return outfile

    def prepare_fwd_rev(self, bamfile, outfile, chrom):
        forward = []
        reverse = []
        infile = pysam.AlignmentFile(bamfile, 'rb')
        for read in infile.fetch(f'chr{chrom}'):
            if read.mapping_quality >= 1:
                start = read.reference_start + 1
                length = read.query_length
                if read.is_reverse:
                    reverse.append((start, length))
                else:
                    forward.append((start, length))
        fwdfile = in_mem_retro(forward, outfile + '.fwd', True)
        revfile = in_mem_retro(reverse, outfile + '.rev', False)
        return fwdfile, revfile

    # def prepare_fwd(self, bamfile, outfile, chrom):
    #     one = [self.samtools, "view", bamfile, "chr{}".format(chrom), "-F", "20", "-q", "1"]
    #     two = [sys.executable, "retro.py", "-fwd"]
    #     p1 = subprocess.Popen(one, stdout=subprocess.PIPE)
    #     with open(outfile, 'w') as out_:
    #         p2 = subprocess.Popen(two, stdin=p1.stdout, stdout=out_)
    #         output = p2.communicate()[0]
    #     return outfile

    # def prepare_rev(self, bamfile, revoutfile, chrom):
    #     one = [self.samtools, "view", bamfile, "chr{}".format(chrom), "-f", "16", "-F", "4", "-q", "1"]
    #     two = [sys.executable, "retro.py", "-rev"]
    #     p1 = subprocess.Popen(one, stdout=subprocess.PIPE)
    #     with open(revoutfile, 'w') as out_:
    #         p2 = subprocess.Popen(two, stdin=p1.stdout, stdout=out_)
    #         output = p2.communicate()[0]
    #     return revoutfile

    def prepare_samples(self, datafolder, rspfolder):
        chromosomes = range(1, 23)
        todo = []
        outfiles = []
        for root, subdir, files in os.walk(datafolder):
            for fname in files:
                if fname.endswith('.bam'):
                    bamfile = os.path.join(root, fname)
                    for chrom in chromosomes:
                        template = os.path.join(rspfolder, fname + f'.{chrom}.start')
                        outfile = template + '.fwd'
                        revoutfile = template + '.rev'
                        if not os.path.isfile(outfile) or not os.path.isfile(revoutfile):
                            run = (bamfile, template, chrom)
                            todo.append(run)

                    logger.debug("prep_samples for {}".format(fname))

        tot = len(todo)
        if tot > 0:
            logger.info('Spawning {} threads'.format(tot))
            with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
                jobs = {}
                for job in todo:
                    sub = executor.submit(self.prepare_fwd_rev, *job)
                    jobs[sub] = job

                for job in concurrent.futures.as_completed(jobs):
                    outfile = job.result()  # this is a pair (fwd, rev)
                    outfiles.extend(outfile)

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
    #rspbuilder.prepare_fwd(infile, '/data/tempff/1.new.fwd', chrom)
    #rspbuilder.prepare_rev(infile, '/data/tempff/1.new.rev', chrom)
    rspbuilder.prepare_fwd_rev(infile, '/data/tempff/1.test', chrom)
