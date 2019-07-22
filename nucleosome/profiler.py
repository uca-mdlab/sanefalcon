import multiprocessing as mp
import logging

import nucleosome.get_profile_parallel as gpp


logger = logging.getLogger(__name__)


class Profiler:
    def __init__(self, fm, tracker):
        self.fm = fm
        self.tracker = tracker

    def compute_profiles(self, mapping):
        data = self.tracker.get_data(mapping)
        input_list = [(chrom, dic, self.fm.profilefolder) for chrom, dic in data.items()]
        logger.info("Launching multiprocessing pool...")
        num_cores = mp.cpu_count()
        with mp.Pool(num_cores) as pool:
            finished = pool.map(gpp.submit_process, input_list)
            logger.info("Done. Result = {}".format(len(finished) == 22))

            # finished = [
            #  { '/var/tmp/testsuite/a': [
            #      ('/var/tmp/testsuite/profiles/p18_084.bam.1.start.fwd.1.fwd',
            #       '/var/tmp/testsuite/profiles/p18_084.bam.1.start.fwd.1.ifwd'),
            #      ('/var/tmp/testsuite/profiles/p18_084.bam.1.start.rev.1.irev',
            #       '/var/tmp/testsuite/profiles/p18_084.bam.1.start.rev.1.rev')
            #     ],
            #    '/var/tmp/testsuite/b': [
            #      ('/var/tmp/testsuite/profiles/p18_085.bam.1.start.fwd.1.fwd',
            #       '/var/tmp/testsuite/profiles/p18_085.bam.1.start.fwd.1.ifwd'),
            #      ('/var/tmp/testsuite/profiles/p18_085.bam.1.start.rev.1.irev',
            #       '/var/tmp/testsuite/profiles/p18_085.bam.1.start.rev.1.rev')
            #     ]
            #   }
            # ]
