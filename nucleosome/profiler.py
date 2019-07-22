import multiprocessing as mp
import logging

import nucleosome.get_profile_parallel as gpp
import nucleosome.combine_profiles as cp


logger = logging.getLogger(__name__)


class Profiler:
    def __init__(self, fm, tracker):
        self.fm = fm
        self.tracker = tracker
        self.profiles = []  # list of defaultdict

    def compute_profiles(self, mapping):
        data = self.tracker.get_data(mapping)
        input_list = [(chrom, dic, self.fm.profilefolder) for chrom, dic in data.items()]
        logger.info("Launching multiprocessing pool...")
        num_cores = mp.cpu_count()
        with mp.Pool(num_cores) as pool:
            self.profiles = pool.map(gpp.submit_process, input_list)
            logger.info("Done. Result = {}".format(len(self.profiles) == 22))

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

    def combine(self):
        outfile = self.fm.config['default']['trainnucl']
        cp.save_streams_to_file(self.fm, outfile)
