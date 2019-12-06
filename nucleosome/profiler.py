import multiprocessing as mp


import nucleosome.get_profile_parallel as gpp
import nucleosome.combine_profiles as cp


from log_setup import setup_logger
logger = setup_logger(__name__, 'logs/nucleosome.log')


class Profiler:
    def __init__(self, fm, tracker):
        self.fm = fm
        self.tracker = tracker
        self.profiles = []  # list of defaultdict

    def compute_profiles(self, mapping, training=True):
        data = self.tracker.get_data(mapping)
        if training:
            profilefolder = self.fm.profilefolder
        else:
            profilefolder = self.fm.testprofilefolder
        input_list = [(chrom, dic, profilefolder) for chrom, dic in data.items()]
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

    def combine(self, training=True):
        if training:
            outfile = self.fm.config['default']['trainnucl']
        else:
            outfile = self.fm.config['default']['testnucl']
        combined, img_file = cp.save_streams_to_file(self.fm, outfile)
        return combined, img_file
