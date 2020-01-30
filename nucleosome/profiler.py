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
        logger.debug(f"{len(input_list)} jobs for computing profiles...")
        logger.info("Launching multiprocessing pool...")
        num_cores = mp.cpu_count()
        with mp.Pool(num_cores) as pool:
            self.profiles = pool.map(gpp.submit_process, input_list)
            logger.info("Done. Result = {}".format(len(self.profiles) == 22))

    def combine(self, training=True):
        if training:
            outfile = self.fm.config['default']['trainnucl']
        else:
            outfile = self.fm.config['default']['testnucl']
        combined, img_file = cp.save_streams_to_file(self.fm, outfile, training)
        return combined, img_file
