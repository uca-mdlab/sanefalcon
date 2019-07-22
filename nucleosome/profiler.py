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
