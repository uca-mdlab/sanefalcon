import os
from log_setup import setup_logger
import pickle
import heapq

logger = setup_logger(__name__, 'logs/merge_utils.log')


class Utils:

    @staticmethod
    def write_to_pickle(data, outfile):
        pickle.dump(data, open(outfile, 'wb'))

    @staticmethod
    def sort_and_write(data, outfile):
        # data = list(set(data))  # remove duplicates and sort
        data.sort()
        Utils.write_to_pickle(data, outfile)
        # pickle.dump(data, open(outfile, 'wb'))
        # to_be_written = map(lambda x: str(x) + "\n", data)
        # with open(outfile, 'w') as out:
        #     out.writelines(to_be_written)

    @staticmethod
    def read_all_files(file_list, is_merge=False):
        data = []
        if not is_merge:
            logger.debug(f"Reading {len(file_list)} plain text files: {file_list[:4]}")
            for f in file_list:
                data.extend([int(line.strip()) for line in open(f, 'r')])  # rsp are in plain text
        else:
            logger.debug(f"Reading {len(file_list)} pickle binary files: {file_list[:4]}")
            data = Utils.read_pickle_data(file_list)
        return data

    @staticmethod
    def read_pickle_data(file_list):
        data = []
        for f in file_list:
            tmp = pickle.load(open(f, 'rb'))
            data.append(tmp)

        return [item for item in heapq.merge(*data)]

    @staticmethod
    def pick_manip_pairs(manips, list_):
        res = {}
        for subdir, manip_list in manips.items():
            for manip in manip_list:
                files = []
                manip_name = os.path.basename(manip)
                for path in list_:
                    base = os.path.basename(path)
                    if os.path.commonprefix([base, manip_name]) == manip_name:
                        files.append(path)
                res[manip] = files
        return res

    @staticmethod
    def prepare_jobs_arguments(files_to_merge):
        runs = []
        for subdir, dic in files_to_merge.items():
            logger.debug(f"prepare_job_arguments for {subdir}")
            logger.debug(f"{dic.items()}")
            for chrom in range(1, 23):
                fwd = dic['fwd'][chrom]
                rev = dic['rev'][chrom]
                files = fwd + rev
                runs.append((files, subdir, chrom))
        return runs
