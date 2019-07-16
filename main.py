import os
import argparse
import logging
import configparser
from merge import merge_all
from nucleosome_detector import create_nucleosome_files
from getProfileParallel import get_data, submit_process
import multiprocessing as mp
from file_manager import FileManager
from combine_profiles import save_streams_to_file

logging.basicConfig(format='%(asctime)s - %(levelname)s: %(message)s',
                    filename='sanefalcon.log', filemode='w', level=logging.DEBUG)

logger = logging.getLogger("main")


def run_profileParallel(fm, training=True):
    if training:
        nucl_stub = fm.anti_file_template
    else:
        nucl_stub = fm.nucl_file_template
    outfolder = fm.profilefolder
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)
        logger.info('run_profileParallel: Created out folder {}'.format(outfolder))

    data = get_data(fm)  # all the available data

    input_list = [(chrom, dic, outfolder) for chrom, dic in data.items()]
    logger.info("Launching multiprocessing pool...")
    num_cores = mp.cpu_count()
    with mp.Pool(num_cores) as pool:
        finished = pool.map(submit_process, input_list)
        logger.info("Done. Result = {}".format(len(finished) == 22))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Launch everything',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-c', '--conffile', nargs='?', type=str, default="./sanefalcon.conf",
                        help='path of the configuration file')
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(args.conffile)

    fm = FileManager(config)

    logger.info("Starting sanefalcon with configuration file {}".format(args.conffile))

    fm.prepare_train_folder()
    logger.info("prepare_folders ok")

    merge_all(fm)
    logger.info("merge_all ok")

    create_nucleosome_files(fm, training=True)
    logger.info("nucleosome ok")

    run_profileParallel(fm, training=True)  # training: fm.anti_file_template
    logger.info("run profile parallel ok")

    save_streams_to_file(fm, fm.trainnuclfile)