import os
import argparse
import logging
import configparser

logging.basicConfig(format='%(asctime)s - %(levelname)s: %(message)s',
                    filename='sanefalcon.log', filemode='w', level=logging.DEBUG)

logger = logging.getLogger("main")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Launch everything',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', dest='training', action='store_true', help='training phase')
    parser.add_argument('-c', dest='conffile', nargs='?', type=str, default="./sanefalcon.conf",
                        help='path of the configuration file')

    args = parser.parse_args()
    is_training = args.training

    config = configparser.ConfigParser()
    config.read(args.conffile)
