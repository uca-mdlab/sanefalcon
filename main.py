import argparse
import configparser
import monitor
import predictor
import os
from log_setup import setup_logger
from test_model import compute_ff
import re

logger = setup_logger(__name__, "logs/sanefalcon.log")


def define_training_and_testing_set(test_file, data_dir):
    """
    Given a list of file for testing, find all the other samples in the cohorte, to be used for training
    :param test_file: file of testing samples
    :return: list of files
    """
    to_remove = ['p2_096_dedup.bam', 'p42_089_dedup.bam']
    group_name = os.path.basename(test_file).split('.')[0]
    data_files = [f for f in os.listdir(data_dir) if f.endswith('.bam') and f not in to_remove][:160] # FIXME remove [:150]
    logger.debug(f'Cohorte consists of {len(data_files)} samples')
    with open(test_file, 'r') as groupfile:
        testing_set = [x.strip().split(';')[0] + '.bam' for x in groupfile.readlines()]

    logger.debug(f'Testing set for {group_name}: {len(testing_set)} samples')
    training_set = [os.path.join(data_dir, f) for f in data_files if f not in testing_set]
    logger.debug(f'Training set for {group_name}: {len(training_set)} samples')
    return group_name, training_set, testing_set


def train(group_name, training_set):
    logger.info(f'Training phase {group_name}')

    outmodel_file = os.path.join(config['folders']['train'], group_name)
    reference_file = config['default']['trainref']
    logger.debug(f'Final model will be saved in {outmodel_file}')
    logger.debug(f'Reference file:  {reference_file}')
    nucleosome_file = monitor.training(config, training_set)

    logger.info('Nucleosome file computed: {}'.format(nucleosome_file))
    modelfileName = predictor.run_model(nucleosome_file, reference_file, outmodel_file)
    logger.info(f'Model trained for {group_name}: {modelfileName}')
    return modelfileName


def test(group_name, outmodel_file, testing_set):
    logger.info(f'Testing phase {group_name}')
    nucleosome_file = monitor.testing(config, testing_set)
    test_profile_dir = config['folders']['testprofiles']
    ffs = compute_ff(outmodel_file, test_profile_dir)
    for k, v in ffs.items():
        print(k, v)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Launch everything',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("testing_group_file", help="List of samples in the testing group")
    parser.add_argument('-c', dest='conffile', nargs='?', type=str, default="./sanefalcon.conf",
                        help='path of the configuration file')

    args = parser.parse_args()
    group_file = args.testing_group_file
    if not os.path.isfile(group_file):
        exit(f'{group_file} not found.')

    config = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
    config.read(args.conffile)
    logger.info('Started sanefalcon with config:')
    for section in config.sections():
        logger.info('=={}=='.format(section))
        for k, v in config[section].items():
            logger.info('{}:\t{}'.format(k, v))

    data_folder = config['folders']['data']

    group_name, training_set, testing_set = define_training_and_testing_set(group_file, data_folder)
    logger.debug(f'TESTING {group_name}, {len(testing_set)}, {testing_set[:4]}')
    logger.debug(f'TRAINING {group_name}, {len(training_set)}, {training_set[:4]}')
    model = train(group_name, training_set)
    test(group_name, model, testing_set)
    logger.info(f'{group_name} terminated.')

    logger.info('Done.')



