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
    data_files = [f for f in os.listdir(data_dir) if f.endswith('.bam') and f not in to_remove]
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
    # parser.add_argument('-t', dest='training', action='store_true', help='training phase')
    # parser.add_argument("testing_group", help="List of samples in the testing group")
    parser.add_argument('-c', dest='conffile', nargs='?', type=str, default="./sanefalcon.conf",
                        help='path of the configuration file')

    args = parser.parse_args()
    # is_training = args.training

    config = configparser.ConfigParser()
    config.read(args.conffile)
    logger.info('Started sanefalcon with config:')
    for section in config.sections():
        logger.info('=={}=='.format(section))
        for k, v in config[section].items():
            logger.info('{}:\t{}'.format(k, v))

    list_testing_dir = config['folders']['listtestingdir']

    # group_files = [os.path.join(list_testing_dir, f) for f in os.listdir(list_testing_dir)
    #                if re.match('list_testing', f)]

    # logger.info(f'Found {len(group_files)} testing groups')
    data_folder = config['folders']['data']

    group_file = os.path.join(list_testing_dir, 'list_testing1.txt')
    # for group_file in group_files:
    group_name, training_set, testing_set = define_training_and_testing_set(group_file, data_folder)
    logger.debug(f'TESTING {group_name}, {len(testing_set)}, {testing_set[:4]}')
    logger.debug(f'TRAINING {group_name}, {len(training_set)}, {training_set[:4]}')
    model = train(group_name, training_set)
    test(group_name, model, testing_set)
    logger.info(f'{group_name} terminated.')

    # if is_training:
    #     logger.info('Training phase...')
    #     outmodel_file = os.path.join(config['folders']['train'], 'out')
    #     reference_file = config['default']['trainref']
    #     nucleosome_file = monitor.training(config)
    #     logger.info('Nucleosome file computed: {}'.format(nucleosome_file))
    #     modelfileName = predictor.run_model(nucleosome_file, reference_file, outmodel_file)
    #     logger.info('Model trained. {}'.format(outmodel_file))
    # else:
    #     logger.info('Testing phase...')
    #     outmodel_file = os.path.join(config['folders']['train'], 'out.model')
    #     training_nucleosome_file = config['default']['trainnucl']
    #     training_reference_file = config['default']['trainref']
    #     logger.info(f'Model file located at {outmodel_file}')
    #     if not os.path.isfile(training_nucleosome_file):
    #         exit('Run training first. -t flag')
    #     nucleosome_file = monitor.testing(config)
    #     reference_file = config['default']['testref']
    #     test_profile_dir = config['folders']['testprofiles']
    #
    #     logger.info('Testing. training _nucleosome_file = {}; _reference_file = {}'.format(training_nucleosome_file, training_reference_file))
    #     logger.info('Testing. testing _nucleosome_file = {}, _reference_file = {}'.format(nucleosome_file, reference_file))
    #     ffs = compute_ff(outmodel_file, test_profile_dir)
    #     for k, v in ffs.items():
    #         print(k, v)
    #     # predictor.run_model(training_nucleosome_file, training_reference_file, outmodel_file, nucleosome_file, reference_file)

    logger.info('Done.')



