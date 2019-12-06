import argparse
import configparser
import monitor
import predictor
import os
from log_setup import setup_logger

logger = setup_logger(__name__, "logs/sanefalcon.log")


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
    logger.info('Started sanefalcon with config:')
    for section in config.sections():
        logger.info('=={}=='.format(section))
        for k, v in config[section].items():
            logger.info('{}:\t{}'.format(k, v))
    exit()
    if is_training:
        logger.info('Training phase...')
        outmodel_file = os.path.join(config['folders']['train'], 'out')
        reference_file = config['default']['trainref']
        nucleosome_file = monitor.training(config)
        logger.info('Nucleosome file computed: {}'.format(nucleosome_file))
        predictor.run_model(nucleosome_file, reference_file, outmodel_file)
        logger.info('Model trained. {}'.format(outmodel_file))
    else:
        outmodel_file = os.path.join(config['folders']['test'], 'out')
        training_nucleosome_file = config['default']['trainnucl']
        training_reference_file = config['default']['trainref']

        if not os.path.isfile(training_nucleosome_file):
            exit('Run training first. -t flag')
        nucleosome_file = monitor.testing(config)
        reference_file = config['default']['testref']
        logger.info('Testing. training _nucleosome_file = {}; _reference_file = {}'.format(training_nucleosome_file, training_reference_file))
        logger.info('Testing. testing _nucleosome_file = {}, _reference_file = {}'.format(nucleosome_file, reference_file))
        # predictor.run_model(training_nucleosome_file, training_reference_file, outmodel_file, nucleosome_file, reference_file)

    logger.info('Done.')



