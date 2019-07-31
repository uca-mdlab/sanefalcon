import argparse
import logging
import configparser
import monitor
import predictor
import os


logging.basicConfig(format='%(asctime)s - %(name)-21s - %(levelname)-8s: %(message)s',
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

    if is_training:
        outmodel_file = os.path.join(config['folders']['train'], 'out')
        reference_file = config['default']['trainref']
        nucleosome_file = monitor.training(config)
        predictor.run_model(nucleosome_file, reference_file, outmodel_file)
    else:
        outmodel_file = os.path.join(config['folders']['test'], 'out')
        training_nucleosome_file = config['default']['trainnucl']
        training_reference_file = config['default']['trainref']
        if not os.path.isfile(training_nucleosome_file):
            exit('Run training first. -t flag')
        nucleosome_file = monitor.testing(config)
        reference_file = config['default']['testref']
        predictor.run_model(training_nucleosome_file, training_reference_file, outmodel_file, nucleosome_file, reference_file)




