import logging
from logging.handlers import RotatingFileHandler
import sys

formatter = logging.Formatter('%(asctime)s - %(name)-21s - %(levelname)-8s: %(message)s')


def setup_logger(name, log_file, level=logging.DEBUG):
    if log_file:
        handler = logging.FileHandler(log_file)
    else:
        handler = logging.StreamHandler(sys.stdout)

    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger

