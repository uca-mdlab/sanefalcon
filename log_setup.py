import logging
# from logging.handlers import RotatingFileHandler

formatter = logging.Formatter('%(asctime)s - %(name)-21s - %(levelname)-8s: %(message)s')


def setup_logger(name, log_file, level=logging.DEBUG):

    #handler = RotatingFileHandler(log_file, maxBytes=1e7, backupCount=10)
    handler = logging.FileHandler(log_file)

    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger

