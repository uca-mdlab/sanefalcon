import logging

formatter = logging.Formatter('%(asctime)s - %(name)-21s - %(levelname)-8s: %(message)s')


def setup_logger(name, log_file, level=logging.DEBUG):

    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger

