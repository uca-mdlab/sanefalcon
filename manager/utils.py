import os
import re
from itertools import islice

class Utils:

    @staticmethod
    def get_data(folder, mask=''):
        pattern = re.compile(mask)
        return [os.path.join(folder, x) for x in os.listdir(folder) if re.search(pattern, x)]

    @staticmethod
    def readfile(file):
        """
        Read a file and returns a list of stripped lines.
        Tries 'int' conversion for read-start positions
        :param file: the file to read
        :return: a list of either 'str' or 'int'
        """
        with open(file) as f:
            lines = [x.strip() for x in f.readlines()]
        f.close()

        try:
            lines = list(map(int, lines))
        except ValueError:
            pass

        return lines

    @staticmethod
    def prepare_batches(lis, size):
        print(lis)
        it = iter(lis)
        return iter(lambda: tuple(islice(it, size)), ())

