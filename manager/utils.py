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
        it = iter(lis)
        return iter(lambda: tuple(islice(it, size)), ())

    @staticmethod
    def compute_batch_size(n):
        divs = [x for x in range(1, n + 1) if n / x == int(n / x)]
        center = (divs[int(len(divs) / 2) - 1], divs[int(len(divs) / 2)])
        return center[1]

