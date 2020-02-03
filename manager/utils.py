import os
import re
from itertools import islice
import pysam
from log_setup import setup_logger

logger = setup_logger(__name__, 'logs/manager.log')


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
    def compute_num_batches(n):
        divs = [x for x in range(1, n + 1) if n / x == int(n / x)]
        center = (divs[int(len(divs) / 2) - 1], divs[int(len(divs) / 2)])
        low = center[0]
        high = center[1]
        return low, high  # trade offs merge/merge_anti_subs

    @staticmethod
    def count_number_of_reads_per_bam(filename):
        to_avoid = ['chrX', 'chrY', 'chrM']
        bamfile = pysam.AlignmentFile(filename)
        stats = bamfile.get_index_statistics()
        res = {}
        for stat in stats:
            if stat.contig not in to_avoid:
                res[stat.contig] = stat.mapped
        return res, sum(res.values())

    @staticmethod
    def balance_batches(dic, numreads):
        thr = 0.05  # 5%
        num_batches = len(dic)
        tot_reads = sum(numreads.values())
        num_reads_per_batch = {}
        average_reads = tot_reads // num_batches
        average_reads_low = average_reads - average_reads * thr
        average_reads_high = average_reads + average_reads * thr
        logger.debug(f"Average reads (lo-hi) = {average_reads}, ({average_reads_low}, {average_reads_high})")
        high = {}
        low = {}
        ok = {}
        processed = set()

        while len(ok) < num_batches:
            high = {}
            low = {}
            for batch_name, l in dic.items():
                count = [sum([numreads[x] for x in sublist]) for sublist in l]
                assert len(count) == len(l)
                num_reads_per_batch[batch_name] = {'count': sum(count), 'groups': [(x, y) for x, y in zip(l, count)]}

            for batch in num_reads_per_batch:
                if num_reads_per_batch[batch]['count'] > average_reads_high:
                    high[batch] = num_reads_per_batch[batch]['groups']
                elif num_reads_per_batch[batch]['count'] < average_reads_low:
                    low[batch] = num_reads_per_batch[batch]['groups']
                # else:
                #     ok[batch] = num_reads_per_batch[batch]['groups']

            hi_batch_name, hi_group = sorted(high.items(), key=lambda item: item[1], reverse=True)[0]
            lo_batch_name, lo_group = sorted(low.items(), key=lambda item: item[1])[0]
            if (lo_batch_name, hi_batch_name) not in processed:
                logger.debug(f"Swapping from {hi_batch_name} ({num_reads_per_batch[hi_batch_name]['count']}) "
                             f"and {lo_batch_name} ({num_reads_per_batch[lo_batch_name]['count']})")

                processed.add((hi_batch_name, lo_batch_name))
                hi_candidate = sorted(hi_group, key=lambda tup: tup[1], reverse=True)[0]
                lo_candidate = sorted(lo_group, key=lambda tup: tup[1])[0]

                hi_group.remove(hi_candidate)
                hi_group.append(lo_candidate)

                lo_group.remove(lo_candidate)
                lo_group.append(hi_candidate)

                dic[hi_batch_name] = [x[0] for x in hi_group]
                dic[lo_batch_name] = [x[0] for x in lo_group]

            else:
                break

        for k, v in high.items():
            tot = sum(x[1] for x in v)
            logger.debug(f"H, {k}, tot = {tot}, diff = {tot - average_reads}")
        for k, v in low.items():
            tot = sum(x[1] for x in v)
            logger.debug(f"L, {k}, tot = {tot}, diff = {tot - average_reads}")
        for k, v in ok.items():
            tot = sum(x[1] for x in v)
            logger.debug(f"O, {k}, tot = {tot}, diff = {tot - average_reads}")

        return dic

