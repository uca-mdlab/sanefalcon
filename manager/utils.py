import os
import re
from itertools import islice
import pysam


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
        print(average_reads, average_reads_low, average_reads_high)
        high = {}
        low = {}
        ok = {}
        for batch_name, l in dic.items():
            count = [sum([numreads[x] for x in sublist]) for sublist in l]
            assert len(count) == len(l)
            num_reads_per_batch[batch_name] = [(x, y) for x, y in zip(l, count)]
            if sum(count) > average_reads_high:
                high[batch_name] = num_reads_per_batch[batch_name]
            elif sum(count) < average_reads_low:
                low[batch_name] = num_reads_per_batch[batch_name]
            else:
                ok[batch_name] = num_reads_per_batch[batch_name]

        return high, low, ok, average_reads

if __name__ == "__main__":
    import pickle
    import sys
    dic = pickle.load(open(sys.argv[1], 'rb'))
    numreads = pickle.load(open(sys.argv[2], 'rb'))
    high, low, ok, average_reads = Utils.balance_batches(dic, numreads)
    print('Average number of reads:', average_reads)
    high_totals = {}
    low_totals = {}
    ok_totals = {}
    for k, v in high.items():
        print('H', k, sum(x[1] for x in v))
        high_totals[k] = sum(x[1] for x in v)
    for k, v in low.items():
        print('L', k, sum(x[1] for x in v))
        low_totals[k] = sum(x[1] for x in v)
    for k, v in ok.items():
        print('O', k, sum(x[1] for x in v))
        ok_totals[k] = sum(x[1] for x in v)

    high_ordered = sorted(high_totals.items(), key=lambda item: item[1], reverse=True)
    low_ordered = sorted(low_totals.items(), key=lambda item: item[1])

    hi_key = high_ordered[0][0]
    lo_key = low_ordered[0][0]
    print('swap groups between (H)', hi_key, 'and (L)', lo_key)

    old_hi = sorted(high[hi_key], key=lambda tup: tup[1], reverse=True)[0]
    old_lo = sorted(low[lo_key], key=lambda tup: tup[1])[0]

    print('Forecast: ', hi_key, ' = ', high_totals[hi_key], high_totals[hi_key] - old_hi[1] + old_lo[1])
    print('Forecast: ', lo_key, ' = ', low_totals[lo_key], low_totals[lo_key] - old_lo[1] + old_hi[1])
    # for k, v in res.items():
    #     print(k, v)

