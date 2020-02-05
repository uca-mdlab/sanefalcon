import os
import re
from itertools import islice
import pysam
from log_setup import setup_logger
import matplotlib.pyplot as plt
import string

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

        sizes = []
        for k, v in dic.items():
            s = 0
            for list_ in v:
                s += sum(numreads[x] for x in list_)
            sizes.append(s)
            logger.debug(f'Batch {k}. Num_reads = {s}')

        return dic  # , sizes

    @staticmethod
    def knapsack_ratio_greedy(number, capacity, weight_cost):
        """Greedy ratio method
        :param number: number of existing items
        :param capacity: the capacity of knapsack
        :param weight_cost: list of tuples like: [(weight, cost), (weight, cost), ...]
        :return: tuple like: (best cost, best combination list(contains 1 and 0))
        """
        ratios = [(index, item[1] / float(item[0])) for index, item in enumerate(weight_cost)]
        ratios = sorted(ratios, key=lambda x: x[1], reverse=True)
        best_combination = [0] * number
        best_cost = 0
        weight = 0
        for index, ratio in ratios:
            if weight_cost[index][0] + weight <= capacity:
                weight += weight_cost[index][0]
                best_cost += weight_cost[index][1]
                best_combination[index] = 1

        return best_cost, best_combination

    @staticmethod
    def balance_w_knapsack(dic, numreads):
        num_reads_per_batch = {}
        for batch_name, l in dic.items():
            count = [sum([numreads[x] for x in sublist]) for sublist in l]
            assert len(count) == len(l)
            num_reads_per_batch[batch_name] = {'count': sum(count), 'groups': [(x, y) for x, y in zip(l, count)]}
        total_numer_of_reads = sum(numreads.values())
        thr = 0.05
        average_reads = total_numer_of_reads // len(num_reads_per_batch)
        knapsack_capacity = average_reads + average_reads * thr

        # print(average_reads, knapsack_capacity)
        groups = [x['groups'] for x in num_reads_per_batch.values()]

        items = []
        for group in groups:
            for tup in group:
                items.append(tup)

        tmp_batches = []
        copy_of_items = items.copy()
        while len(copy_of_items) > 0:
            number = len(items)
            weight_cost = [(t[1], 1) for t in copy_of_items]
            best_cost, best_combination = Utils.knapsack_ratio_greedy(number, knapsack_capacity, weight_cost)
            indexes = [index for index, v in enumerate(best_combination) if v == 1]
            tmp = []
            for index in indexes:
                tmp.append(copy_of_items[index][0])

            indexes.sort(reverse=True)
            tmp_batches.append(tmp)
            for i in indexes:
                del copy_of_items[i]

        # tmp_batches.append([x[0] for x in copy_of_items])

        sizes = []
        for newb in tmp_batches:
            size = 0
            for l in newb:
                for path in l:
                    size += numreads[path]
            # print('groups in batch', len(newb), size)
            sizes.append(size)

        return tmp_batches, sizes

    @staticmethod
    def balance_new(dic, numreads):
        batches = list(dic.values())
        newbatches = string.ascii_lowercase
        total_numer_of_reads = sum(numreads.values())
        avg = total_numer_of_reads // len(batches)
        groups = []
        for batch in batches:
            groups.extend([sublist for sublist in batch])

        count_per_group = []
        for group in groups:
            s = sum(numreads[path] for path in group)
            count_per_group.append(s)

        new = {}
        c = 0
        sizes = []
        while groups:
            batch_name = newbatches[c]
            batch = []
            size_of_batch = 0
            used = []
            for i, g in enumerate(groups):
                size = count_per_group[i]
                if size_of_batch + size > avg:
                    break
                else:
                    batch.append(g)
                    size_of_batch += size
                    used.append(i)

            used.sort(reverse=True)
            for index in used:
                del groups[index]
                del count_per_group[index]

            new[batch_name] = batch
            # print(batch_name, size_of_batch)
            sizes.append(size_of_batch)
            c += 1

        return new, sizes

if __name__ == '__main__':
    import pickle
    dic = pickle.load(open('/home/mmilanes/Downloads/sanefalcon/dic.pkl', 'rb'))
    numreads = pickle.load(open('/home/mmilanes/Downloads/sanefalcon/numreads.pkl', 'rb'))

    # Method 1
    old, sizes1 = Utils.balance_batches(dic, numreads)

    # Method 2
    kn, sizes2 = Utils.balance_w_knapsack(dic, numreads)

    # Method 3
    new, sizes3 = Utils.balance_new(dic, numreads)


    x = list(old.keys())
    plt.plot(x, sizes1, 'b', label='old')
    x = list(new.keys())
    plt.plot(x, sizes2, 'r', label='knap')
    x = list(new.keys())
    plt.plot(x, sizes3, 'g', label='new')
    plt.legend()
    plt.title('3 ways of balancing reads')
    plt.ylabel('# reads')
    plt.xlabel('batch')
    plt.savefig('balancing_batches.png')

