import sys

binsize = 1000000
minShift = 4
threshold = 4
chrColumn = 2
startColumn = 3
sequenceColumn = 9


def retro(flag=True):
    """
    flag = True -> fwd files
    flag = False -> rev files
    """

    # Flush the current stack of reads
    def flush(read_buffer):
        stairSize = len(read_buffer)
        if stairSize <= threshold or threshold < 0:
            for line in read_buffer:
                if flag:  # .fwd files: only start position
                    print(line[0])
                else:  # .rev files: start + length of sequence
                    print(line[1])

    read_buffer = []
    prev_start = 0
    prev_chr = ''

    for line in sys.stdin:
        words = line.strip().split()
        cur_start, cur_chr = int(words[startColumn]), words[chrColumn]
        length = len(words[sequenceColumn]) - 1
        rev_l = cur_start + length  # (start + len sequence) for .rev files
        if not ((cur_chr == prev_chr) and (minShift >= cur_start - prev_start)):
            flush(read_buffer)
            read_buffer = []
        if len(read_buffer) == 0 or cur_start != read_buffer[-1][0]:
            # Normal ndups will be appended here
            read_buffer.append((cur_start, rev_l))

        prev_start = cur_start
        prev_chr = cur_chr

    flush(read_buffer)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-fwd', dest='feature', action='store_true')
    parser.add_argument('-rev', dest='feature', action='store_false')
    parser.set_defaults(feature=True)

    args = parser.parse_args()

    flag = args.feature
    retro(flag)

