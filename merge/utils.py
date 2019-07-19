import os


class Utils:

    @staticmethod
    def sort_and_write(data, outfile):
        # data = list(set(data))  # remove duplicates and sort
        data.sort()
        to_be_written = map(lambda x: str(x) + "\n", data)
        with open(outfile, 'w') as out:
            out.writelines(to_be_written)

    @staticmethod
    def read_all_files(file_list):
        data = []
        for f in file_list:
            data.extend([int(line.strip()) for line in open(f, 'r')])
        return data

    @staticmethod
    def pick_manip_pairs(manips, list_):
        res = {}
        for subdir, manip_list in manips.items():
            for manip in manip_list:
                files = []
                manip_name = os.path.basename(manip)
                for path in list_:
                    base = os.path.basename(path)
                    if os.path.commonprefix([base, manip_name]) == manip_name:
                        files.append(path)
                res[manip] = files
        return res

    @staticmethod
    def prepare_jobs_arguments(files_to_merge):
        runs = []
        for subdir, dic in files_to_merge.items():
            for chrom in range(1, 23):
                fwd = dic['fwd'][chrom]
                rev = dic['rev'][chrom]
                files = fwd + rev
                runs.append((files, subdir, chrom))
        return runs