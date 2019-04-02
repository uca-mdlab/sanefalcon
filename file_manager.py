import os
import string


class Utils:
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


class FileManager:
    def __init__(self, config):
        self.datafolder = config['default']['datafolder']
        self.trainfolder = config['default']['trainfolder']
        self.rspfolder = config['default']['rspfolder']
        self.nucl_file_template = config['default']['nucltemplate']
        self.anti_file_template = self.nucl_file_template + '_anti'
        self.bamlist = self.set_bamlist(config['default']['bamlist'])

    def list_files_to_use(self):
        if not self.bamlist:
            exit('Unable to find the list of bam files to use. Aborting.')
        list_files = [os.path.basename(x) for x in self.bamlist]
        files_to_link = []
        manip_list = []
        for root, subdir, files in os.walk(self.datafolder):
            for f in files:
                if f in list_files:
                    files_to_link.append(os.path.join(root, f))
                    manip_list.append(os.path.split(root)[1])

        return files_to_link, manip_list

    def set_bamlist(self, bamlistfile):
        """
        Get the list of bam files to use from the bam list file
        :param bamlistfile:
        :return: list of paths
        """
        return Utils.readfile(bamlistfile)

    def prepare_train_folder(self):
        letters = list(string.ascii_lowercase)
        files_to_link, manip_list = self.list_files_to_use()

        exit()
        batches = {}
        for num_batch, batch in enumerate(prepare_batches(manip_list)):
            batches[letters[num_batch]] = batch
            logger.debug("Batch {}: {}".format(letters[num_batch], batch))

        for batch_name, batch_list in batches.items():
            workingdir = os.path.join(traindir, batch_name)
            try:
                os.mkdir(workingdir)
                logger.debug("Created folder: {}".format(workingdir))
            except FileExistsError:
                logger.warning("Folder {} exists, skipping...".format(workingdir))

            for manip in batch_list:
                manip_regex = re.compile(manip)
                files = [f for f in files_to_link if re.search(manip_regex, f)]
                logger.debug("Batch {}: {}".format(batch_name, files))
                for fname in files:
                    run = fname.split(bamdir)[1].split('/')[-2]
                    runpath = os.path.join(workingdir, run)
                    if not os.path.isdir(runpath):
                        os.mkdir(runpath)
                    os.symlink(fname, os.path.join(runpath, fname.split('/')[-1]))
                    fname = fname + ".bai"
                    os.symlink(fname, os.path.join(runpath, fname.split('/')[-1]))
            logger.info("Batches created with symlinks")


if __name__ == '__main__':
    import configparser

    config = configparser.ConfigParser()
    config.read('my.conf')
    f = FileManager(config)
    f.prepare_train_folder()
