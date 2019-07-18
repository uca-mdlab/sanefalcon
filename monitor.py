import configparser

from manager.file_manager import FileManager
from manager.utils import Utils
from manager.rsp import RspBuilder

if __name__ == '__main__':
    config = configparser.ConfigParser()
    config.read('local.conf')
    f = FileManager(config)
    f.check_paths()
    # print(f.get_data())
    # print(f.get_data(mask='.bam$'))
    bamlist = Utils.readfile(config['training']['bamlist'])
    batchsize = int(config['training']['batchsize'])
    batches = f.prepare_train_folder(bamlist, batchsize)