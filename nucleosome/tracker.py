import os
import re
from collections import defaultdict

import nucleosome.nucleosome_detector as nd
from multi.multiprocess import launch_multiprocess as lp


class Tracker:
    def __init__(self, fm):
        self.fm = fm
        self.nucleosome_tracks = defaultdict(list)

    def build_runs(self, folder, outfilestub, training):
        runs = []
        if training:
            files = [os.path.join(folder, f) for f in os.listdir(folder) if re.match('^anti', f)]
        else:
            files = [os.path.join(folder, f) for f in os.listdir(folder) if re.match('^merge', f)]
        for f in files:
            pattern = re.compile(r'\d{1,2}$')
            chrom = re.search(pattern, f).group()
            outfile = outfilestub + '.{}'.format(chrom)
            runs.append((chrom, f, os.path.join(folder, outfile)))
        return runs

    def create_nucleosome_files(self, folder, outfilestub, training=True):
        runs = self.build_runs(folder, outfilestub, training)
        self.nucleosome_tracks[folder] = lp(runs, nd.create_nucleosome_file)

    def create_tracks(self, training=True):
        if training:
            to_exclude = [self.fm.config['folders']['profiles'], self.fm.config['folders']['rsp']]
            folders = [os.path.join(self.fm.trainfolder, o) for o in os.listdir(self.fm.trainfolder)
                       if os.path.isdir(os.path.join(self.fm.trainfolder, o)) and
                       not (os.path.join(self.fm.trainfolder, o) in to_exclude)]
            nucl_track_stub = self.fm.config['default']['nucltemplate'] + '_anti'
        else:
            folders = [self.fm.testfolder]
            nucl_track_stub = self.fm.config['default']['nucltemplate']

        for folder in folders:
            self.create_nucleosome_files(folder, nucl_track_stub, training)

    def get_data(self, mapping):
        pack = defaultdict(dict)
        for subdir, dic in mapping.items():
            for name, files in dic.items():
                fwd = [f for f in files if f.endswith('.fwd')]
                rev = [f for f in files if f.endswith('.rev')]
                for chrom in range(1, 23):
                    pattern = re.compile(r'\.{}\.'.format(chrom))
                    nucl_pattern = re.compile(r'\.{}$'.format(chrom))
                    nucl_track = [f for f in self.nucleosome_tracks[subdir] if re.search(nucl_pattern, f)][0]
                    fwd_c = list(filter(lambda x: re.search(pattern, x), fwd))
                    rev_c = list(filter(lambda x: re.search(pattern, x), rev))
                    pack[chrom].update({subdir: {'fwd': fwd_c, 'rev': rev_c, 'nucl_file': nucl_track}})
        return pack


if __name__ == '__main__':
    import configparser
    from manager.file_manager import FileManager

    config = configparser.ConfigParser()
    config.read('../sanefalcon.conf')
    filemanager = FileManager(config)
    t = Tracker(filemanager)
    nucl_track_stub = t.fm.config['default']['nucltemplate'] + '_anti'
    folder = os.path.join(t.fm.trainfolder, 'a')
    runs = t.build_runs(folder, nucl_track_stub, training=True)
    print(len(runs))
    print(runs[0])


