import os
import re
from collections import defaultdict

import nucleosome.nucleosome_detector as nd
from multi.multiprocess import launch_multiprocess as lp


class Tracker:
    def __init__(self, fm):
        self.fm = fm
        self.nucleosome_tracks = defaultdict(list)

    def create_nucleosome_files(self, folder, outfilestub, training=True):
        runs = []
        if training:
            files = [os.path.join(folder, f) for f in os.listdir(folder) if re.match('^anti', f)]
        else:
            files = [os.path.join(folder, f) for f in os.listdir(folder) if re.match('^merge', f)]
        for f in files:
            pattern = re.compile(r'\d{1,2}$')
            chrom = re.search(pattern, f).group()
            outfile = outfilestub + '.{}'.format(chrom)
            runs.append((folder, chrom, f, os.path.join(folder, outfile)))

        self.nucleosome_tracks[folder] = lp(runs, nd.create_nucleosome_file)

    def create_tracks(self, training=True):
        if training:
            to_exclude = [self.fm.config['folders']['profiles'], self.fm.config['folders']['rsp']]
            folders = [os.path.join(self.fm.trainfolder, o) for o in os.listdir(self.fm.trainfolder)
                       if os.path.isdir(os.path.join(self.fm.trainfolder, o)) and
                       not (os.path.join(self.fm.trainfolder, o) in to_exclude)]
            nucl_track_stub = self.fm.config['default']['nucltemplate'] + '_anti'
        else:
            folders = None  # FIXME where to put test merge files?
            nucl_track_stub = self.fm.config['default']['nucltemplate']

        for folder in folders:
            self.create_nucleosome_files(folder, nucl_track_stub, training)


