import nucleosome.get_profile_parallel as gpp


class Profiler:
    def __init__(self, fm, tracks):
        self.fm = fm
        self.tracks = tracks

    def compute_profiles(self):
        data = gpp.get_data(self.fm)