import pickle
from collections import Counter
import time
from nucleosome.nucleosome_detect import prepare_chunk
import concurrent.futures

pool_size = 100


def work(list_of_center, positions):
    res = []
    for center in list_of_center:
        pos = positions.index(center)
        sublist = positions[pos - 200: pos + 200]
        tmp = list(filter(lambda x: x < 100, [abs(center - rs) for rs in sublist]))
        res.extend(tmp)
    return res


def get_profile(positions, nucleosome):
    distances = []
    s = time.time()
    i = 0
    for list_of_center in prepare_chunk(list(nucleosome.keys()), 100):
        with concurrent.futures.ThreadPoolExecutor(max_workers=pool_size) as executor:
            jobs = {}
            sub = executor.submit(work, list_of_center, positions)
            jobs[sub] = list_of_center

            for job in concurrent.futures.as_completed(jobs):
                res = job.result()
                distances.extend(res)

        print('list', i, time.time() - s)
        i += len(list_of_center)
        s = time.time()

    return Counter(distances)


if __name__ == "__main__":
    mergefile = '/data/tempff/anti.22'
    nucleosomefile = 'nucl.anti.22'
    positions = pickle.load(open(mergefile, 'rb'))
    nucleosome = pickle.load(open(nucleosomefile, 'rb'))
    print(len(nucleosome))
    profile = get_profile(positions, nucleosome)
    pickle.dump(profile, open('profile.22', 'wb'))

