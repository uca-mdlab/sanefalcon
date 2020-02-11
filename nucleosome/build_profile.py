import pickle
from collections import Counter
import time
from nucleosome.nucleosome_detect import prepare_chunk
import concurrent.futures
from log_setup import setup_logger

pool_size = 50

logger = setup_logger('nucl-profile', 'logs/nucldetect.log')


def work(list_of_center, positions):
    res = []
    for center in list_of_center:
        pos = positions.index(center)
        sublist = positions[pos - 110: pos + 110]
        tmp = list(filter(lambda x: x < 100, [abs(center - rs) for rs in sublist]))
        res.extend(tmp)
    return res


def get_profile(positions, nucleosome):
    logger.debug(f'Found {len(nucleosome)} centers to loop on...')
    distances = []
    s0 = time.time()
    s = s0
    i = 0
    for list_of_center in prepare_chunk(list(nucleosome.keys()), pool_size):
        with concurrent.futures.ThreadPoolExecutor(max_workers=pool_size) as executor:
            jobs = {}
            sub = executor.submit(work, list_of_center, positions)
            jobs[sub] = list_of_center

            for job in concurrent.futures.as_completed(jobs):
                res = job.result()
                distances.extend(res)

        if i % 1000 == 0:
            logger.debug(f'{i} centers profiled, time: {round(time.time() - s)} sec.')
            s = time.time()

        i += len(list_of_center)

    profile = Counter(distances)
    logger.debug(f'Profile computed time: {round(time.time() - s0)} sec.')
    return profile


if __name__ == "__main__":
    mergefile = '/data/tempff/anti.22'
    nucleosomefile = 'nucl.anti.22'
    positions = pickle.load(open(mergefile, 'rb'))
    nucleosome = pickle.load(open(nucleosomefile, 'rb'))
    # print(len(nucleosome))
    profile = get_profile(positions, nucleosome)
    pickle.dump(profile, open('profile.22', 'wb'))

