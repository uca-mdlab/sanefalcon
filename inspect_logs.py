# Utility tool for tracking advancement by diving into logs/nucleosome.log
import os
import subprocess
from collections import Counter
import glob
import re


trainingdir = '/home/mdlab/storage/sanefalcon/training'
profiledir = os.path.join(trainingdir, 'profiles')
logdir = './logs'

genlog = os.path.join(logdir, 'sanefalcon.log')
nucllog = os.path.join(logdir, 'nucleosome.log')
multilog = os.path.join(logdir, 'multi.log')
managerlog = os.path.join(logdir, 'manager.log')


def launch_grep_on_file(hook, filename):
    cmd = f"grep -e {hook} {filename}"
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    out, err = p.communicate()
    if out:
        return out.decode('utf-8').strip()
    else:
        return None


def get_app_summary():
    name, trainingsamples, testingsamples = None, None, None
    hook = 'TRAINING'
    res = launch_grep_on_file(hook, genlog)
    try:
        tmp = res.split()[8:10]
        name, trainingsamples = tmp[0].strip(','), int(tmp[1].strip(','))
    except:
        print(f'Found nothing on {hook} {genlog}')

    hook = 'TESTING'
    res = launch_grep_on_file(hook, genlog)
    try:
        tmp = res.split()[9]
        testingsamples = int(tmp.strip(','))
    except:
        print(f'Found nothing on {hook} {genlog}')

    return name, trainingsamples, testingsamples


def get_batches():
    hook = 'Batch\ .'
    res = launch_grep_on_file(hook, managerlog)
    tmp = [x[x.index('Batch'):] for x in res.split('\n')]
    num_batches = len(tmp)
    for row in tmp:
        batch, runs, num_samples = row.replace(' - ', '').split(',')
        print(batch, runs, num_samples)
    return num_batches


name, trainingsamples, testingsamples = get_app_summary()
print('Run: ', name)
print('Training samples : ', trainingsamples)
print('Testing samples : ', testingsamples)
print()
num_batches = get_batches()
print('Num batches : ', num_batches)
print('----')


expected_training_profiles = trainingsamples * 4 * 22
expected_testing_profiles = testingsamples * 4 * 22
p = subprocess.Popen("tail -1 {}".format(multilog), stdout=subprocess.PIPE, shell=True)
out, err = p.communicate()
res = out.decode('utf-8').strip()
if re.search('runs_left', res):
    print('Multi - Last command: ', ' '.join(res.split(':')[-2:]).strip())
print('----')

# Nucleosome tracks
p = subprocess.Popen("grep -e 'saved\ .*\/training/' {}".format(nucllog), stdout=subprocess.PIPE, shell=True)
out, err = p.communicate()

if out:
    res = out.decode('utf-8').strip().split('\n')
    subdirs = [row.partition(trainingdir)[2].split('/')[1] for row in res]
    subdirs.sort()
    c = Counter(subdirs)
    if all([x == 22 for x in c.items()]):
        print('Nucleosome Tracks terminated')
    else:
        print('Nucleosome Tracks:')
        done = 0
        for k, v in c.items():
            if v == 22:
                done += 1
            else:
                print(f'Batch {k}: tracks saved: {v}')
        print(f"Batches complete: {done}/{num_batches}")


# Nucleosome Profiles Forward
p = subprocess.Popen("grep 'End of forward' {}".format(nucllog), stdout=subprocess.PIPE, shell=True)
out, err = p.communicate()
if out:
    res = out.decode('utf-8').strip().split('\n')
    chroms = [row.partition(' phase for chrom')[2].strip().split()[0] for row in res]
    chroms = [int(x) for x in chroms]
    chroms.sort()
    c = Counter(chroms)
    if all([x == trainingsamples for x in c.items()]):
        print('Nucleosome profiles: Forward complete')
    else:
        print('Nucleosome profiles: Forward')
        done = 0
        for k, v in c.items():
            if v == trainingsamples:
                done += 1
            else:
                print('chrom', k, ' - count: ', v)
        print('Done:', done)

# Nucleosome Profiles Reverse
p = subprocess.Popen("grep 'End of reverse' {}".format(nucllog), stdout=subprocess.PIPE, shell=True)
out, err = p.communicate()
if out:
    res = out.decode('utf-8').strip().split('\n')
    chroms = [row.partition(' phase for chrom')[2].strip().split()[0] for row in res]
    chroms = [int(x) for x in chroms]
    chroms.sort()
    c = Counter(chroms)
    if all([x == trainingsamples for x in c.items()]):
        print('Nucleosome profiles: Reverse complete')
    else:
        print('Nucleosome profiles: Reverse')
        done = 0
        for k, v in sorted(c.items()):
            if v == trainingsamples:
                done += 1
            else:
                print('chrom', k, ' - count: ', v)
        print('Done:', done)

print('----')
profiles = glob.glob("{}/*.*".format(profiledir))
print('Profiles saved: {0:.2f}%'.format((len(profiles) / expected_training_profiles) * 100))

# grep - e 'saved\ .*\/testing/'