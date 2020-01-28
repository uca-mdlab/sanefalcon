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


p = subprocess.Popen("grep 'TRAINING' {}".format(genlog), stdout=subprocess.PIPE, shell=True)
out, err = p.communicate()
res = out.decode('utf-8').strip().split()[8:10]
name, trainingsamples = res[0].strip(','), int(res[1].strip(','))


p = subprocess.Popen("grep 'TESTING' {}".format(genlog), stdout=subprocess.PIPE, shell=True)
out, err = p.communicate()
res = out.decode('utf-8').strip().split()[9]
testingsamples = int(res.strip(','))

print('Run: ', name)
print('Training samples : ', trainingsamples)
print('Testing samples : ', testingsamples)
print('----')

expected_training_profiles = trainingsamples * 4 * 22
expected_testing_profiles = testingsamples * 4 * 22
p = subprocess.Popen("tail -1 {}".format(multilog), stdout=subprocess.PIPE, shell=True)
out, err = p.communicate()
res = out.decode('utf-8').strip()
if re.search('runs_left', res):
    print('Multi - Last command: ', ' '.join(res.split(':')[-2:]).strip())

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
        done = 0
        for k, v in c.items():
            if v == 22:
                done += 1
            else:
                print(k, v)
        print('Done: ', done)


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
        print('Forward complete')
    else:
        print('Forward')
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
        print('Reverse complete')
    else:
        print('Reverse')
        done = 0
        for k, v in sorted(c.items()):
            if v == trainingsamples:
                done += 1
            else:
                print('chrom', k, ' - count: ', v)
        print('Done:', done)

profiles = glob.glob("{}/*.*".format(profiledir))
print('Profiles saved: {0:.2f}%'.format((len(profiles) / expected_training_profiles) * 100))

# grep - e 'saved\ .*\/testing/'