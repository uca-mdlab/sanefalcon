# Utility tool for tracking advancement by diving into logs/nucleosome.log
import os
import subprocess
from collections import Counter
import glob


trainingdir = '/home/mdlab/storage/sanefalcon/training'
profiledir = os.path.join(trainingdir, 'profiles')
logdir = './logs'

genlog = os.path.join(logdir, 'sanefalcon.log')
nucllog = os.path.join(logdir, 'nucleosome.log')

trainingsamples = 0
p = subprocess.Popen("grep 'TRAINING' {}".format(genlog), stdout=subprocess.PIPE, shell=True)
out, err = p.communicate()
res = out.decode('utf-8').strip().split()[8:10]
name, trainingsamples = res[0].strip(','), int(res[1].strip(','))
expected_profile_count = trainingsamples * 4 * 22

print('Run: ', name)
print('Training samples : ', trainingsamples)

# Nucleosome tracks
p = subprocess.Popen("grep 'saved' {}".format(nucllog), stdout=subprocess.PIPE, shell=True)
out, err = p.communicate()

if out:
    res = out.decode('utf-8').strip().split('\n')
    # res = [x.split()[14] for x in gres]
    subdirs = [row.partition('/home/mdlab/storage/sanefalcon/training')[2].split('/')[1] for row in res]
    subdirs.sort()
    c = Counter(subdirs)
    if all([x == 22 for x in c.items()]):
        print('tracks terminated')
    else:
        for k, v in c.items():
            if v == 22:
                print(k, 'Done')
            else:
                print(k, v)
else:
    print('Tracks completed')

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
        for k, v in c.items():
            if v == trainingsamples:
                print('chrom', k, 'Done')
            else:
                print('chrom', k, ' - count: ', v)

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
        for k, v in sorted(c.items()):
            if v == trainingsamples:
                print('chrom', k, 'Done')
            else:
                print('chrom', k, ' - count: ', v)


profiles = glob.glob("{}/*.*".format(profiledir))
print('Profiles saved: {0:.2f}%'.format((len(profiles) / expected_profile_count)*100))
