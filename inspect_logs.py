# Utility tool for tracking advancement by diving into logs/nucleosome.log
import os
import subprocess
from collections import Counter

trainingdir = '/home/mdlab/storage/sanefalcon/training'
logdir = './logs'
nucllog = os.path.join(logdir, 'nucleosome.log')

trainingsamples = 0
for root, subd, files in os.walk(trainingdir):
    for f in files:
        if f.endswith('.bam'):
            trainingsamples += 1

print('Samples total', trainingsamples)

# Nucleosome tracks
p = subprocess.Popen("grep 'saved' {}".format(nucllog), stdout=subprocess.PIPE, shell=True)
out, err = p.communicate()

if out:
    res = out.decode('utf-8').strip().split('\n')
    # res = [x.split()[14] for x in gres]
    subdirs = [row.partition('/home/mdlab/storage/sanefalcon/training')[2].split('/')[1] for row in res]
    c = Counter(subdirs)
    if all([x == 22 for x in c.items()]):
        print('tracks terminated')
    else:
        for k, v in sorted(c.items()):
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
    chroms = [int(row.partition(' phase for chrom')[2].strip().split()[0]) for row in res]
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

