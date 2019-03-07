import re
from collections import defaultdict
import statistics


fname = 'nucleosome.log'

pattern = re.compile(r"DEBUG:\[p:[0-9]{5} f:")
fnames = ['merge.{}'.format(x) for x in range(1, 23)]

with open(fname, 'r') as infile:
    debuglines = [l for l in infile.readlines() if re.search(pattern, l)]

times = defaultdict(float)
for l in debuglines:
    f, val = l.split('f:')[1].split(']')
    times[f] += float(val.split('time = ')[1])


# for k, v in times.items():
#     print(k, '{:.1f}'.format(v))


v = list(times.values())
print('min {0:.2f}'.format(min(v) / 60.0))
print('max {0:.2f}'.format(max(v) / 60.0))
print('mean {0:.2f}'.format(statistics.mean(v) / 60.0))
print('std {0:.2f}'.format(statistics.stdev(v) / 60.0))
