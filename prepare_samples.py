# old: prepSamples.sh + retro.py

import configparser
import os
import subprocess
import sys
import logging


logging.basicConfig(format='%(asctime)s - %(levelname)s: %(message)s',
                    filename='sanefalcon.log', filemode='w', level=logging.DEBUG)

logger = logging.getLogger("prepare_samples")


def prepare_samples(datafolder, trainfolder, samtools):
    chromosomes = range(1, 23)

    for root, subdir, files in os.walk(datafolder):
        for fname in files:
            if fname.endswith('.bam'):
                bamfile = os.path.join(root, fname)
                outfile_name = fname.split('.')[0] + '.sort.bam'
                for chrom in chromosomes:
                    outfile = os.path.join(trainfolder, outfile_name + '.{}.start.fwd'.format(chrom))
                    if not os.path.isfile(outfile):
                        one = [samtools, "view", bamfile, "chr{}".format(chrom), "-F", "20", "-q" "1"]
                        two = [sys.executable, "retro.py"]
                        three = ["awk", "{print $4}"]
                        p1 = subprocess.Popen(one, stdout=subprocess.PIPE)
                        p2 = subprocess.Popen(two, stdin=p1.stdout, stdout=subprocess.PIPE)
                        p3 = subprocess.Popen(three, stdin=p2.stdout, stdout=open(outfile, 'w'))
                        output = p3.communicate()[0]

                        revoutfile = os.path.join(trainfolder, outfile_name + '.{}.start.rev'.format(chrom))
                        one = [samtools, "view", bamfile, "chr{}".format(chrom), "-f", "16", "-F", "4", "-q" "1"]
                        three = ["awk", "{print ($4 + length($10) - 1)}"]
                        p1 = subprocess.Popen(one, stdout=subprocess.PIPE)
                        p2 = subprocess.Popen(two, stdin=p1.stdout, stdout=subprocess.PIPE)
                        p3 = subprocess.Popen(three, stdin=p2.stdout, stdout=open(revoutfile, 'w'))
                        output = p3.communicate()[0]

                logger.debug("prep_samples for {}".format(fname))


if __name__ == "__main__":
    conf_file = 'sanefalcon.conf'
    config = configparser.ConfigParser()
    config.read(conf_file)

    samtools = config['default']['samtools']
    datafolder = config['default']['datafolder']
    trainfolder = config['default']['trainfolder']
    nucleosomefolder = config['default']['nucleosomefolder']

    prepare_samples(datafolder, trainfolder, samtools)
