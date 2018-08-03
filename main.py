import os
import argparse
import string
import re

import time
import datetime

import random
import logging
from prepare_folders import prepare_train_folder

logging.basicConfig(format='%(asctime)s - %(levelname)s: %(message)s',
                    filename='sanefalcon.log', filemode='w', level=logging.DEBUG)

logger = logging.getLogger(__name__)


def retro_function(binsize=1000000, retdist=4, retthres=4, colchr=2, colpos=3):
    testline = "5J1WT:06455:14655       0       chr1    10034   6       66M1D20M1D9M2S  *       0       0       CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTACCCCTAACCTAACCCAACCCTAACCCTACCCTAACCTA    98*8<6;:379/55*55/5<3;<69;1//)//)/95::2:9)55*5/)/80/.).8).....(4:1728:5::37279/665;;4;:990))()())        ZP:B:f,0.00788558,0.00494094,0.000700576        ZM:B:s,280,-16,240,-16,-10,268,-32,258,212,490,238,492,-8,58,262,248,52,274,32,236,236,256,220,192,228,-20,24,-4,796,40,-30,42,260,456,680,42,246,388,648,22,228,10,70,-4,440,10,766,12,254,50,2,450,14,770,18,46,238,0,14,386,852,62,18,114,248,518,778,-12,236,494,844,-24,210,0,-28,4,592,58,776,36,220,52,16,588,0,768,-22,86,170,-42,-18,282,816,22,-16,132,248,508,616,-6,224,494,720,30,116,74,-30,-44,462,36,776,-8,176,124,40,440,-20,646,40,104,212,44,-20,232,738,22,64,78,216,390,492,8,146,376,616,46,98,10,-96,-58,528,92,726,76,204,72,52,374,10   ZF:i:28 MD:Z:59A6^C20^A9        NM:i:3  AS:i:77 XM:i:97 XA:Z:map4-1     XS:i:73 RG:Z:5J1WT.IonXpress_065      PG:Z:tmap"

    chr_names = ['X', 'Y']
    chr_names.extend([str(x) for x in range(1, 23)])
    chromosomes = dict.fromkeys(chr_names, [0])

    prevWords = ['0'] * 10
    readBuff = []
    fullBuff = []

    curWords = testline.split()
    print(len(curWords))
    if not ((curWords[colchr] == prevWords[colchr]) and (
            retdist >= (int(curWords[colpos]) - int(prevWords[colpos])))):
        readBuff = []
        fullBuff = []

    if len(readBuff) == 0 or int(curWords[colpos]) != readBuff[-1][1]:
        # Normal ndups will be appended here
        readBuff.append([curWords[colchr], int(curWords[colpos])])
        fullBuff.append(testline)

    stairSize = len(readBuff)
    if stairSize <= retthres or retthres < 0:
        for line in fullBuff:
            print(line.strip())
    print(readBuff)


if __name__ == "__main__":
    # parser = argparse.ArgumentParser(description='Prepare the environment and start the plugin',
    #                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #
    # parser.add_argument('bamdir', type=str, default="/results/analysis/output/Home", help='path of the manipulation')
    # parser.add_argument('traindir', type=str, default="/tmp/sanefalcontrain2", help='path of the train subtree')
    #
    # args = parser.parse_args()
    # bamdir = args.bamdir
    # traindir = args.traindir
    #
    # prepare_train_folder(bamdir, traindir)
    # print(traindir)
    retro_function()