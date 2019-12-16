#!/bin/env python
#
# @note Copyright (C) 2019, Atgenomix Incorporated. All Rights Reserved.
#       This program is an unpublished copyrighted work which is proprietary to
#       Atgenomix Incorporated and contains confidential information that is not to
#       be reproduced or disclosed to any other person or entity without prior
#       written consent from Atgenomix, Inc. in each and every instance.
#
# @warning Unauthorized reproduction of this program as well as unauthorized
#          preparation of derivative works based upon the program or distribution of
#          copies by sale, rental, lease or lending are violations of federal copyright
#          laws and state trade secret laws, punishable by civil and criminal penalties.
#
# @file    FASTQTrimmer.py
#
# @brief   Trimming barcodes from IonProton's data
#
# @author  Chung-Tsai Su(chungtsai_su@atgenomix.com)
#
# @date    2019/12/16
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os
import gzip
from collections import defaultdict

# CONSTANT
NUM_TRIMMING = 25
MIN_LENGTH = 60
MAX_LENGTH = 180


def usage():
    print("FASTQTrimmer.py -i <Input FASTQ.gz file > -t <Number of the first bps trimmed>"
          " -l <minimal length> -o <output FASTQ.gz file>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: input file (fq.gz)")
    print("\t-t: Number of the first bps trimmed (Default: %d)" % NUM_TRIMMING)
    print("\t-l: minimal length (Default: %d)" % MIN_LENGTH)
    print("\t-m: maximal length (Default: %d)" % MAX_LENGTH)
    print("\t-o: output file")
    print("Usage:")
    print("\tpython ~/src/github/AnomeINDEL/scripts/FASTQTrimmer.py -i ./input.fq.gz -t 25 -l 60 -m 180 -o output.fq.gz")
    print("\tpython ./FASTQTrimmer.py -i ../data/YG/user_YGProton1-883-chip22_PLPK_mdel_mdb1.fq.gz -t 25 -l 60 -m 180 "
          " -o ../data/YG/user_YGProton1-883-chip22_PLPK_mdel_mdb1.t25.l60.m180.fq.gz")

    return


def trimmer(ifile, num_trimming, min_len, max_len, ofile):
    ofd = gzip.open(ofile, "wb")
    ifd = gzip.open(ifile, "rt")
    buf = ""
    idx = 0
    skipped = False
    total = 0
    cnt = 0

    for line in ifd:
        if idx % 4 == 0:
            buf = line
            total += 1
        elif idx % 4 == 1:
            if len(line) < num_trimming + min_len + 1 or len(line) > num_trimming + max_len + 1:
                skipped = True
            else:
                skipped = False
                cnt += 1
                buf += line[num_trimming:]
        elif not skipped:
            if idx % 4 == 2:
                buf += line
            else:
                buf += line[num_trimming:]
                ofd.write(buf)
        idx += 1

    ifd.close()
    ofd.close()
    print("%d / %d reads are passed" % (cnt, total))
    return


def main(argv):
    num_trimming = NUM_TRIMMING
    min_len = MIN_LENGTH
    max_len = MAX_LENGTH
    ofile = ""
    ifile = ""

    try:
        opts, args = getopt.getopt(argv, "hi:t:l:m:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt == '-i':
            ifile = arg
        elif opt == '-t':
            num_trimming = int(arg)
        elif opt == '-l':
            min_len = int(arg)
        elif opt == '-m':
            max_len = int(arg)
        elif opt == '-o':
            ofile = arg

    if ifile == "":
        print("Error: '-i' is required")
        usage()
        sys.exit(2)
    elif not os.path.isfile(ifile):
        print("Error: input file(%s) is not existed" % ifile)
        usage()
        sys.exit(3)

    if ofile == "":
        ofile = "output.t%d.l%d.fq.gz" % (num_trimming, min_len)

    trimmer(ifile, num_trimming, min_len, max_len, ofile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
