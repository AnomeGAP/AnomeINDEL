#!/bin/env python
#
# @note Copyright (C) 2018, Atgenomix Incorporated. All Rights Reserved.
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
# @file    FASTQanalyzer.py
#
# @brief   Parsing FASTQ file to analyze the distribution of bad quality reads
#
# @author  Chung-Tsai Su(chungtsai_su@atgenomix.com)
#
# @date    2018/12/24
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os.path
import re
from collections import defaultdict
import fileinput

# Default Parameter
MAX_LENGTH = 151
MIN_QUALITY = 20


def usage():
    print("zcat <Input FASTQ> | FASTQanalyzer.py > quality_distribution.tsv")
    print("Argument:")
    print("\t-h: Usage")
    print("Usage:")
    print("\tzcat NTUH-PBMC-01_S99_L999_R1_001.fastq.gz | python ./FASTQanalyzer.py > NTUH-PBMC-01_R1.qual.tsv")
    print("\tzcat NTUH-PBMC-01_S99_L999_R2_001.fastq.gz | python ./FASTQanalyzer.py > NTUH-PBMC-01_R2.qual.tsv")

    return


def analyzer():
    idx = 0
    quality = [0] * MAX_LENGTH

    for line in sys.stdin:
        if idx % 4 == 3:
            # print(line)
            cnt = 0
            for i in line.rstrip():
                q = ord(i) - ord('!')
                if q < MIN_QUALITY:
                    cnt += 1
                    # print("%s\t%d\t%d" % (i, q, cnt))
            quality[cnt] += 1
            # break
        idx += 1

    s = 0
    for i in range(MAX_LENGTH):
        s += quality[i]
        print("%d\t%d\t%d" % (i, quality[i], s))

    return


def main(argv):

    try:
        opts, args = getopt.getopt(argv, "hi")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()

    analyzer()

    return


if __name__ == '__main__':
    main(sys.argv[1:])
