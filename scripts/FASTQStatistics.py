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
# @file    FASTQStatistics.py
#
# @brief   Parsing FASTQ file to analyze the distribution of read length
#
# @author  Chung-Tsai Su(chungtsai_su@atgenomix.com)
#
# @date    2019/01/21
#
# @version 1.0
#
# @remark
#

import sys
import getopt

# CONSTANT
#ORI_LENGTH = 151
#UNIT = 100
#NUM_ELEMENT = 10
#UNIT = 1000
#NUM_ELEMENT = 50
ORI_LENGTH = 24
UNIT = 1
NUM_ELEMENT = 400

MAX_LENGTH = UNIT * NUM_ELEMENT


def usage():
    print("gzcat <Input FASTQ.gz> | FASTQStatistics.py > length_distribution.tsv")
    print("Argument:")
    print("\t-h: Usage")
    print("Usage:")
    print("\tgzcat ntuh-assembly-result.fastq.gz | python ~/src/github/AnomeINDEL/scripts/FASTQStatistics.py > ntu-assembly-result.fastq.tsv")
    print("\tgzcat ~/NA12878-novaseq/v1.0.2/result-1.0.2-qual.fq.gz | python ~/src/github/AnomeINDEL/scripts/FASTQStatistics.py > ~/NA12878-novaseq/v1.0.2/result-1.0.2-qual.tsv")
    print("\tgzcat ~/NA12878-novaseq/v1.0.2/result-qual-fix.fq.gz | python ~/src/github/AnomeINDEL/scripts/FASTQStatistics.py > ~/NA12878-novaseq/v1.0.2/result-qual-fix.fq.tsv")

    return


def analyzer():
    idx = 0
    cnt = 0
    a_len = [0] * (NUM_ELEMENT+1)
    max_length = 0
    total_bp = 0
    total_ori = 0

    for line in sys.stdin:
        if idx % 4 == 1:
            # print(line)
            l = len(line.strip())
            total_bp += l
            if max_length < l:
                max_length = l
            if l >= MAX_LENGTH:
                a_len[NUM_ELEMENT] += 1
            else:
                a_len[int(l/UNIT)] += 1
            if l <= ORI_LENGTH:
                total_ori += 1
            cnt += 1
        idx += 1

    print("#reads\t%d" % cnt)
    print("#bps\t%ld" % total_bp)
    print("average_legnth\t%f" % (float(total_bp/cnt)))
    print("max_length\t%d" % max_length)
    print("#singletons(<=%d)\t%d" % (ORI_LENGTH,total_ori))

    for i in range(NUM_ELEMENT+1):
        print("%d\t%d" % (i*UNIT, a_len[i]))

    return


def main(argv):

    try:
        opts, args = getopt.getopt(argv, "h")
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
