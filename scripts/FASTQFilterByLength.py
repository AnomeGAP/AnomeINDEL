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
# @file    FASTQFilterByLength.py
#
# @brief   Select sequence with minimal length from FASTQ file
#
# @author  Chung-Tsai Su(chungtsai_su@atgenomix.com)
#
# @date    2019/06/14
#
# @version 1.0
#
# @remark
#

import sys
import getopt

# CONSTANT
MIN_LENGTH = 500
BIN_SIZE = 500
BIN_NUM = 20
MAX_LENGTH = BIN_SIZE * BIN_NUM


def usage():
    print("zcat <Input FASTQ.gz> | FASTQFilterByLength.py -l <minimal length> -o <output FASTQ file>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-l: minimal length (Default: 500)")
    print("\t-o: output file")
    print("Usage:")
    print("\ttime gzcat ~/NA12878-novaseq/v1.0.2/result-1.0.2-qual-fix-2.fq.gz | python ~/src/github/AnomeINDEL/scripts/FASTQFilterByLength.py -l 1000 -o ~/result-1.0.2-qual-fix-2.1000.fq")
    print("\ttime zcat /seqslab/atsai/data/NA12878/result-1.0.2-qual-fix-2.fq.gz | python /seqslab/atsai/script/FASTQFilterByLength.py -l 1000 -o /seqslab/atsai/data/NA12878/result-1.0.2-qual-fix-2.1000.fq > /seqslab/atsai/data/NA12878/result-1.0.2-qual-fix-2.fq.gz.length-distribution.log")

    return


def analyzer(min_len, ofn):
    idx = 0
    total_bps = 0
    total_reads = 0
    pass_bps = 0
    pass_reads = 0
    enabled = False
    data = ""
    ofd = open(ofn, "w")
    a_reads = [0] * (BIN_NUM+1)
    a_bps = [0] * (BIN_NUM+1)

    for line in sys.stdin:
        if idx % 4 == 1:
            l = len(line.strip())
            if l >= min_len:
                pass_bps += l
                pass_reads += 1
                enabled = True
                data += line
            total_bps += l
            total_reads += 1
            if l >= MAX_LENGTH:
                for i in range(0, BIN_NUM+1):
                    a_reads[i] += 1
                    a_bps[i] += l
            else:
                for i in range(0, int(l / BIN_SIZE)+1):
                    a_reads[i] += 1
                    a_bps[i] += l
        elif idx % 4 == 0:
            if enabled:
                # ofd.write("%s" % data)
                enabled = False
            data = line
        elif enabled:
            data += line

        idx += 1

    print("Pass %d/%d reads with %d/%d bps" % (pass_reads, total_reads, pass_bps, total_bps))
    ofd.close()

    for i in range(BIN_NUM+1):
        print(">=%d\t%d\t%d" % (i*BIN_SIZE, a_reads[i], a_bps[i]))

    return


def main(argv):
    min_len = MIN_LENGTH
    ofile = "output.%d.fq" % min_len

    try:
        opts, args = getopt.getopt(argv, "hl:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt == '-l':
            min_len = int(arg)
            ofile = "output.%d.fq" % min_len
        elif opt == '-o':
            ofile = arg

    analyzer(min_len, ofile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
