#!/bin/env python3
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
# @file    CalErrorDistribution.py
#
# @brief   Calculate the distribution of sequencing errors generated by PairReadGenerator.py
#
# @author  A-Tsai Su(chungtsai_su@atgenomix.com)
#
# @date    2018/08/17
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os.path
import re

# Parameter setting
READ_LENGTH = 151


def _usage():
    """ Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    """
    print(
        "CalErrorDistribution.py -i <Error file> -o <Output file> -l <read length> ")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <Error file> (TSV format)")
    print("\t-o: <Output file>")
    print("\t-l: <read length>")

    print("Usage:")
    print("\tpython3 ~/src/github/AnomeINDEL/scripts/CalErrorDistribution.py "
          "-i ~/data/hg19/chr22_error.tsv "
          "-o ~/data/hg19/chr22_error.dist "
          "-l 151")

    return


def cal_error_distribution(ifn, ofn):
    ofd = open(ofn + "_error.tsv", "w")
    ifd = open(ifn, "r")
    lcount = [0 for x in range(READ_LENGTH)]
    icount = [0 for x in range(READ_LENGTH)]
    dcount = [0 for x in range(READ_LENGTH)]

    cnt = 0
    cnt_ins = 0
    cnt_del = 0

    # write header
    ofd.write("#POS\t#Substitutes\t#Insertions\t#Deletions\n")
    for line in ifd:
        if re.match("#", line):
            continue
        items = line.split('\t')
        if items[2] == '-': # insertion
            cnt_ins += 1
            icount[int(items[1])] += 1
        elif items[3] == '-':
            cnt_del += 1
            dcount[int(items[1])] += 1
        else:
            lcount[int(items[1])] += 1
            cnt += 1
        if cnt % 1000 == 0:
            print("%d errors are processed" % cnt)

    for i in range(READ_LENGTH):
        ofd.write("%d\t%d\t%d\t%d\n" % (i, lcount[i], icount[i], dcount[i]))
    ofd.write("Total\t%d\t%d\t%d\n" % (cnt, cnt_ins, cnt_del))

    print("Total: %d substitutes, %d insertions and %d deletions" % (cnt, cnt_ins, cnt_del))
    ifd.close()
    ofd.close()

    return


def main(argv):
    infile = ""
    outfile = ""
    global READ_LENGTH

    try:
        opts, args = getopt.getopt(argv, "hi:o:l:")
    except getopt.GetoptError:
        _usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            _usage()
            sys.exit()
        elif opt in "-i":
            infile = arg
            outfile = infile + ".dist"
        elif opt in "-o":
            outfile = arg
        elif opt in "-l":
            READ_LENGTH = int(arg)
    if not infile or not os.path.isfile(infile):
        print("Error: input file(%s) is not existed" % (infile))
        _usage()
        sys.exit(2)

    # Main Function
    cal_error_distribution(infile, outfile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
