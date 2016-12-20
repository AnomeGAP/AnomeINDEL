#!/bin/env python3
#
# @note Copyright (C) 2016, Anome Incorporated. All Rights Reserved.
#       This program is an unpublished copyrighted work which is proprietary to
#       Anome Incorporated and contains confidential information that is not to
#       be reproduced or disclosed to any other person or entity without prior
#       written consent from Anome, Inc. in each and every instance.
#
# @warning Unauthorized reproduction of this program as well as unauthorized
#          preparation of derivative works based upon the program or distribution of
#          copies by sale, rental, lease or lending are violations of federal copyright
#          laws and state trade secret laws, punishable by civil and criminal penalties.
#
# @file    CalErrorCorrectionRate.py
#
# @brief   An calculator for Error Correction Rate
#
# @author  A-Tsai Su(chungtsai_su@anome.com.tw)
#
# @date    2016/12/19
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os.path
import re
import random
from collections import defaultdict

# Parameter setting


def _usage():
    """ Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    """
    print(
        "CalErrorCorrectionRate.py -i <FASTQ file> -a <Answer file> -o <Output file>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <FASTQ file>")
    print("\t-a: <Answer file>")
    print("\t-o: <Output file>")

    print("Usage:")
    print("\tpython3 ~/src/github/AnomeINDEL/scripts/CalErrorCorrectionRate.py -i ~/data/hg19/chr22.corrected.pp.fq "
          "-a ~/data/hg19/chr22_error.tsv -o ~/data/hg19/chr22.corrected.pp.fq.ans")
    print("\tpython3 ~/src/github/AnomeINDEL/scripts/CalErrorCorrectionRate.py -i ~/data/hg19/chr22.corrected.pp.fq.2k "
          "-a ~/data/hg19/chr22_error.tsv.20 -o ~/data/hg19/chr22.corrected.pp.fq.ans")

    return


def cal_ec(ifn, afn, ofn):
    ofd = open(ofn, "w")
    afd = open(afn, "r")
    ifd = open(ifn, "r")
    total = 0
    num_checked = 0
    num_corrected = 0
    vhash = defaultdict(list)

    for line in afd:
        if re.match("#", line):
            continue
        items = line.strip().split("\t")
        vhash[items[0]+"\t"+items[1]].append(items[2] + ":" + items[3] + ":" + items[4] + ":" + items[5])
        total += 1
        #print("%s:%s" % (items[0], items[1]))
        #print("%s:%s:%s:%s" % (items[2], items[3], items[4], items[5]))

    # write header
    ofd.write("#ID\t#pair\tPOS\tREF\tALT\tREF_POS\tCorrected\n")
    b_seq = False
    id = 0

    for line in ifd:
        if re.match("^@", line):
            b_seq = True
            items = line.strip().split("/")
            id = items[0][5:] + "\t" + items[1]
            #print("%s" % id)
            continue
        if not b_seq:
            continue

        if id in vhash:
            seq = list(line.strip())
            #print("%s" % seq)
            for val in vhash[id]:
                #print("%s" % val)
                items = val.split(":")
                if seq[int(items[0])] == items[1]:
                    num_corrected += 1
                    # #ID\t#pair\tPOS\tREF\tALT\tREF_POS\tCorrected\n
                    ofd.write("%s\t%s\t%s\t%s\t%s\t%d\n" % (id, items[0], items[1], items[2], items[3], 1))
                else:
                    ofd.write("%s\t%s\t%s\t%s\t%s\t%d\n" % (id, items[0], items[1], items[2], items[3], 0))

                num_checked += 1

        b_seq = False

    print("Total=%d" % total)
    print("Corrected=%d" % num_corrected)
    print("Checked=%d" % num_checked)
    print("Error Correction Rate = %.2f%%" % (num_corrected*100/num_checked))

    ifd.close()
    afd.close()
    ofd.close()

    return


def main(argv):
    infile = ""
    ansfile = ""
    outfile = ""

    try:
        opts, args = getopt.getopt(argv, "hi:a:o:")
    except getopt.GetoptError:
        _usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            _usage()
            sys.exit()
        elif opt in "-i":
            infile = arg
            outfile = infile + ".out"
        elif opt in "-a":
            ansfile = arg
        elif opt in "-o":
            outfile = arg
    if not infile or not os.path.isfile(infile):
        print("Error: input file(%s) is not existed" % (infile))
        _usage()
        sys.exit(2)
    elif not ansfile or not os.path.isfile(ansfile):
        print("Error: input file(%s) is not existed" % (ansfile))
        _usage()
        sys.exit(3)

    # Main Function
    cal_ec(infile, ansfile, outfile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
