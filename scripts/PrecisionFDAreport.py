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
# @file    PrecisionFDAreport.py
#
# @brief   Generating the report for PrecisionFDA challenge
#
# @author  A-Tsai Su(chungtsai_su@atgenomix.com.tw)
#
# @date    2018/10/02
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

HELLO_WORLD = "Hello_World"
HEADER = "Reference_Genome_ID\tHello_World\tC01\tC02\tC03\tC04\tC05\tC06\tC07\tC08\tC09\tC10\tC11\tC12\tC13\tC14\t" \
              "C15\tC16\tC17\tC18\tC19\tC20\tC21"
NUM_SAMPLES = 21
UNKNOWN = "Unclassified"


def _usage():
    """ Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    """
    print(
        "PrecisionFDAreport.py -i <Input Folder>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <Input Folder>")

    print("Usage:")
    print("\tpython3 ./PrecisionFDAreport.py -i ../result/CR517/1 ")

    return


def generate_report(ifn, suffix):
    res = []
    # collect
    for i in range(NUM_SAMPLES+1):
        if i == 0:
            fn = "%s/%s_%s" % (ifn, HELLO_WORLD, suffix)
        else:
            fn = "%s/%02d_%s" % (ifn, i, suffix)
        ifd = open(fn, "r")
        idx = 0
        for line in ifd:
            if line.startswith("#"):
                continue
            items = line.strip().split("\t")
            if i == 0:
                res.append(line.strip())
            else:
                res[idx] += "\t" + items[1]
            if idx < 10:
                print("%d %d %s" % (i, idx, res[idx]))
            idx += 1
        ifd.close()

    # identify NTC and Positive Control
    for i in range(len(res)):
        items = str(res[i]).split("\t")
    # output
    ofd = open(ifn + "/" + suffix, "w")
    ofd.write("%s\n" % HEADER)
    for i in range(len(res)):
        if str(res[i]).startswith("CR_162"):
            continue
        ofd.write("%s\n" % res[i])
    ofd.close()

    return


def report(ifn):
    generate_report(ifn, "Genome_Identification.tsv")
    generate_report(ifn, "Genome_Quantification.tsv")

    return


def main(argv):
    infile = ""

    try:
        opts, args = getopt.getopt(argv, "hi:")
    except getopt.GetoptError:
        _usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            _usage()
            sys.exit()
        elif opt in "-i":
            infile = arg
    if not infile or not os.path.isdir(infile):
        print("Error: input folder(%s) is not existed" % infile)
        _usage()
        sys.exit(2)

    # Main Function
    report(infile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
