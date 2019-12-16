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
# @file    LengthFilter.py
#
# @brief   FASTA filter by Length
#
# @author  Chung-Tsai Su(chungtsai.su@atgenomix.com)
#
# @date    2019/07/08
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os.path

# CONST

# Default Parameter
SHORT_THRESHOLD = 200
LONG_THRESHOLD = 1000


def usage():
    print("LengthFilter.py -i <Input FASTA file>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: Input BLAST file")
    print("Usage:")
    print("\tpython ./LengthFilter.py -i ~/NA12878-novaseq/v1.0.2/SV-L50-Q1.tsv_deletion.q1.fa ")
    print("\tpython ./LengthFilter.py -i ~/NA12878-novaseq/v1.0.2/SV-L50-Q1.tsv_deletion.q1.fa ")
    print("\tpython ./LengthFilter.py -i ~/NA12878-novaseq/v1.0.2/SV-L50-Q1.tsv_softclipping.q1.fa ")
    print("\tpython ./LengthFilter.py -i ~/NA12878-novaseq/v1.0.2/U0.fa ")
    print("\tpython ./LengthFilter.py -i ../data/NA12878/SV-L50-Q1.tsv_insertion.q1.fa ")
    print("\tpython ./LengthFilter.py -i ../data/NA12878/SV-L50-Q1.tsv_deletion.q1.fa ")
    print("\tpython ./LengthFilter.py -i ../data/NA12878/SV-L50-Q1.tsv_softclipping.q1.fa ")
    print("\tpython ./LengthFilter.py -i ../data/NA12878/U0.fa ")

    return


def filter_fun(ifile):
    (num_short, num_medium, num_long) = (0, 0, 0)

    ofd1 = open(ifile+"-s.fa", "w")
    ofd2 = open(ifile+"-m.fa", "w")
    ofd3 = open(ifile+"-l.fa", "w")
    ifd = open(ifile, "r")
    for item in ifd:
        if item.startswith(">"):
            name = item
        else:
            length = len(item.strip())
            if length < SHORT_THRESHOLD:
                ofd1.write("%s%s" % (name, item))
                num_short += 1
            elif length < LONG_THRESHOLD:
                ofd2.write("%s%s" % (name, item))
                num_medium += 1
            else:
                ofd3.write("%s%s" % (name, item))
                num_long += 1

    print("Total\tShort\tMedium\tLong")
    print("%d\t%d\t%d\t%d" % (num_short + num_medium + num_long, num_short, num_medium, num_long))
    ifd.close()
    ofd3.close()
    ofd2.close()
    ofd1.close()
    return


def main(argv):
    ifile = ""

    try:
        opts, args = getopt.getopt(argv, "hi:")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt == "-i":
            ifile = arg

    # error handling for input parameters
    if ifile == "" :
        print("Error: '-i' is required")
        usage()
        sys.exit(2)
    elif not os.path.isfile(ifile):
        print("Error: input file(%s) is not existed" % ifile)
        usage()
        sys.exit(3)

    # Main Function
    filter_fun(ifile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
