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
# @file    BAManalyzer.py
#
# @brief   Parsing BAM file to analyze the distribution of soft-clipping reads for 10X Genomics data
#
# @author  Chung-Tsai Su(chungtsai_su@atgenomix.com)
#
# @date    2018/11/23
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os.path
from collections import defaultdict

# CONST

# Default Parameter


def usage():
    print("VCF2BED.py -i <Input VCF> -o <Output BED>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: Input VCF file")
    print("\t-o: Output BED")
    print("Usage:")
    print("\tpython ./VCF2BED.py -i ./NA12878.sorted.vcf -o ./NA12878.sorted.bed ")

    return


def analyzer(ifn, ofn):
    cnt = 0
    ofd = open(ofn, "w")
    ifd = open(ifn, "r")
    for line in ifd:
        if line.startswith("#"):
            continue
        items = line.strip().split('\t')
        infos = items[7].split(";")
        for i in infos:
            if i.startswith("END"):
                # print("%s\t%s\t%s\t%s\t%s" % (items[0], items[1], i[4:], items[4], items[6]))
                ofd.write("%s\t%s\t%s\t%s\t%s\n" % (items[0], items[1], i[4:], items[4], items[6]))
                break
        cnt += 1
    print("Total: %d records" % cnt)
    ofd.close()
    ifd.close()
    return


def main(argv):
    ifile = ""
    ofile = ""

    try:
        opts, args = getopt.getopt(argv, "hi:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt == "-i":
            ifile = arg
        elif opt == "-t":
            type_setting = arg

    # error handling for input parameters
    if ifile == "":
        print("Error: '-i' is required")
        usage()
        sys.exit(2)
    elif not os.path.isfile(ifile):
        print("Error: input file(%s) is not existed" % ifile)
        usage()
        sys.exit(3)

    if ofile == "":
        ofile = "%s.bed" % ifile

    # Main Function
    analyzer(ifile, ofile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
