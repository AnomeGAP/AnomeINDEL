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
# @file    VCFCloseHeterozygous.py
#
# @brief   An example to filter the close variants
#
# @author  A-Tsai Su(chungtsai_su@anome.com.tw)
#
# @date    2018/06/13
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os.path
import re

MIN_GAP = 100


def usage():
    """
        Description: Program Usage
        Argument:    NONE
        Return:	     NONE """
    print("VCFCloseHeterozygous.py -i <VCF file> -g <integer>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <VCF file>")
    print("\t-g: <gap length")
    print("Usage:")
    print("\t./VCFCloseHeterozygous.py -i twsz009b3-heterozygous.vcf -g 100")

    return


def vcfCloseVariant(ifn, max_gap):
    ifd = open(ifn, "r")
    # Read sequence and segment them with length of slength
    prev = 0
    buf = ""
    for line in ifd:
        if re.match("^#", line):
            # print "Skip comment: %s" % (line)
            continue

        line = line.strip()
        items = line.split('\t')
        # print(items[1]) #POS
        pos = int(items[1])
        if abs(pos - prev) <= max_gap:
            print("%s\n%s" % (buf, line))

        prev = pos
        buf = line
    ifd.close()

    return


def main(argv):
    infile = ""
    gap = MIN_GAP

    try:
        opts, args = getopt.getopt(argv, "hi:g:")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in "-i":
            infile = arg
        elif opt in "-g":
            gap = int(arg)
    if not infile:
        usage()
        sys.exit(2)
    elif not os.path.isfile(infile):
        print("Error: input file(%s) is not existed" % (infile))
        Usage()
        sys.exit(3)

    # Main Function
    vcfCloseVariant(infile, gap)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
