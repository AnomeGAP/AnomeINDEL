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
# @file    CloudBlastFilter.py
#
# @brief   Blast Result Filtering\
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
STR_HUMAN = "Homo sapiens"
STR_BAC = "Homo sapiens BAC"
STR_NUI = "non-reference unique insertion sequence"
TYPE_UNKNOW = 0
TYPE_HUMAN = 1
TYPE_BAC = 2
TYPE_NUI = 3

# Default Parameter
COVERAGE_THRESH = 0.8
BITSCORE = 200

def usage():
    print("CloudBlastFilter.py -i <Input BLAST log file> -o <Output BLAST result>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: Input BLAST file")
    print("\t-o: Output BLAST result")
    print("Usage:")
    print("\tpython ./CloudBlastFilter.py -i ~/NA12878-novaseq/v1.0.2/U1000.alt.fa.log "
          "-o ~/NA12878-novaseq/v1.0.2/U1000.alt.fa.log.filtered")

    return


def filter_fun(ifile, ofile):
    (num_nui, num_bac, num_human, num_unknown) = (0, 0, 0, 0)
    t = TYPE_UNKNOW

    ofd = open(ofile, "w")
    ifd = open(ifile, "r")
    contig = ""
    ofd.write("QUERY\tLength\tBIT-Score\tCONTIG\tLength\t#Matches\n")
    for item in ifd:
        items = item.strip().split("\t")
        if len(items) == 2:
            # print(items)
            if t == TYPE_NUI:
                num_nui += 1
            elif t == TYPE_BAC:
                num_bac += 1
            elif t == TYPE_HUMAN:
                num_human += 1
            else:
                num_unknown += 1
            if contig != "":
                ofd.write("%s" % contig)
            else:
                ofd.write("%s\t*\t0\t*\t*\t0\n" % (items[0]))
            t = TYPE_UNKNOW
            contig = ""
            continue

        length = int(items[1])
        matches = int(items[5])
        if float(matches) >= float(length) * COVERAGE_THRESH:
            if items[3].find(STR_NUI):
                t = TYPE_NUI
                contig = item
            elif items[3].find(STR_BAC) and t < TYPE_BAC:
                t = TYPE_BAC
                contig = item
            elif items[3].find(STR_HUMAN) and t < TYPE_HUMAN:
                t = TYPE_HUMAN
                contig = item

    print("Total\tNUI\tBAC\tHUMAN\tUnknown")
    print("%d\t%d\t%d\t%d\t%d" % (num_nui+num_bac+num_unknown, num_nui, num_bac, num_human, num_unknown))
    ifd.close()
    ofd.close()
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
        elif opt == "-o":
            ofile = arg

    if ofile == "":
        ofile = "%s.blast.log" % ifile

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
    filter_fun(ifile, ofile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
