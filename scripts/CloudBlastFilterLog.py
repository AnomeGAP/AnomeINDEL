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
# @file    CloudBlastFilterLog.py
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
STR_HUMAN = "Human"
STR_BAC = "Homo sapiens"
STR_NUI = "non-reference unique insertion sequence"
TYPE_UNKNOW = 0
TYPE_HUMAN = 1
TYPE_NUI = 2

# Default Parameter
COVERAGE_THRESH = 1.0
# COVERAGE_THRESH = 0.8
# COVERAGE_THRESH = 0.6
# COVERAGE_THRESH = 0.4
# COVERAGE_THRESH = 0.2
BITSCORE = 200


def usage():
    print("CloudBlastFilterLog.py -i <Input BLAST log file> [-f <FASTA file>] -o <Output BLAST result>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: Input BLAST file")
    print("\t-f: Input FASTA file for reference check (Optional)")
    print("\t-o: Output BLAST result")
    print("Usage:")
    print("\tpython ./CloudBlastFilterLog.py -i ~/NA12878-novaseq/v1.0.2/U1000.alt.fa.blast"
          " -f ~/NA12878-novaseq/v1.0.2/U1000.fa"
          " -o ~/NA12878-novaseq/v1.0.2/U1000.alt.fa.blast.filtered")
    print("\tpython ./CloudBlastFilterLog.py -i ~/NA12878-novaseq/v1.0.2/SV-L50-Q1.tsv_insertion.q1.fa-l.fa.blast"
          " -o ~/NA12878-novaseq/v1.0.2/SV-L50-Q1.tsv_insertion.q1.fa-l.fa.blast.filtered")
    print("\tpython ./CloudBlastFilterLog.py -i ~/NA12878-novaseq/v1.0.2/SV-L50-Q1.tsv_softclipping.q1.fa-l.fa1.blast"
          " -o ~/NA12878-novaseq/v1.0.2/SV-L50-Q1.tsv_softclipping.q1.fa-l.fa1.blast.filtered")

    return


def filter_fun(ifile, ffile, ofile):
    (num_nui, num_human, num_unknown, num_missing) = (0, 0, 0, 0)

    ref = []
    if ffile != "":
        ffd = open(ffile, "r")
        for item in ffd:
            if item.startswith(">"):
                ref.append(item.strip()[1:])
        ffd.close()
        print("%s has %d contigs" % (ffile, len(ref)))
    ofd = open(ofile, "w")
    ifd = open(ifile, "r")
    contig = ""
    name = ""
    idx = 0
    t = TYPE_UNKNOW
    ofd.write("QUERY\tLength\tBIT-Score\tCONTIG\tLength\t#Matches\n")
    for item in ifd:
        items = item.strip().split("\t")
        #print(items)
        if items[0] == "QUERY":
            continue

        if name == "":
            name = items[0]
        elif name != items[0]:
            # print(items)
            if t == TYPE_NUI:
                num_nui += 1
            elif t == TYPE_HUMAN:
                num_human += 1
            else:
                num_unknown += 1
            if contig != "":
                ofd.write("%s" % contig)
            else:
                ofd.write("%s\t*\t0\t*\t*\t0\n" % name)
                print("%s\t Unmatched" % name)
            t = TYPE_UNKNOW
            contig = ""
            while idx < len(ref) and ref[idx] != items[0]:
                if ref[idx] != name:
                    print("%s\tNo Record" % ref[idx])
                    num_missing += 1
                idx += 1
            name = items[0]

        length = int(items[1])
        matches = int(items[5])
        if float(matches) >= float(length) * COVERAGE_THRESH:
            if items[3].find(STR_NUI) > -1 and t < TYPE_NUI:
                t = TYPE_NUI
                contig = item
            elif (items[3].find(STR_HUMAN) > -1 or items[3].find(STR_BAC) > -1) and t < TYPE_HUMAN:
                t = TYPE_HUMAN
                contig = item

    if t == TYPE_NUI:
        num_nui += 1
    elif t == TYPE_HUMAN:
        num_human += 1
    else:
        num_unknown += 1
    if contig != "":
        ofd.write("%s" % contig)
    else:
        ofd.write("%s\t*\t0\t*\t*\t0\n" % name)
        print("%s\t Unmatched" % name)
    while idx < len(ref):
        if ref[idx] != name:
            print("%s\tNo Record" % ref[idx])
            num_missing += 1
        idx += 1

    if ffile != "":
        print("Total\tNUI\tHUMAN\tUnknown\tMissing")
        print("%d\t%d\t%d\t%d\t%d" % (num_nui + num_human + num_unknown + num_missing, num_nui, num_human, num_unknown,
                                      num_missing))
    else:
        print("Total\tNUI\tHUMAN\tUnknown")
        print("%d\t%d\t%d\t%d" % (num_nui + num_human + num_unknown, num_nui, num_human, num_unknown))
    ifd.close()
    ofd.close()
    return


def main(argv):
    ifile = ""
    ffile = ""
    ofile = ""

    try:
        opts, args = getopt.getopt(argv, "hi:f:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt == "-i":
            ifile = arg
        elif opt == "-f":
            ffile = arg
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
    filter_fun(ifile, ffile, ofile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
