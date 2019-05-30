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
import re
import pysam  # http://pysam.readthedocs.org/en/latest/api.html
from collections import defaultdict

# CONST
TYPE_SEQ = "SEQ"
TYPE_CIGAR = "CIGAR"

# Default Parameter
READ1_LEN = 128
READ2_LEN = 151
READ1_CIGAR = "%dM" % READ1_LEN
READ2_CIGAR = "%dM" % READ2_LEN
MAX_INSERT_SIZE = 1000


def usage():
    print("BAQuery.py -i <Input BAM> -t <Type> -q <keyword> -o <Output SAM>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: Input BAM file")
    print("\t-t: [%s|%s] (Default: %s)" % (TYPE_SEQ, TYPE_CIGAR, TYPE_SEQ))
    print("\t-q: keyword")
    print("\t-o: Output SAM file which matches the query")
    print("Usage:")
    print("\tpython ./BAQuery.py -i ./NA12878.sorted.bam "
          "-t %s -q AAGAGCACACGTCTGAACTCCAGTCACTGCTGTAAATCTCGTATGCCGTCTTCTGCTTGAAAAA -o ./NA12878.query.sorted.sam "
          "> output.log" % TYPE_SEQ)

    return


def analyzer(ifn, ofn, query):
    samfile = pysam.AlignmentFile(ifn, "rb")
    output = pysam.AlignmentFile(ofn, "wh", template=samfile)
    cnt = 0
    print("Number of mapped reads = %d" % samfile.mapped)
    print("Number of unmapped reads = %d" % samfile.unmapped)
    print("Unmapped ratio = %.2f%%" % (float(samfile.unmapped * 100) / (samfile.mapped + samfile.unmapped)))
    print("nocoordinate=%d" % samfile.nocoordinate)
    print("")

    for read in samfile.fetch():
        if str(read.query_sequence).find(query) == 0:
            output.write(read)
            print("%s\t%d\t%d\t%d\t%s\t%s\t%s" %(read.reference_name, read.reference_start, read.is_read1, read.is_proper_pair,
                                             read.cigarstring, read.query_name, read.query_sequence))
            cnt += 1

    print("Total: %d reads" % cnt)
    output.close()
    samfile.close()
    return


def main(argv):
    ifile = ""
    ofile = ""
    type_setting = TYPE_SEQ
    query = ""

    try:
        opts, args = getopt.getopt(argv, "hi:t:q:o:")
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
            type = arg
        elif opt == "-q":
            query = arg
        elif opt == "-o":
            ofile = arg

    # error handling for input parameters
    if ifile == "":
        print("Error: '-i' is required")
        usage()
        sys.exit(2)
    elif not os.path.isfile(ifile):
        print("Error: input file(%s) is not existed" % ifile)
        usage()
        sys.exit(3)
    elif query == "":
        print("Error: '-q' is required")
        usage()
        sys.exit(4)

    if ofile == "":
        ofile = "%s.unmapped.sam" % ifile

    # Main Function
    if type_setting == TYPE_SEQ:
        analyzer(ifile, ofile, query)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
