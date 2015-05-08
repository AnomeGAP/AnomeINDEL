#!/bin/env python
#
# @note Copyright (C) 2015, Anome Incorporated. All Rights Reserved.
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
# @file    SegmentedSuffix.py
#
# @brief   Given a DNA sequence (FASTA file), segmenting them by a user defined length and 
#          then generate all of suffix for each segment.
#
# @author  Chung-Tsai Su(chungtsai_su@anome.com.tw)
#
# @date    2015/05/15
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os.path
import math
import re
from collections import defaultdict

SEGMENT_LENGTH=200
MINIMAL_LENGTH=10

def Usage():
    print("SegmentedSuffix.py -i <Input FASTA file> -l <segment length> -d <minimal length of suffix> -o <Output CSV file>") 
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: Input FASTA file")
    print("\t-l: Segment length [Default:%d]" % (SEGMENT_LENGTH))
    print("\t-d: Minimal length of suffix [Default:%d]" % (MINIMAL_LENGTH))
    print("\t-o: Output CSV file")
    print("Usage:")
    print("\t./SegmentedSuffix.py -i ../example/chr1.fa -l 200 -d 10 -o ../example/chr1.suffix.csv")

    return

def SuffixExtend(items, segment, min_len):
    length = len(segment)-min_len+1
    for i in range(0,length):
        item = segment[i:]
        items.append(item)
        print item
    return

def SegmentedSuffix(ifn, slength, min_len, ofn):
    ifd = open(ifn, "r")
    total_length = 0
    idx = 0
    seq = ""
    items = []

    for line in ifd:
        if re.match("^>", line):
           print "Skip comment: %s" % (line)
           continue

        line = line.strip()
        print "line : %s" % (line)
        seq += line.upper()
        length  = len(seq)
        while length > slength:
            segment = seq[:slength]
            SuffixExtend(items, segment, min_len) 
            seq = seq[slength:]
            length = len(seq)

    items.append(seq)
    print "add %s" % (seq)
    ifd.close()
    items.sort()
    
    ofd = open(ofn, "w")
    for item in items:
        ofd.write("%s\n" % item)

    ofd.close()
    return

def main(argv):
    ifile = ""
    slength = SEGMENT_LENGTH
    min_len = MINIMAL_LENGTH
    ofile = ""

    try:
        opts, args = getopt.getopt(argv,"hi:l:d:o:")
    except getopt.GetoptError:
        Usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            Usage()
            sys.exit()
        elif opt in ("-i"):
            ifile = arg
        elif opt in ("-d"):
            min_len = int(arg)
        elif opt in ("-o"):
            ofile = arg

    #error handling for input parameters
    if (ifile == ""):
        print("Error: '-i' is required")
        Usage()
        sys.exit(2)
    elif (os.path.isfile(ifile) == False): 
        print("Error: input file(%s) is not existed" % (ifile))
        Usage()
        sys.exit(3)

    if (ofile == ""):
        ofile = "%s.suffix.csv" % (ifile)

    #Main Function
    SegmentedSuffix(ifile, slength, min_len, ofile)

    return

if __name__ == '__main__':
    main(sys.argv[1:])

