#!/bin/env python3
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
# @file    DeletionSequence.py
#
# @brief   Attach sequence based on FASTA header (for long deletions generated from BAMSVStats.py
#
# @author  A-Tsai Su(chungtsai.su@atgenomix.com)
#
# @date    2019/06/11
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

READ_LENGTH = 151
SEGMENT_LENGTH = 70
DEVIATION_SEGMENT_LENGTH = 150
COVERAGE = 30
SEQUENCE_ERROR_RATE = 0.001


def _usage():
    """ Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    """
    print("DeletionSequence.py -i <FASTA file> -o <Output file> -r <Reference file>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <FASTA file>")
    print("\t-o: <Output file>")
    print("\t-r: <Reference file>")
    print("Usage:")
    print("\tpython3 ~/src/github/AnomeINDEL/scripts/DeletionSequence.py "
          "-i ~/NA12878-novaseq/v1.0.2/SV.tsv_deletion.q1.header "
          "-r ~/data/hg38/ref.fa "
          "-o ~/NA12878-novaseq/v1.0.2/SV.tsv_deletion.q1.fa")
    print("\tpython3 ./DeletionSequence.py "
          "-i ../data/NA12878//SV-L50-Q1.tsv_deletion.q1.header "
          "-r ../data/hg38/ref.fa "
          "-o ../data/NA12878//SV-L50-Q1.tsv_deletion.q1.fa")

    return


def offset_by_chromosome(idx):
    if idx == 'chrX':
        return 22 * 1000000000
    elif idx == 'chrY':
        return 23 * 1000000000
    elif idx == 'chrM':
        return 24 * 1000000000
    else:
        return (int(idx[3:])-1) * 1000000000


def deletion_sequence_generator(ifn, rfn, ofn):
    l_start = []
    h_header =  defaultdict(str)
    h_len = defaultdict(int)

    ifd = open(ifn, "r")
    for line in ifd:
        if re.match(">", line):
            items = line.strip().split(" ")
            (chrom, start) = items[2].split(":")
            pos = offset_by_chromosome(chrom) + int(start)
            # print("[%s-%s] => %d" % (items[2], items[1], pos))
            l_start.append(pos)
            h_header[pos] = line.strip()
            h_len[pos] = int(items[1])
    ifd.close()

    print("There are %d long deletions" % len(l_start))

    ref = ""
    start = 0
    idx = 0
    offset = 0
    h_seq = defaultdict(str)
    rfd = open(rfn, "r")
    ofd = open(ofn, "w")
    for line in rfd:
        if re.match(">", line):
            if len(h_seq) > 0:
                print("[ERROR] h_len should be empty!! %s" % line)
                return
            if line[1:].split(" ")[0] == "chrM" or len(line[1:].split(" ")[0]) > 5:
                # print(line[1:].split(" ")[0])
                ref = ""
                offset = -1
            else:
                offset = offset_by_chromosome(line[1:].split(" ")[0])
        elif offset != -1:
            ref += line.strip().upper()
            while idx < len(l_start) and offset + len(ref) > l_start[idx]:
                # print("add %d : %d %d" % (idx, l_start[idx], h_len[l_start[idx]]))
                h_seq[l_start[idx]] = ""
                idx += 1
            min_key = 25000000000
            l_del = []
            for key in h_seq.keys():
                if min_key > key:
                    min_key = key
                if key + h_len[key] <= offset + len(ref):
                    ofd.write("%s\n%s\n" % (h_header[key], ref[key-offset:key-offset+h_len[key]]))
                    # print("save %d " % key)
                    l_del.append(key)
            for key in l_del:
                h_seq.pop(key)
                # print("del %d when ref=%d" % (key, offset + len(ref)))
            if min_key != 25000000000 and min_key > offset:
                # print("mov ref from %d to %d" % (offset, min_key))
                ref = ref[min_key-offset:]
                offset = min_key

    rfd.close()
    ofd.close()

    return


def main(argv):
    infile = ""
    outfile = ""
    reffile = ""

    try:
        opts, args = getopt.getopt(argv, "hi:o:r:")
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
        elif opt in "-o":
            outfile = arg
        elif opt in "-r":
            reffile = arg

    if not infile or not os.path.isfile(infile):
        print("Error: input file(%s) is not existed" % infile)
        _usage()
        sys.exit(2)
    if not reffile or not os.path.isfile(reffile):
        print("Error: reference file(%s) is not existed" % reffile)
        _usage()
        sys.exit(3)

    # Main Function
    deletion_sequence_generator(infile, reffile, outfile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
