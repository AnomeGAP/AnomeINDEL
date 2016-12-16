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
# @file    PairReadGenerator.py
#
# @brief   An generator for pair-end reads
#
# @author  A-Tsai Su(chungtsai_su@anome.com.tw)
#
# @date    2016/12/16
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

# Parameter setting


READ_LENGTH = 101
SEGMENT_LENGTH = 400
DEVIATION_SEGMENT_LENGTH = 100
COVERAGE = 30
SEQUENCE_ERROR_RATE=0.01


def _usage():
    """ Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    """
    print(
        "PairReadGenerator.py -i <FASTA file> -o <Output file> -l <read length> -s <segment length> "
        "-d <deviation of segment length> -c <coverage> -e <rate of sequence error>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <FASTA file>")
    print("\t-o: <Output file>")
    print("\t-l: <read length>")
    print("\t-s: <segment length>")
    print("\t-d: <deviation of segment length>")
    print("\t-c: <coverage>")
    print("\t-e: <rate of sequence error>")

    print("Usage:")
    print("\tpython3 ~/src/github/AnomeINDEL/scripts/PairReadGenerator.py -i ~/data/hg19/chr22.fa -o ~/data/hg19/chr22 "
          "-l 101 -s 400 -d 100 -c 30 -e 0.01")

    return


def mutate(ori):
    bases = ('A', 'C', 'G', 'T')
    i = random.randrange(4)
    while bases[i] == ori:
        i = random.randrange(4)

    return bases[i]


def reverse_complement(read):
    bases = ('A', 'C', 'G', 'T')
    slist = list(reversed(read))
    for ndx, base in enumerate(slist):
        if slist[ndx] == 'A':
            slist[ndx] = 'T'
        elif slist[ndx] == 'T':
            slist[ndx] = 'A'
        elif slist[ndx] == 'C':
            slist[ndx] = 'G'
        elif slist[ndx] == 'G':
            slist[ndx] = 'C'

    return "".join(slist)


def generator_reads(o1fd, o2fd, efd, idx, start, ref):

    for i in range(int(COVERAGE/2)):
        # random
        s1 = random.randrange(DEVIATION_SEGMENT_LENGTH)
        s2 = s1 + SEGMENT_LENGTH - READ_LENGTH + random.randrange(-DEVIATION_SEGMENT_LENGTH, DEVIATION_SEGMENT_LENGTH)
        r1 = ref[s1:s1+READ_LENGTH]
        # Force read2 to reverse complement
        r2 = reverse_complement(ref[s2:s2+READ_LENGTH])

        # add noise
        l1 = list(r1)
        l2 = list(r2)
        efd.write("#ID\t#pair\tPOS\tREF\tALT\tREF_POS\n")
        for j in range(READ_LENGTH):
            if random.randrange(int(1/SEQUENCE_ERROR_RATE)) == 0:
                rep = mutate(l1[j])
                efd.write("%d\t1\t%d\t%s\t%s\t%d\n" % (idx, j, l1[j], rep, start+s1+j))
                l1[j] = rep
        for j in range(READ_LENGTH):
            if random.randrange(int(1 / SEQUENCE_ERROR_RATE)) == 0:
                rep = mutate(l2[j])
                efd.write("%d\t2\t%d\t%s\t%s\t%d\n" % (idx, j, l2[j], rep, start+s2+READ_LENGTH-j))
                l2[j] = rep

        r1 = "".join(l1)
        o1fd.write(">READ%d:1\t%d\n%s\n+\n%s\n" % (idx, start+s1, r1, "="*READ_LENGTH))
        r2 = "".join(l2)
        o2fd.write(">READ%d:2\t%d\n%s\n+\n%s\n" % (idx, start+s2+READ_LENGTH, r2, "="*READ_LENGTH))
        idx += 1
        if idx % 1000 == 0:
            print("%d\t%d" % (idx, start))

    return idx


def pair_read_generator(ifn, ofn):
    o1fd = open(ofn + "_1.fa", "w")
    o2fd = open(ofn + "_2.fa", "w")
    efd = open(ofn + "_error.tsv", "w")
    ifd = open(ifn, "r")
    ref = ""
    start = 0
    idx = 1

    for line in ifd:
        if re.match(">", line):
            continue
        ref += line.strip().upper()
#        print("%6d %s" % (start, ref))
        # remove 'N'
        i = 0
        while i < len(ref) and ref[i] == "N":
            i += 1
        if i > 0:
            start += i
            ref = ref[i:]
        if len(ref) < SEGMENT_LENGTH+2*DEVIATION_SEGMENT_LENGTH:
            continue
        # skip single "N" region
        i = ref.rfind("N")
        if i > 0:
            start += i
            ref = ref[i:]
            continue
#        print("%6d %s" % (start, ref))
        idx = generator_reads(o1fd, o2fd, efd, idx, start, ref)
        ref = ref[int((SEGMENT_LENGTH+2*DEVIATION_SEGMENT_LENGTH)/2):]
        start += int((SEGMENT_LENGTH+2*DEVIATION_SEGMENT_LENGTH)/2)

    ifd.close()
    efd.close()
    o2fd.close()
    o1fd.close()

    return


def main(argv):
    infile = ""
    outfile = ""
    global READ_LENGTH, SEGMENT_LENGTH, DEVIATION_SEGMENT_LENGTH, COVERAGE, SEQUENCE_ERROR_RATE

    try:
        opts, args = getopt.getopt(argv, "hi:o:l:s:d:c:e:")
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
        elif opt in "-l":
            READ_LENGTH = int(arg)
        elif opt in "-s":
            SEGMENT_LENGTH = int(arg)
        elif opt in "-d":
            DEVIATION_SEGMENT_LENGTH = int(arg)
        elif opt in "-c":
            COVERAGE = int(arg)
        elif opt in "-e":
            SEQUENCE_ERROR_RATE = float(arg)
    if not infile or not os.path.isfile(infile):
        print("Error: input file(%s) is not existed" % (infile))
        _usage()
        sys.exit(2)

    # Main Function
    pair_read_generator(infile, outfile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
