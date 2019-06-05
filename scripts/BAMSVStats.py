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
# @file    BAMSVStats.py
#
# @brief   Parsing BAM file to extract large SV for Connected-Reads and estmiate the amount of them
#
# @author  Chung-Tsai Su(chungtsai.su@atgenomix.com)
#
# @date    2019/06/03
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

# Default Parameter
MIN_LENGTH = 50
DIVISOR = 50
NUM = 100


def usage():
    print("BAMSVStats.py -i <Input BAM> -l <Minimal Length> -o <Output tsv list>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: Input BAM file")
    print("\t-l: minimal length (default: %d )" % MIN_LENGTH)
    print("\t-o: Output tsv list")
    print("Usage:")
    print("\tpython ./BAMSVStats.py -i ~/NA12878-novaseq/v1.0.2/result-1.0.2.primary.sorted.bam -l %d "
          "-o ~/NA12878-novaseq/v1.0.2/SV.tsv" % MIN_LENGTH)

    return

#     | POS | SEQ
# ----+-----+-----
#   M |  +  |  +
#   I |     |  +
#   D |  +  |
#   S |     |  +
#   H |     |


def analyzer(ifn, min_length, ofn):
    samfile = pysam.AlignmentFile(ifn, "rb")
    output = open(ofn, "w")

    total_contigs = 0
    num_insertions = 0
    num_deletions = 0
    num_others = 0
    num_cnv = 0
    i_count = [0] * NUM
    d_count = [0] * NUM
    hs_count = [0] * NUM
    h_hits = defaultdict(int)

    max_hit = 0
    max_pos = ""

    # mapped reads for Insertion
    for read in samfile.fetch():
        if read.reference_name == 'chrM':
            break

        # if any(re.findall(r'I|D|H|S', str(read.cigarstring), re.IGNORECASE)):
        if any(re.findall(r'I', str(read.cigarstring), re.IGNORECASE)) and read.mapping_quality > 0:
            # print("cigar=%s" % read.cigarstring)
            items = re.split("([0-9]+\S)", str(read.cigarstring))
            # print(items)
            idx = 0
            is_matched = 0

            for i in range(len(items)):
                # print("%d:[%s]" % (i, items[i]))
                if items[i] == "":
                    continue
                # print(items[i][:-1])
                l = int(items[i][:-1])
                if "I" in items[i] and l >= min_length:
                    is_matched += 1
                    if max_hit < l:
                        max_hit = l
                        max_pos = "%s\t%s:%d" % (read.query_name, read.reference_name, read.reference_start)
                    key = "%s:%d" % (read.reference_id, int((read.reference_start+idx) / DIVISOR))
                    if h_hits[key] == 0:
                        output.write("I\t%s\t%s\t%d\n" % (read.query_name, read.reference_name,
                                                          int(read.reference_start) + idx))
                        if l > DIVISOR * NUM:
                            i_count[NUM-1] += 1
                        else:
                            i_count[int(l/DIVISOR)] += 1
                        h_hits[key] += 1

                if "M" in items[i] or "D" in items[i]:
                        idx += l

            if is_matched > 0:
                num_insertions += 1
        total_contigs += 1

    print("There are %d contigs" % total_contigs)
    print("There are %d long insertions" % num_insertions)
    print("There are %d distinct long insertions" % len(h_hits))
    print("Maximal length of hit : %d at %s" % (max_hit, max_pos))

    max_hit = 0
    # mapped reads for Deletion
    for read in samfile.fetch():
        if read.reference_name == 'chrM':
            break

        if any(re.findall(r'D', str(read.cigarstring), re.IGNORECASE)) and read.mapping_quality > 0:
            items = re.split("([0-9]+\S)", str(read.cigarstring))
            idx = 0
            is_matched = 0

            for i in range(len(items)):
                if items[i] == "":
                    continue
                l = int(items[i][:-1])
                if "D" in items[i] and l >= min_length:
                    is_matched += 1
                    if max_hit < l:
                        max_hit = l
                        max_pos = "%s\t%s:%d" % (read.query_name, read.reference_name, read.reference_start)
                    key = "%s:%d" % (read.reference_id, int((read.reference_start+idx) / DIVISOR))
                    if h_hits[key] == 0:
                        output.write("D\t%s\t%s\t%d\n" % (read.query_name, read.reference_name,
                                                          int(read.reference_start) + idx))
                        if l > DIVISOR * NUM:
                            d_count[NUM - 1] += 1
                        else:
                            d_count[int(l / DIVISOR)] += 1
                        h_hits[key] += 1

                if "M" in items[i] or "D" in items[i]:
                    idx += l

            if is_matched > 0:
                num_deletions += 1
        total_contigs += 1

    print("There are %d long deletions" % num_deletions)
    print("There are %d distinct long insertions + long deletions" % len(h_hits))
    print("Maximal length of hit : %d at %s" % (max_hit, max_pos))

    h_clipping = defaultdict(int)

    max_hit = 0
    # mapped reads for Deletion
    for read in samfile.fetch():
        if read.reference_name == 'chrM':
            break

        if any(re.findall(r'H|S', str(read.cigarstring), re.IGNORECASE)) and read.mapping_quality > 0:
            # print("%s" % read.cigarstring)
            items = re.split("([0-9]+\S)", str(read.cigarstring))
            idx = 0
            is_matched = 0

            for i in range(len(items)):
                if items[i] == "":
                    continue
                l = int(items[i][:-1])
                # print("%d\t%d\t%s" % (len(items), i, items[i]))
                #if ("H" in items[i] or "S" in items[i]) and l >= min_length and i == (len(items) - 2):
                if ("H" in items[i] or "S" in items[i]) and l >= min_length:
                    # print("YES\t%d\t%s" % (i, items[i]))
                    is_matched += 1
                    if i == 1:
                        if max_hit < l:
                            max_hit = l
                            max_pos = "%s\t%s:%d" % (read.query_name, read.reference_name, read.reference_start)
                        key = "%s:%d" % (read.reference_id, int((read.reference_start+idx) / DIVISOR))
                        if h_hits[key] == 0:
                            output.write("H|S\t%s\t%s\t%d\n" % (read.query_name, read.reference_name,
                                                                int(read.reference_start) + idx))
                            if l > DIVISOR * NUM:
                                hs_count[NUM - 1] += 1
                            else:
                                hs_count[int(l / DIVISOR)] += 1
                            h_hits[key] += 1

                if "M" in items[i] or "D" in items[i]:
                    idx += l

            if is_matched > 0:
                num_others += 1
                if is_matched > 1:
                    num_cnv += 1
                    output.write("CNV\t%s\t%s\t%d\n" % (read.query_name, read.reference_name, read.reference_start))

        total_contigs += 1

    print("There are %d long num_others" % num_others)
    print("There are %d potential cnv" % num_cnv)
    print("There are %d distinct long insertions + long deletions + long others" % len(h_hits))
    print("Maximal length of hit : %d at %s" % (max_hit, max_pos))

    for i in range(NUM):
        print("%d\t%d\t%d\t%d" % (i, i_count[i], d_count[i], hs_count[i]))
    output.close()
    samfile.close()
    return


def main(argv):
    ifile = ""
    ofile = ""
    min_length = MIN_LENGTH

    try:
        opts, args = getopt.getopt(argv, "hi:l:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt == "-i":
            ifile = arg
        elif opt == "-l":
            min_length = int(arg)
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
    if ofile == "":
        ofile = "%s.tsv" % ifile

    # Main Function
    analyzer(ifile, min_length, ofile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
