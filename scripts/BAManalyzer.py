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

# Default Parameter
READ1_LEN = 128
READ2_LEN = 151
READ1_CIGAR = "%dM" % READ1_LEN
READ2_CIGAR = "%dM" % READ2_LEN
MAX_INSERT_SIZE = 1000


def usage():
    print("BAManalyzer.py -i <Input BAM> -o <Output SAM>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: Input BAM file")
    print("\t-o: Output SAM file which contains soft-clipping reads")
    print("Usage:")
    print("\tpython ./BAManalyzer.py -i ./NA12878.sorted.bam -o ./NA12878.sc.sorted.sam > output.log")

    return


def analyzer(ifn, ofn):
    samfile = pysam.AlignmentFile(ifn, "rb")
    output = pysam.AlignmentFile(ofn, "wh", template=samfile)

    print("Number of mapped reads = %d" % samfile.mapped)
    print("Number of unmapped reads = %d" % samfile.unmapped)
    print("Unmapped ratio = %.2f%%" % (float(samfile.unmapped * 100) / (samfile.mapped + samfile.unmapped)))
    print("nocoordinate=%d" % samfile.nocoordinate)
    print("")
    num_pairs = 0
    num_qualified = 0
    insertsize = [0] * MAX_INSERT_SIZE
    insertsize4sc = [0] * MAX_INSERT_SIZE
    qual1 = [0] * 100
    qual2 = [0] * 100
    qualALL1 = [0] * 100
    qualALL2 = [0] * 100
    (num_softclip1,num_softclip2) = (0, 0)
    (num_sc1head, num_sc1tail, num_sc2head, num_sc2tail) = (0, 0, 0, 0)
    sc1 = [0] * 200
    sc2 = [0] * 200
    h_1fh = defaultdict(int)
    h_1ft = defaultdict(int)
    h_1rh = defaultdict(int)
    h_1rt = defaultdict(int)
    h_2fh = defaultdict(int)
    h_2ft = defaultdict(int)
    h_2rh = defaultdict(int)
    h_2rt = defaultdict(int)

    for read in samfile.fetch():
        if read.is_read1:
            if read.is_proper_pair and read.cigarstring == READ1_CIGAR and read.get_tag("MC") == READ2_CIGAR:
                # print("query_name = %s" % read.query_name)
                # print("cigarstring = %s" % read.cigarstring)
                # print("get_tag(MC) = %s" % read.get_tag("MC"))
                # print("reference_start = %d" % read.reference_start)
                # print("next_reference_start = %d" % read.next_reference_start)
                # print("query_alignment_length = %d" % read.query_alignment_length)
                if read.reference_start <= read.next_reference_start:
                    diff = read.next_reference_start + READ2_LEN - read.reference_start
                else:
                    diff = read.reference_start + READ1_LEN - read.next_reference_start
                if diff >= MAX_INSERT_SIZE:
                    diff = MAX_INSERT_SIZE - 1
                if diff < READ2_LEN:
                    diff = READ2_LEN
                insertsize[diff] += 1
                qual1[read.mapping_quality] += 1
                # print("insert size=%d" % diff )

                # output.write(read)
                num_qualified += 1
                # exit()
            elif read.is_proper_pair :
                if str(read.cigarstring).find('S') >= 0:
                    num_softclip1 += 1
                    items = re.split("([0-9]+S)", str(read.cigarstring))
                    # print("cigar=%s" % read.cigarstring)
                    output.write(read)
                    for i in range(len(items)):
                        # print("%d:[%s]" % (i, items[i]))
                        if 'S' in items[i]:
                            if i <= 1:
                                num_sc1head += 1
                                n = int(items[i].split("S")[0])
                                if n >= 20:
                                    if read.is_reverse:
                                        h_1rh["1RH-%s" % read.query_sequence[:n]] += 1
                                    else:
                                        h_1fh["1FH-%s" % read.query_sequence[:n]] += 1
                                    #print("1H\t%d\t%s" % (read.is_reverse, read.query_sequence[:n]))
                            elif i >= len(items) - 2:
                                num_sc1tail += 1
                                n = int(items[i].split("S")[0])
                                if n >= 20:
                                    if read.is_reverse:
                                        h_1rt["1RT-%s" % read.query_sequence[:n]] += 1
                                    else:
                                        h_1ft["1FT-%s" % read.query_sequence[:n]] += 1
                                    #print("1T\t%d\t%s" % (read.is_reverse, read.query_sequence[:n]))
                                sc1[int(items[i].split("S")[0])] += 1
                if read.reference_start and read.next_reference_start and read.has_tag("MC") and str(read.cigarstring)\
                        and (str(read.cigarstring).find('S') >= 0 or str(read.get_tag("MC")).find('S') >= 0):
                    items = re.split("([0-9]+S)", str(read.get_tag("MC")))
                    s = 0
                    for i in range(len(items)):
                        if 'S' in items[i]:
                            s += int(items[i].split("S")[0])
                    end2 = read.next_reference_start + 151 - s
                    items = re.split("([0-9]+S)", str(read.cigarstring))
                    s = 0
                    for i in range(len(items)):
                        if 'S' in items[i]:
                            s += int(items[i].split("S")[0])
                    end1 = read.reference_start + READ1_LEN - s
                    if end1 >= end2:
                        if read.reference_start <= read.next_reference_start:
                            diff = end1 - read.reference_start
                        else:
                            diff = end1 - read.next_reference_start
                    else:
                        if read.reference_start <= read.next_reference_start:
                            diff = end2 - read.reference_start
                        else:
                            diff = end2 - read.next_reference_start
                    if diff >= MAX_INSERT_SIZE:
                        diff = MAX_INSERT_SIZE - 1
                    insertsize4sc[diff] += 1
                    # if diff == READ1_LEN:
                    #     print("%s\t%s\t%d\t%d" % (read.cigarstring, read.get_tag("MC"), read.reference_start, read.next_reference_start))
            num_pairs += 1
            qualALL1[read.mapping_quality] += 1
        elif not read.is_read1:
            if read.is_proper_pair and read.cigarstring == READ2_CIGAR and read.get_tag("MC") == READ1_CIGAR:
                qual2[read.mapping_quality] += 1
            elif str(read.cigarstring).find('S') >= 0:
                num_softclip2 += 1
                items = re.split("([0-9]+S)", str(read.cigarstring))
                # print("cigar=%s" % read.cigarstring)
                output.write(read)
                for i in range(len(items)):
                    # print("%d:[%s]" % (i, items[i]))
                    if 'S' in items[i]:
                        if i <= 1:
                            num_sc2head += 1
                            n = int(items[i].split("S")[0])
                            if n >= 20:
                                if read.is_reverse:
                                    h_2rh["2RH-%s" % read.query_sequence[:n]] += 1
                                else:
                                    h_2fh["2FH-%s" % read.query_sequence[:n]] += 1
                                #print("2H\t%d\t%s" % (read.is_reverse, read.query_sequence[:n]))
                        elif i >= len(items) - 2:
                            num_sc2tail += 1
                            n = int(items[i].split("S")[0])
                            if n >= 20:
                                if read.is_reverse:
                                    h_2rt["2RT-%s" % read.query_sequence[:n]] += 1
                                else:
                                    h_2ft["2FT-%s" % read.query_sequence[:n]] += 1
                                #print("2T\t%d\t%s" % (read.is_reverse, read.query_sequence[:n]))
                        sc2[int(items[i].split("S")[0])] += 1

            qualALL2[read.mapping_quality] += 1
        # debug
        # if num_softclip2 > 1000000:
        #     break
    # debug for adaptor
    topN = 100
    i = 0
    for key, value in sorted(h_1fh.items(), key=lambda item: (item[1], item[0]), reverse=True):
        if i >= topN:
            break
        print("%s\t%d" % (key, value))
        i += 1
    i = 0
    for key, value in sorted(h_1rh.items(), key=lambda item: (item[1], item[0]), reverse=True):
        if i >= topN:
            break
        print("%s\t%d" % (key, value))
        i += 1
    i = 0
    for key, value in sorted(h_1ft.items(), key=lambda item: (item[1], item[0]), reverse=True):
        if i >= topN:
            break
        print("%s\t%d" % (key, value))
        i += 1
    i = 0
    for key, value in sorted(h_1rt.items(), key=lambda item: (item[1], item[0]), reverse=True):
        if i >= topN:
            break
        print("%s\t%d" % (key, value))
        i += 1
    i = 0
    for key, value in sorted(h_2fh.items(), key=lambda item: (item[1], item[0]), reverse=True):
        if i >= topN:
            break
        print("%s\t%d" % (key, value))
        i += 1
    i = 0
    for key, value in sorted(h_2rh.items(), key=lambda item: (item[1], item[0]), reverse=True):
        if i >= topN:
            break
        print("%s\t%d" % (key, value))
        i += 1
    i = 0
    for key, value in sorted(h_2ft.items(), key=lambda item: (item[1], item[0]), reverse=True):
        if i >= topN:
            break
        print("%s\t%d" % (key, value))
        i += 1
    i = 0
    for key, value in sorted(h_2rt.items(), key=lambda item: (item[1], item[0]), reverse=True):
        if i >= topN:
            break
        print("%s\t%d" % (key, value))
        i += 1

    print("Number of read pairs = %d, Number of unqualified read pairs=%d, Number of qualified reads=%d (%.2f%%)" % (
        num_pairs, num_pairs - num_qualified, num_qualified, float(num_qualified * 100) / num_pairs))
    print("Number of clipped reads in read1 and read2:\t%d\t%d" % (
        num_softclip1, num_softclip2))
    print("Number of head-clipped reads in read1 and read2:\t%d\t%d" % (
        num_sc1head, num_sc2head))
    print("Number of tail-clipped reads in read1 and read2:\t%d\t%d" % (
        num_sc1tail, num_sc2tail))

    print("Table 1 - Insertion size distribution for Perfectly Mapped reads and soft clipping reads")
    print("insert size\tPerfectly Mapped\tSoft Clipping")
    for i in range(MAX_INSERT_SIZE):
        print("%d\t%d\t%d" % (i, insertsize[i], insertsize4sc[i]))

    print("Table 2 - Quality distribution of Read 1 and Read 2:")
    print("Quality\tProper_Mapped_Read1\tProper_Mapped_Read2\tALL_Read1\tALL_Read2")
    for i in range(100):
        print("%d\t%d\t%d\t%d\t%d" % (i, qual1[i], qual2[i], qualALL1[i], qualALL2[i]))

    print("Table 3 - Length distribution of soft clipping in Read 1 and Read 2:")
    print("Length\tRead1\tRead2")
    for i in range(200):
        print("%d\t%d\t%d" % (i, sc1[i], sc2[i]))

    output.close()
    samfile.close()
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
        ofile = "%s.unmapped.sam" % ifile

    # Main Function
    analyzer(ifile, ofile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
