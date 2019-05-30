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
# @file    GraaphSeqStats.py
#
# @brief   An statistics tool to analyze the distribution of suffix and edges by partition
#
# @author  A-Tsai Su(chungtsai_su@anome.com.tw)
#
# @date    2018/4/20
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os.path
import collections
from collections import defaultdict

def Usage():
    ''' Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    '''
    print("GraaphSeqStats.py -i <GraphSeq stats file> -n <normalization value> -o <output file>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <GraphSeq stats file>")
    print("\t-n: <normalization value> [default: 100]")
    print("\t-o: <output file>")
    print("Usage:")
    print("\t./GraaphSeqStats.py -i summary.txt -n 100 -o summary.100.out")

    return


def graaphseq_stats(ifn, ofn, normalization):

    ifd = open(ifn, "r")
    ofd = open(ofn, "w")

    shash = defaultdict(int)
    ehash = defaultdict(int)
    s_total = 0
    e_total = 0
    max_suffix = 0
    min_suffix = 10000000
    max_edges = 0
    min_edges = 10000000
    avg_suffix = 0
    avg_edges = 0
    max_suffix_raw = ""
    min_suffix_raw = ""
    max_edges_raw = ""
    min_edges_raw = ""

    for line in ifd:
        line = line.strip()
        if not line.startswith("#"):
            #print("%s\n" % line)
            items = line.split('\t')
            num_suffix = int(items[2])
            num_edges = int(items[3])

            # batchID partitionID #suffix #edges
            shash[int(num_suffix / normalization)] += 1
            ehash[int(num_edges / normalization)] += 1
            s_total += num_suffix
            e_total += num_edges
            avg_suffix += 1
            avg_edges += 1
            if max_suffix < num_suffix:
                max_suffix = num_suffix
                max_suffix_raw = line
            if max_edges < num_edges:
                max_edges = num_edges
                max_edges_raw = line
            if min_suffix > num_suffix:
                min_suffix = num_suffix
                min_suffix_raw = line
            if min_edges > num_edges:
                min_edges = num_edges
                min_edges_raw = line

    avg_suffix = s_total / avg_suffix
    avg_edges = e_total / avg_edges

    for k in sorted(shash.keys()):
        ofd.write("%d\t%d\n" % (k*normalization, shash[k]))
    ofd.write("\n")

    for k in sorted(ehash.keys()):
        ofd.write("%d\t%d\n" % (k*normalization, ehash[k]))

    ifd.close()
    ofd.close()

    print("Total suffix = %d" % s_total)
    print("MAX(suffix) = %d at %s" % (max_suffix, max_suffix_raw))
    print("MIN(suffix) = %d at %s" % (min_suffix, min_suffix_raw))
    print("AVG(suffix) = %f" % avg_suffix)

    print("Total edges = %d" % e_total)
    print("MAX(edges) = %d at %s" % (max_edges, max_edges_raw))
    print("MIN(edges) = %d at %s" % (min_edges, min_edges_raw))
    print("AVG(edges) = %f" % avg_edges)

    return


def main(argv):
    infile = ""
    outfile = "output.txt"
    normalization = 1000

    try:
        opts, args = getopt.getopt(argv, "hi:o:n:")
    except getopt.GetoptError:
        Usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            Usage()
            sys.exit()
        elif opt in "-i":
            infile = arg
        elif opt in "-n":
            normalization = int(arg)
        elif opt in "-o":
            outfile = arg

    if not infile:
        Usage()
        sys.exit(2)
    elif not os.path.isfile(infile):
        print("Error: input file(%s) is not existed" % infile)
        Usage()
        sys.exit(3)

    # Main Function
    graaphseq_stats(infile, outfile, normalization)

    return

if __name__ == '__main__':
    main(sys.argv[1:])
