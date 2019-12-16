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
# @file    BAMExtractor.py
#
# @brief   Parsing BAM file to extract long insertion or deletions for Connected-Reads
#
# @author  Chung-Tsai Su(chungtsai.su@atgenomix.com)
#
# @date    2019/05/23
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
DIVISOR = 50
NUM = 100


def usage():
    print("BAMExtractor.py -i <Input BAM> -t [I|D|S|H|U] -l <Minimal Length> -o <Output tsv list>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: Input BAM file")
    print("\t-t: [I|D|S|H|U] (default: I )")
    print("\t-l: minimal length (default: 200 )")
    print("\t-o: Output tsv list")
    print("Usage:")
    print("\ttime python3 ./BAMExtractor.py -i ../data/NA12878/result-1.0.2-qual-fix-6.primary_alt.sorted.bam -t I -l 50 -o ../data/NA12878/I50.tsv")
    print("\ttime python3 ./BAMExtractor.py -i ../data/NA12878/result-1.0.2-qual-fix-6.primary_alt.sorted.bam -t U -l 200 -o ../data/NA12878/U200.tsv")
    print("\ttime python3 ./BAMExtractor.py -i ~/NA12878-novaseq/v1.0.2/result-1.0.2-qual-fix-6.primary_alt.sorted.bam -t D -l 200 -o ~/NA12878-novaseq/v1.0.2/longDeletion.tsv")
    print("\ttime python3 ./BAMExtractor.py -i ~/NA12878-novaseq/v1.0.2/result-1.0.2-qual-fix-6.primary_alt.sorted.bam -t U -l 200 -o ~/NA12878-novaseq/v1.0.2/unmapped.fa > ~/NA12878-novaseq/v1.0.2/unmapped.tsv")

    return


def analyzer(ifn, vtype, min_length, ofn):
    samfile = pysam.AlignmentFile(ifn, "rb")
    output = open(ofn, "w")

    print("#mapped-reads\t%d" % samfile.mapped)
    print("#unmapped-reads\t%d" % samfile.unmapped)
    print("Unmapped-Ratio\t%.2f%%" % (float(samfile.unmapped * 100) / (samfile.mapped + samfile.unmapped)))
    print("#nocoordinate\t%d" % samfile.nocoordinate)
    print("")
    num_contigs = 0
    num_hits = 0
    count = [0] * NUM
    max_hit = 0

    if vtype == "U":
        # unmapped reads
        for read in samfile.fetch(contig="*"):
            num_contigs += 1
            if read.is_unmapped:
                l = len(read.seq)
                if max_hit < l:
                    max_hit = l
                if l >= DIVISOR * NUM:
                    count[NUM - 1] += 1
                else:
                    count[int(l / DIVISOR)] += 1
                if l >= min_length:
                    output.write(">%s\n" % read.query_name)
                    output.write("%s\n" % read.seq)
                    num_hits += 1

    else:
        # mapped reads
        for read in samfile.fetch():
            if str(read.cigarstring).find(vtype) >= 0:
                # print("cigar=%s" % read.cigarstring)
                items = re.split("([0-9]+\S)", str(read.cigarstring))
                # print(items)
                idx = 0

                for i in range(len(items)):
                    # print("%d:[%s]" % (i, items[i]))
                    if items[i] != "" and vtype in items[i]:
                        l = int(items[i].split(vtype)[0])
                        if max_hit < l:
                            max_hit = l
                        if l >= DIVISOR * NUM:
                            count[NUM-1] += 1
                        else:
                            count[int(l/DIVISOR)] += 1
                        if l >= min_length:
                            num_hits += 1
                            # print(read.reference_name)
                            # print(read.reference_id)
                            # print(read.reference_start)
                            # print(read.cigarstring)
                            # print(read.seq)
                            if read.seq:
                                output.write("%s\t%s\t%d\t%s\t%d\t%s\n" % (read.query_name, read.reference_name, read.reference_start, read.cigarstring, l, read.seq[idx:idx+l]))
                                num_contigs += 1
                            # else:
                            #     output.write("%s\t%s\t%d\t%s\t%d\n" % (read.query_name, read.reference_name, read.reference_start, read.cigarstring, l))
                        idx += l

    print("#Matched\t%d" % num_contigs)
    print("#Hits\t%d" % num_hits)
    print("MaxLength\t%d" % max_hit)
    print("")
    print("Length\tAmount")
    for i in range(NUM):
        print("%d\t%d" % (i, count[i]))
    output.close()
    samfile.close()
    return


def main(argv):
    ifile = ""
    ofile = ""
    vtype = "I"
    min_length = 200

    try:
        opts, args = getopt.getopt(argv, "hi:t:l:o:")
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
            vtype = arg
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
    analyzer(ifile, vtype, min_length, ofile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
