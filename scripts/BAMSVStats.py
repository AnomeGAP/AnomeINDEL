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
MIN_MAPQ = 0


def usage():
    print("BAMSVStats.py -i <Input BAM> -l <Minimal Length> -q <Minimal MAPQ> -o <Output tsv list>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: Input BAM file")
    print("\t-l: minimal length (default: %d )" % MIN_LENGTH)
    print("\t-q: minimal MAPQ (default: %d )" % MIN_MAPQ)
    print("\t-o: Output tsv list")
    print("Usage:")
    print("\ttime python3 ./BAMSVStats.py -i ~/NA12878-novaseq/v1.0.2/result-1.0.2.primary.sorted.bam -l %d -q %d "
          "-o ~/NA12878-novaseq/v1.0.2/SV.tsv > ~/NA12878-novaseq/v1.0.2/SV.log" % (MIN_LENGTH, MIN_MAPQ))
    print("\ttime python3 ./BAMSVStats.py -i ../data/NA12878/result-1.0.2-qual-fix-6.primary.sorted.bam -l %d -q %d "
          "-o ../data/NA12878/SV-L%d-Q%d.tsv > ../data/NA12878/SV-L%d-Q%d.log" % (MIN_LENGTH, MIN_MAPQ, MIN_LENGTH,
                                                                                  MIN_MAPQ, MIN_LENGTH, MIN_MAPQ))
    print("\ttime python3 ./BAMSVStats.py -i ~/NA12878-novaseq/v1.0.2/result-1.0.2-qual-fix-6.primary_alt.sorted.bam "
          "-l 50 -q 1 -o ~/NA12878-novaseq/v1.0.2/SV-L50-Q1.tsv > ~/NA12878-novaseq/v1.0.2/SV-L50-Q1.SVStats.log")
    return

#     | POS | SEQ
# ----+-----+-----
#   M |  +  |  +
#   I |     |  +
#   D |  +  |
#   S |     |  +
#   H |     |


def analyzer(ifn, min_length, min_mapq, ofn):
    samfile = pysam.AlignmentFile(ifn, "rb")
    output = open(ofn, "w")
    output_fa = open(ofn + "_insertion.q" + str(min_mapq) + ".fa", "w")

    total_contigs = 0
    num_insertions = 0
    num_deletions = 0
    num_others = 0
    num_cnv = 0
    num_tra = 0
    i_count = [0] * NUM
    d_count = [0] * NUM
    hs_count = [0] * NUM
    h_hits = defaultdict(int)
    h_ins = defaultdict(int)
    h_del = defaultdict(int)
    h_cli = defaultdict(int)
    max_hit = 0
    max_pos = ""

    # mapped reads for Insertion
    for read in samfile.fetch():
        if read.reference_name == 'chrM':
            break

        # if any(re.findall(r'I|D|H|S', str(read.cigarstring), re.IGNORECASE)):
        if any(re.findall(r'I', str(read.cigarstring), re.IGNORECASE)) and read.mapping_quality >= min_mapq:
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
                        output.write("I\t%s\t%s\t%d\t%d\t%d\n" % (read.query_name, read.reference_name,
                                                                  int(read.reference_start) + idx,
                                                                  read.mapping_quality, l))
                        if read.query_sequence:
                            output_fa.write(">%s %d %s:%d\n%s\n" % (read.query_name, l, read.reference_name,
                                                                    read.reference_start + idx,
                                                                    read.query_sequence[idx:(idx+l)]))
                        h_ins[read.reference_name] += 1
                        if l >= DIVISOR * NUM:
                            i_count[NUM-1] += 1
                        else:
                            i_count[int(l/DIVISOR)] += 1
                        h_hits[key] += 1

                if "M" in items[i] or "D" in items[i]:
                        idx += l

            if is_matched > 0:
                num_insertions += 1
        total_contigs += 1

    output_fa.close()

    print("There are %d contigs" % total_contigs)
    print("There are %d long insertions" % num_insertions)
    print("There are %d distinct long insertions" % len(h_hits))
    print("Maximal length of hit : %d at %s" % (max_hit, max_pos))

    output_fa = open(ofn + "_deletion.q" + str(min_mapq) + ".header", "w")
    max_hit = 0
    # mapped reads for Deletion
    for read in samfile.fetch():
        if read.reference_name == 'chrM':
            break

        if any(re.findall(r'D', str(read.cigarstring), re.IGNORECASE)) and read.mapping_quality >= min_mapq:
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
                        output.write("D\t%s\t%s\t%d\t%d\t%d\n" % (read.query_name, read.reference_name,
                                                                  int(read.reference_start) + idx,
                                                                  read.mapping_quality, l))
                        output_fa.write(">%s %d %s:%d\n" % (read.query_name, l, read.reference_name,
                                                                read.reference_start + idx))
                        h_del[read.reference_name] += 1
                        if l >= DIVISOR * NUM:
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
    output_fa.close()

    h_clipping = defaultdict(str)
    h_clipping_pos = defaultdict(str)
    h_del2 = defaultdict(int)

    output_tra = open(ofn + "_tra.q" + str(min_mapq) + ".tsv", "w")
    output_fa = open(ofn + "_softclipping.q" + str(min_mapq) + ".fa", "w")
    max_hit = 0
    # mapped reads for Deletion
    for read in samfile.fetch():
        if read.reference_name == 'chrM':
            break

        if any(re.findall(r'H|S', str(read.cigarstring), re.IGNORECASE)) and read.mapping_quality >= min_mapq:
            # print("%s" % read.cigarstring)
            items = re.split("([0-9]+\S)", str(read.cigarstring))
            idx = 0
            is_matched = 0
            is_tra = 0

            for i in range(len(items)):
                if items[i] == "":
                    continue
                l = int(items[i][:-1])
                # print("%d\t%d\t%s" % (len(items), i, items[i]))
                #if ("H" in items[i] or "S" in items[i]) and l >= min_length and i == (len(items) - 2):
                if ("H" in items[i] or "S" in items[i]) and l >= min_length:
                    # print("YES\t%d\t%s" % (i, items[i]))
                    is_matched += 1
                    if h_clipping[read.query_name] == "":
                        key = "%s:%d" % (read.reference_id, int((read.reference_start + idx) / DIVISOR))
                        h_del2[read.query_name] = int(read.reference_start) + idx
                        if max_hit < l:
                            max_hit = l
                            max_pos = "%s\t%s:%d" % (read.query_name, read.reference_name, read.reference_start)
                        if h_hits[key] == 0:
                            output.write("H|S\t%s\t%s\t%d\t%d\t%d\n" % (read.query_name, read.reference_name,
                                                                        int(read.reference_start) + idx,
                                                                        read.mapping_quality, l))
                            if "S" in items[i] and read.query_sequence:
                                output_fa.write(">%s %d %s:%d\n%s\n" % (read.query_name, l, read.reference_name,
                                                                        read.reference_start + idx,
                                                                        read.query_sequence[idx:(idx + l)]))

                            h_cli[read.reference_name] += 1
                            if l >= DIVISOR * NUM:
                                hs_count[NUM - 1] += 1
                            else:
                                hs_count[int(l / DIVISOR)] += 1
                            h_hits[key] += 1
                    elif h_clipping[read.query_name] != read.reference_name:
                        is_tra = 1
                    elif h_del2[read.query_name] != read.reference_start + idx:
                        output.write("DEL2\t%s\t%s\t%d\t%d\n" % (read.query_name, read.reference_name,
                                                                 h_del2[read.query_name], read.reference_start + idx))
                        h_del2[read.query_name] = read.reference_start + idx

                if "M" in items[i] or "D" in items[i]:
                    idx += l

            if is_tra > 0:
                num_tra += 1
                output_tra.write("%s\t%s\t%d\t%s:%d\t%s\n" % (read.query_name, h_clipping_pos[read.query_name],
                                                              read.mapping_quality, read.reference_name,
                                                              read.reference_start, str(read.cigarstring)))
                output.write("TRA\t%s\t%s\t%s\t%d\n" % (read.query_name, h_clipping[read.query_name],
                                                        read.reference_name, read.reference_start))
            if is_matched > 0:
                h_clipping[read.query_name] = read.reference_name
                h_clipping_pos[read.query_name] = "%d\t%s:%d\t%s" % (read.mapping_quality, read.reference_name, read.reference_start,
                                                                     str(read.cigarstring))
                num_others += 1
                if is_matched > 1:
                    num_cnv += 1
                    output.write("CNV1\t%s\t%s\t%d\n" % (read.query_name, read.reference_name, read.reference_start))

        total_contigs += 1

    print("There are %d long num_others" % num_others)
    print("There are %d potential CNV" % num_cnv)
    print("There are %d potential TRA" % num_tra)
    print("There are %d distinct long insertions + long deletions + long others" % len(h_hits))
    print("Maximal length of hit : %d at %s" % (max_hit, max_pos))
    output_fa.close()
    output_tra.close()

    print("CHR\t#INS\t#DEL\t#CLIPPING")
    for i in h_ins.keys():
        print("%s\t%d\t%d\t%d" % (i, h_ins[i], h_del[i], h_cli[i]))

    print("\n")
    for i in range(NUM):
        print("%d\t%d\t%d\t%d" % (i, i_count[i], d_count[i], hs_count[i]))
    output.close()
    samfile.close()
    return


def main(argv):
    ifile = ""
    ofile = ""
    min_length = MIN_LENGTH
    min_mapq = MIN_MAPQ

    try:
        opts, args = getopt.getopt(argv, "hi:l:q:o:")
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
        elif opt == "-q":
            min_mapq = int(arg)
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
    analyzer(ifile, min_length, min_mapq, ofile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
