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
# @file    BAMPileupScan.py
#
# @brief   Diff the coverage of each region in genome by BAMPileupScan.py
#
# @author  Chung-Tsai Su(chungtsai_su@atgenomix.com)
#
# @date    2019/02/18
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os.path

# CONST
chrID = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
         'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
         'chr21', 'chr22', 'chrX', 'chrY']
chrSize = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422,
           135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167,
           46709983, 50818468, 156040895, 57227415]

# Default Parameter
DEBUG = 0
OUTPUT_FORMAT = 1   # 0:IGV format; 1:tsv


def usage():
    print("BAMPileupScan.py -i <Input Pileup Result 1> -j <Input Pileup Result 2> "
          "-o <Output missing regions of Input file 1>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: Input Pileup Result 1")
    print("\t-j: Input Pileup Result 2")
    print("\t-o: Output missing regions of Input file 1")
    print("Usage:")
    print("\ttime python ./BAMPileupDiff.py -i ~/1.txt -j ~/2.txt -o ~/3.log")

    return


def pileup_diff(ifn1, ifn2, ofn):
    missing2 = []
    ifd2 = open(ifn2, "r")
    for line in ifd2:
        items = line.strip().split("\t")
        missing2.append(items)
    ifd2.close()
    cnt = len(missing2)
    idx = 0
    num_regions = 0
    num_bps = 0
    total_regions = 0
    total_bps = 0
    pre_chr = ""
    ofd = open(ofn, "w")
    ifd1 = open(ifn1, "r")
    for line in ifd1:
        items = line.strip().split("\t")
        if DEBUG:
            print("READ:%s, pre_chr=%s, idx=%d, missing2=%s" % (items, pre_chr, idx, missing2[idx]))

        while idx < cnt and missing2[idx][0] != items[0] and pre_chr == missing2[idx][0]:
            if DEBUG:
                print("[1]idx++")
            idx += 1

        if pre_chr == "":
            pre_chr = items[0]
        elif pre_chr != items[0]:
            print("%s\t%d\t%d" % (pre_chr, num_bps, num_regions))
            total_regions += num_regions
            total_bps += num_bps
            num_regions = num_bps = 0
            pre_chr = items[0]

        if DEBUG and idx < cnt:
            print("idx=%d, cnt=%d, items=%s, missing2[idx]=%s" %(idx, cnt, items, missing2[idx]))
        if idx >= cnt or missing2[idx][0] != items[0]:
            # all missings can be covered by ifd2
            if DEBUG:
                print("[1]%s:%s-%s\t%s\n" % (items[0], items[1], items[2], items[3]))
            if OUTPUT_FORMAT == 0:
                ofd.write("%s:%s-%s\t%s\n" % (items[0], items[1], items[2], items[3]))
            else:
                ofd.write("%s\t%s\t%s\t%s\n" % (items[0], items[1], items[2], items[3]))
            num_regions += 1
            num_bps += int(items[2]) - int(items[1]) + 1
        else:
            while idx < cnt and missing2[idx][0] != items[0]:
                if DEBUG:
                    print("[2]idx++")
                idx += 1
            while idx < cnt and missing2[idx][0] == items[0] and int(missing2[idx][2]) < int(items[1]):
                if DEBUG:
                    print("[3]idx++")
                idx += 1
            if idx >= cnt:
                if DEBUG:
                    print("[2]%s:%s-%s\t%s" % (items[0], items[1], items[2], items[3]))
                if OUTPUT_FORMAT == 0:
                    ofd.write("%s:%s-%s\t%s\n" % (items[0], items[1], items[2], items[3]))
                else:
                    ofd.write("%s\t%s\t%s\t%s\n" % (items[0], items[1], items[2], items[3]))
                num_regions += 1
                num_bps += int(items[2]) - int(items[1]) + 1
            elif int(missing2[idx][1]) > int(items[2]):
                if DEBUG:
                    print("[3]%s:%s-%s\t%s\n" % (items[0], items[1], items[2], items[3]))
                if OUTPUT_FORMAT == 0:
                    ofd.write("%s:%s-%s\t%s\n" % (items[0], items[1], items[2], items[3]))
                else:
                    ofd.write("%s\t%s\t%s\t%s\n" % (items[0], items[1], items[2], items[3]))
                num_regions += 1
                num_bps += int(items[2]) - int(items[1]) + 1
            elif int(missing2[idx][2]) >= int(items[2]):
                if int(missing2[idx][1]) > int(items[1]):
                    if DEBUG:
                        print("[4]%s:%s-%d\t%s\n" % (items[0], items[1], int(missing2[idx][1])-1, int(missing2[idx][1]) - int(items[1])))
                    if OUTPUT_FORMAT == 0:
                        ofd.write("%s:%s-%d\t%s\n" % (items[0], items[1], int(missing2[idx][1])-1, int(missing2[idx][1]) - int(items[1])))
                    else:
                        ofd.write("%s\t%s\t%d\t%s\n" % (items[0], items[1], int(missing2[idx][1]) - 1, int(missing2[idx][1]) - int(items[1])))
                    num_regions += 1
                    num_bps += int(missing2[idx][1]) - int(items[1])
            else:
                start = int(items[1]) + 1
                while idx < cnt and int(missing2[idx][2]) <= int(items[2]) and missing2[idx][0] == items[0]:
                    if DEBUG:
                        print("start: %d\t%s" %(idx, missing2[idx]))
                    if int(missing2[idx][1]) > int(items[1]):
                        if DEBUG:
                            print("[5]%s:%s-%d\t%d\n" % (items[0], start, int(missing2[idx][1])-1, int(missing2[idx][1]) - start))
                        if OUTPUT_FORMAT == 0:
                            ofd.write("%s:%s-%d\t%d\n" % (items[0], start, int(missing2[idx][1])-1, int(missing2[idx][1]) - start))
                        else:
                            ofd.write("%s\t%s\t%d\t%d\n" % (
                                items[0], start, int(missing2[idx][1]) - 1, int(missing2[idx][1]) - start))
                        num_regions += 1
                        num_bps += int(missing2[idx][1]) - start
                    elif missing2[idx][0] != items[0]:
                        if DEBUG:
                            print("[6]%s:%s-%s\t%d\n" % (items[0], start, items[2], int(items[2]) - start + 1))
                        if OUTPUT_FORMAT == 0:
                            ofd.write("%s:%s-%s\t%d\n" % (items[0], start, items[2], int(items[2]) - start + 1))
                        else:
                            ofd.write("%s\t%s\t%s\t%d\n" % (items[0], start, items[2], int(items[2]) - start + 1))
                        num_regions += 1
                        num_bps += int(items[2]) - start + 1
                    start = int(missing2[idx][2]) + 1
                    idx += 1

    ifd1.close()
    ofd.close()
    print("%s\t%d\t%d" % (pre_chr, num_bps, num_regions))
    total_regions += num_regions
    total_bps += num_bps
    print("Total\t%d\t%d" % (total_bps, total_regions))
    return


def main(argv):
    ifile1 = ""
    ifile2 = ""
    ofile = ""

    try:
        opts, args = getopt.getopt(argv, "hi:j:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt == "-i":
            ifile1 = arg
        elif opt == "-j":
            ifile2 = arg
        elif opt == "-o":
            ofile = arg

    # error handling for input parameters
    if ifile1 == "":
        print("Error: '-i' is required")
        usage()
        sys.exit(2)
    elif not os.path.isfile(ifile1):
        print("Error: input file(%s) is not existed" % ifile1)
        usage()
        sys.exit(3)
    elif ifile2 == "":
        print("Error: '-j' is required")
        usage()
        sys.exit(4)
    elif not os.path.isfile(ifile1):
        print("Error: input file(%s) is not existed" % ifile2)
        usage()
        sys.exit(5)

    if ofile == "":
        ofile = "%s.diff.log" % ifile1

    # Main Function
    pileup_diff(ifile1, ifile2, ofile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
