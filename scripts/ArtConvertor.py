#!/bin/env python3
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
# @file    ArtConvertor.py
#
# @brief   Convert ART dataset to SeqGraph Error Correction dataset
#
# @author  A-Tsai Su(chungtsai_su@atgenomix.com)
#
# @date    2018/08/17
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

READ_LENGTH = 151
MAX_ERRORPERREAD = 2

def _usage():
    """ Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    """
    print(
        "ArtConvertor.py -i <Read1 fq file> -j <Read2 fq file> -m <Read1 aln file> -n <Read2 aln file> -s <SAM file>"
        "-o <Output prefix file>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <Read1 fq file>")
    print("\t-j: <Read2 fq file>")
    print("\t-m: <Read1 aln file>")
    print("\t-n: <Read2 aln file>")
    print("\t-s: <SAM file>")
    print("\t-o: <Output prefix file>")

    print("Usage:")
    print("\tpython3 ~/src/github/AnomeINDEL/scripts/ArtConvertor.py "
          "-i ~/data/art/100_dat1.fq -j ~/data/art/100_dat2.fq "
          "-m ~/data/art/100_dat1.aln -n ~/data/art/100_dat2.aln "
          "-s ~/data/art/100_dat.sam -o ~/data/hg19/test")

    return


def art_convertor(ifq1, ifq2, ialn1, ialn2, isam, ofn):
    ofd = open(ofn + "_error.tsv", "w")
    step = 0
    skip = 1
    idx = 0
    ref1 = ""
    ref2 = ""
    start1 = 0
    start2 = 0
    direction1 = 1
    direction2 = 1
    numMultipleErrorPerRead = 0

    ofd.write("#ID\tPOS\tREF\tALT\tREF_POS\n")
    noerrorfd = open(ofn + "_noerror.fq", "w")

    with open(ialn1, 'r') as f1, open(ialn2, 'r') as f2:
        for x, y in zip(f1, f2):
            if skip:
                if str(x).startswith('##Header End'):
                    skip = 0
                continue
            if str(x).startswith('>'):
                args1 = x.split('\t')
                args2 = y.split('\t')
                start1 = int(args1[2])
                start2 = int(args2[2])
                if args1[3] == "+":
                    direction1 = 1
                else:
                    direction1 = -1
                if args2[3] == "+":
                    direction2 = 1
                else:
                    direction2 = -1
                step = 0
            elif step == 0:
                step = 1
                ref1 = []
                ref1 += str(x).strip()
                ref2 = []
                ref2 += str(y).strip()
                seq = str(x).strip().replace("-", "")
                noerrorfd.write("@%09d %09d\n%s\n+\n%s\n" % (idx, idx, seq, "A"*len(seq)))
                seq = str(y).strip().replace("-", "")
                noerrorfd.write("@%09d %09d\n%s\n+\n%s\n" % (idx+1, idx+1, seq, "A" * len(seq)))
            else:
                step = 2
                read1 = []
                read1 += str(x).strip()
                read2 = []
                read2 += str(y).strip()
                cnt = 0
                length = len(str(x).strip())
                for i in range(length):
                    if ref1[i] != read1[i]:
                        ofd.write("%d\t%d\t%s\t%s\t%d\n" % (idx, i, ref1[i], read1[i], start1+direction1*i))
                        cnt += 1
                if cnt >= MAX_ERRORPERREAD:
                    # print("read %d has %d errors(%d)" % (idx, cnt, start1))
                    numMultipleErrorPerRead += 1
                cnt = 0
                length = len(str(y).strip())
                for i in range(length):
                    if ref2[i] != read2[i]:
                        ofd.write("%d\t%d\t%s\t%s\t%d\n" % (idx+1, i, ref2[i], read2[i], start2 + direction2 * i))
                        cnt += 1
                if cnt >= MAX_ERRORPERREAD:
                    # print("read %d has %d errors(%d)" % (idx+1, cnt, start2))
                    numMultipleErrorPerRead += 1
                idx += 2
    total = idx
    ofd.close()
    noerrorfd.close()

    print("There are %d reads having not less than %d errors" % (numMultipleErrorPerRead, MAX_ERRORPERREAD))

    ofd1 = open(ofn + "_1.fq", "w")
    ofd2 = open(ofn + "_2.fq", "w")
    ofd = open(ofn + "_paired.fq", "w")
    idx = 0
    step = 0
    str1 = ""
    str2 = ""
    with open(ifq1, 'r') as f1, open(ifq2, 'r') as f2:
        for x, y in zip(f1, f2):
            if step % 4 == 0:
                if not str(x).startswith('@'):
                    print("ERROR: %s Line %d doesn't start with @" % (ifq1, idx))

                ofd1.write("%s" % str1)
                ofd.write("%s" % str1)
                str1 = '%s %09d\n' % (str(x).strip(), idx)
                idx += 1
                if not str(y).startswith('@'):
                    print("ERROR: %s Line %d doesn't start with @" % (ifq2, idx))
                ofd2.write("%s" % str2)
                ofd.write("%s" % str2)
                str2 = '%s %09d\n' % (str(y).strip(), idx)
                idx += 1
            else:
                str1 = '%s%s' % (str1, str(x))
                str2 = '%s%s' % (str2, str(y))

            step += 1
    ofd1.write("%s" % str1)
    ofd.write("%s" % str1)
    ofd2.write("%s" % str2)
    ofd.write("%s" % str2)

    ofd.close()
    ofd1.close()
    ofd2.close()

    if total != idx:
        print("ERROR: total=%d, idx=%d" % (total, idx))
    else:
        print("Total = %d reads" % idx)

    idx = 0
    ofd = open(ofn + "_paired.sam", "w")
    sam = open(isam, "r")
    for line in sam:
        if line.startswith("@"):
            continue
        items = line.strip().split("\t")
        items[0] = str(idx)
        ofd.write("%s\n" % "\t".join(items))
        idx += 1

    ofd.close()
    sam.close()

    if total != idx:
        print("ERROR: total=%d, idx=%d (SAM file)" % (total, idx))
    else:
        print("Total = %d reads in SAM file" % idx)

    return


def main(argv):
    infq1 = ""
    infq2 = ""
    inaln1 = ""
    inaln2 = ""
    insam = ""
    outfileprefix = ""

    try:
        opts, args = getopt.getopt(argv, "hi:j:m:n:s:o:")
    except getopt.GetoptError:
        _usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            _usage()
            sys.exit()
        elif opt in "-i":
            infq1 = arg
            outfileprefix = infq1 + "_convert"
        elif opt in "-j":
            infq2 = arg
        elif opt in "-m":
            inaln1 = arg
        elif opt in "-n":
            inaln2 = arg
        elif opt in "-s":
            insam = arg
        elif opt in "-o":
            outfileprefix = arg
        elif opt in "-l":
            READ_LENGTH = int(arg)
    if not infq1 or not os.path.isfile(infq1):
        print("Error: input file(%s) is not existed" % (infq1))
        _usage()
        sys.exit(2)
    elif not infq2 or not os.path.isfile(infq2):
        print("Error: input file(%s) is not existed" % (infq2))
        _usage()
        sys.exit(3)
    elif not inaln1 or not os.path.isfile(inaln1):
        print("Error: input file(%s) is not existed" % (inaln1))
        _usage()
        sys.exit(4)
    elif not inaln2 or not os.path.isfile(inaln2):
        print("Error: input file(%s) is not existed" % (inaln2))
        _usage()
        sys.exit(5)
    elif not insam or not os.path.isfile(insam):
        print("Error: input file(%s) is not existed" % (insam))
        _usage()
        sys.exit(6)

    # Main Function
    art_convertor(infq1, infq2, inaln1, inaln2, insam, outfileprefix)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
