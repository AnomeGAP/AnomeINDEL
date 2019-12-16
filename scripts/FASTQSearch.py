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
# @file    FASTQSearch.py
#
# @brief   Search FASTQ file by readname, sequence or quality
#
# @author  Chung-Tsai Su(chungtsai_su@atgenomix.com)
#
# @date    2019/06/14
#
# @version 1.0
#
# @remark
#

import sys
import getopt

TYPE_NAME = "NAME"
TYPE_SEQ = "SEQ"
TYPE_QUALITY = "QUALITY"

def usage():
    print("zcat <Input FASTQ.gz> | FASTQSearch.py -t [%s|%s|%s] -k <keyword> -o <output file>" % (TYPE_NAME, TYPE_SEQ, TYPE_QUALITY))
    print("Argument:")
    print("\t-h: Usage")
    print("\t-t: [%s|%s|%s]" % (TYPE_NAME, TYPE_SEQ, TYPE_QUALITY))
    print("\t-k: <Keyword>")
    print("\t-o: <Output file>")
    print("Usage:")
    print("\ttime zcat /seqslab/NA12878-novaseq/NA12878-novqseq_r1.fastq.gz | python /seqslab/atsai/script/FASTQSearch.py -t SEQ -k AAAAAAAAAAAAAAAAAAAAAATTTTCAAACGATTTCAGTTCCCAACAGGAGTACTATTAGATAGGGAATGAGTTAAATTTAATTTCTGTTTTCCTCCCAATAAAAAGAAGTGGATTGCAAAATGTGGTTTATGTACTTGTAATAA -o /seqslab/atsai/data/NA12878/AAAAAAAAAAAAAAAAAAAAAATTTTCAAACGATTTCAGTTCCCAACAGGAGTACTATTAGATAGGGAATGAGTTAAATTTAATTTCTGTTTTCCTCCCAATAAAAAGAAGTGGATTGCAAAATGTGGTTTATGTACTTGTAATAA.log")

    return


def search(t, keyword, ofn):
    idx = 0
    data = ""
    enabled = False

    ofd = open(ofn, "w")
    for line in sys.stdin:
        if idx % 4 == 0:
            data = line
            if t == TYPE_NAME and keyword in line.strip():
                enabled = True
        elif idx % 4 == 1:
            data += line
            if t == TYPE_SEQ and keyword in line.strip():
                enabled = True
        elif idx % 4 == 3:
            if t == TYPE_QUALITY and keyword in line.strip():
                enabled = True
            if enabled:
                ofd.write("%s%s" % (data, line))
                enabled = False
        idx += 1

    ofd.close()
    return


def main(argv):
    t = "NAME"
    keyword = ""
    ofile = ""

    try:
        opts, args = getopt.getopt(argv, "ht:k:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt == '-t':
            t = arg
        elif opt == '-k':
            keyword = arg
        elif opt == '-o':
            ofile = arg

    if t != "NAME" and t != "SEQ" and t != "QUALITY":
        print("Error: '-t %s' is undefined. Should be NAME, SEQ or QUALITY" % t)
        usage()
        sys.exit(2)

    if keyword == "":
        print("Error: '-k' is undefined")
        usage()
        sys.exit(3)

    if ofile == "":
        print("Error: '-o' is undefined")
        usage()
        sys.exit(4)

    search(t, keyword, ofile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
