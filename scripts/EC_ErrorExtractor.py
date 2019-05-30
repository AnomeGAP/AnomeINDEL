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
# @file    EC_ErrorExtractor.py
#
# @brief   An format transformation for Error Correction Fastq data
#
# @author  A-Tsai Su(chungtsai_su@anome.com.tw)
#
# @date    2018/09/12
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


def _usage():
    """ Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    """
    print(
        "EC_ErrorExtractor.py -i <raw FASTQ file> -j <EC FASTQ file> -o <Output file>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <FASTQ file>")
    print("\t-j: <FASTQ file>")
    print("\t-o: <Output file>")

    print("Usage:")
    print("\tpython3 ./EC_ErrorExtractor.py -i new_chr22.fq "
          "-j new_chr22.ec.fq -o new_chr22_sga_error.tsv")

    return


def cal_ec(ifn, i2fn, ofn):
    step = 0
    idx = 0
    num_corrected = 0
    ofd = open(ofn, "w")
    # write header
    ofd.write("#ID\tPOS\tERROR\tCORRECT\tREF_POS\n")

    with open(ifn, 'r') as f1, open(i2fn, 'r') as f2:
        for x, y in zip(f1, f2):
            if step % 4 == 1:
                xseq = list(str(x))
                yseq = list(str(y))
                for i in range(len(xseq)):
                    if xseq[i] != yseq[i]:
                        ofd.write("%d\t%d\t%s\t%s\t%d\n" % (idx, i, xseq[i], yseq[i], -1))
                        num_corrected += 1
                idx += 1

            step += 1
        f1.close()
        f2.close()

    print("Total=%d" % idx)
    print("Corrected=%d" % num_corrected)

    ofd.close()

    return


def main(argv):
    infile = ""
    in2file = ""
    outfile = ""

    try:
        opts, args = getopt.getopt(argv, "hi:j:o:")
    except getopt.GetoptError:
        _usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            _usage()
            sys.exit()
        elif opt in "-i":
            infile = arg
            outfile = infile + "_error.tsv"
        elif opt in "-j":
            in2file = arg
        elif opt in "-o":
            outfile = arg
    if not infile or not os.path.isfile(infile):
        print("Error: input file(%s) is not existed" % infile)
        _usage()
        sys.exit(2)

    # Main Function
    cal_ec(infile, in2file, outfile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
