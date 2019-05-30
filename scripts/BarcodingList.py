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
# @file    BarcodingList.py
#
# @brief   Listing 10X barcode information
#
# @author  A-Tsai Su(chungtsai_su@atgenomix.com)
#
# @date    2018/10/25
#
# @version 1.0
#
# @remark
#

import sys
import getopt
from collections import defaultdict


def Usage():
    """ Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    """
    print("BarcodingList.py -f 1")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-f: format")
    print("\t\t1: index sequence")
    print("\t\t2: header")
    print("Usage:")
    print("\tzcat ~/data/NA12878_WGS_v2/NA12878_WGS_v2_S1_L00*_I1_001.fastq.gz | python ./BarcodingList.py -f 1 > log")
    print("\tzcat ~/data/wfu/HN*/read-I1_si* | python ./BarcodingList.py -f 1 > log.wfu")

    return


def barcoding_list(fmt):
    cnt = 0
    h_barcode = defaultdict(int)
    code = ""

    for line in sys.stdin:
        cnt += 1
        if cnt % 10000000 == 0:
            print("%d" % cnt)

        if fmt == 1 and cnt % 4 == 2:
            code = str(line).strip()
            h_barcode[code] += 1
        elif fmt == 2 and cnt % 4 == 1:
            items = str(line).strip().split(":")
            # print("%d\t%s" % (len(items),items))
            code = items[len(items)-1]
            if len(code) < 8:
                print("%d\t%s" % (cnt, code))
            h_barcode[code] += 1

    for key, value in sorted(h_barcode.items(), key=lambda item: (item[1], item[0]), reverse=True):
        print("%s\t%d" % (key, value))
    return


def main(argv):
    fmt = 1

    try:
        opts, args = getopt.getopt(argv, "hf:")
    except getopt.GetoptError:
        Usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            Usage()
            sys.exit()
        elif opt in "-f":
            fmt = int(arg)

    # Main Function
    barcoding_list(fmt)

    return

if __name__ == '__main__':
    main(sys.argv[1:])
