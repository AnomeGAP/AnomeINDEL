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
# @file    BarcodingExtractor.py
#
# @brief   Extractor a specific barcode reads from 10X barcoding FASTQ files
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


def Usage():
    """ Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    """
    print("BarcodingExtractor.py -e: <Barcode> ")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-e: <Barcode>")
    print("Usage:")
    print("\tzcat ~/data/NA12878_WGS_v2/NA12878_WGS_v2_S1_L001_R1_001.fastq.gz | "
          "python ./BarcodingExtractor.py -e TGCTGTAG > "
          "~/data/NA12878_WGS_v2/NA12878_WGS_v2_S1_L001_R1_001.TGCTGTAG.fastq")

    return


def barcoding_extractor(ptn):
    cnt = 0
    skip = 1

    for line in sys.stdin:
        cnt += 1
        # if cnt % 100000000 == 0:
        #     print("%d" % cnt)

        if cnt % 4 == 1:
            items = str(line).strip().split(":")
            if items[len(items)-1] == ptn:
                skip = 0
            else:
                skip = 1
        if skip == 0:
            print("%s" % str(line).strip())

    return


def main(argv):
    pattern = ""

    try:
        opts, args = getopt.getopt(argv, "he:")
    except getopt.GetoptError:
        Usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            Usage()
            sys.exit()
        elif opt in "-e":
            pattern = str(arg)
    if pattern == "":
        print("Error: Please specify which barcode you want")
        Usage()
        sys.exit(2)

    # Main Function
    barcoding_extractor(pattern)

    return

if __name__ == '__main__':
    main(sys.argv[1:])
