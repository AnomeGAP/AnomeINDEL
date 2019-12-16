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
# @file    KallistoExtracter.py
#
# @brief   Extract TPM from Kalliso Output folder
#
# @author  Chung-Tsai Su(chungtsai_su@atgenomix.com)
#
# @date    2019/10/03
#
# @version 1.0
#
# @remark
#

im
import sys
import getopt
import os.path
from collections import defaultdict


def usage():
    """ Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    """
    print("KallistoExtractor.py -i <Kallisto output folder> -t <Transcript ID>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <Kallisto output folder>")
    print("\t-t: <Transcript ID>")
    print("Usage:")
    print("\tpython ./KallistoExtractor.py -i ../data/RNA -t ENSMUST00000001631.6")

    return


def extractor(ifn1, tid):
    h_hits = defaultdict(float)
    cmd = "grep %s %s/*/abundance.tsv" % (tid, ifn1)
    for i in os.popen(cmd):
        items = i.strip().split("\t")
        end = items[0].find("/abundance.tsv")
        sample = items[0][len(ifn1)+1:end]
        h_hits[sample] = float(items[4])

    output = tid
    for i in range(1, 11):
        output = output + "\t" + str(h_hits["sample" + str(i)])
    print("%s" % output)

    return


def main(argv):
    infile1 = ""
    tid = ""

    try:
        opts, args = getopt.getopt(argv, "hi:t:")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in "-i":
            infile1 = arg
        elif opt in "-t":
            tid = arg
    if not infile1:
        usage()
        sys.exit(2)
    elif tid == "":
        print("Please specify the transcript ID you want to query")
        usage()
        sys.exit(3)

    # Main Function
    extractor(infile1, tid)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
