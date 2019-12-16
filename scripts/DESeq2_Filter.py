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
# @file    DESeq2_Filter.py
#
# @brief   Filter genes from DESeq2
#
# @author  Chung-Tsai Su(chungtsai.su@atgenomix.com)
#
# @date    2019/10/23
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os.path
from collections import defaultdict


def usage():
    """ Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    """
    print("DESeq2_Filter.py -i <DESeq2 normalized readcount file>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <DESeq2 readcount file>")
    print("Usage:")
    print("\tpython ./DESeq2-to-GSEA.py -i ~/ballgown-ALL/gene_count_matrix.csv ")

    return


def blocks(files, size=65536):
    while True:
        b = files.read(size)
        if not b: break
        yield b


def filtering(ifn):
    cnt_up = 0
    cnt_down = 0

    ofd = open(ifn+"_epithelial.csv", "w")
    ifd = open(ifn, "r")
    ofd.write("Gene\tParental\tPrimary2\tPrimary3\tCTC1\tCTC2\tSecondary1\tSecondary2\n")
    for line in ifd:
        if line.startswith("gene_id") or line.startswith("MSTRG"):
            continue

        items = line.strip().split(",")
        # Samples: Parental, Primary2, Primary3, CTC1, CTC2, Secondary1, Secondary2
        # idx:         1         3        4        5    6       8            9
        # Parental > Secondary > CTC
        # Primary > Secondary > CTC
        Parental = int(items[1])
        P2 = int(items[3])
        P3 = int(items[4])
        C1 = int(items[5])
        C2 = int(items[6])
        S1 = int(items[8])
        S2 = int(items[9])
        if Parental > S1 and Parental > S2 and P2 > S1 and P3 > S1 and P2 > S2 and P3 > S2 \
        and S1 > C1 and S2 > C2 and S1 > C2 and S2 > C1 and (S1+S2) > 2*(C1+C2) and (P2+P3) > 2*(S1+S2):
            ofd.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" % (items[0].upper(), Parental, P2, P3, C1, C2, S1, S2))
            cnt_up += 1
    ifd.close()
    ofd.close()

    ofd = open(ifn+"_mesenchymal.csv", "w")
    ifd = open(ifn, "r")
    ofd.write("Gene\tParental\tPrimary2\tPrimary3\tCTC1\tCTC2\tSecondary1\tSecondary2\n")
    for line in ifd:
        if line.startswith("gene_id") or line.startswith("MSTRG"):
            continue

        items = line.strip().split(",")
        # Samples: Parental, Primary2, Primary3, CTC1, CTC2, Secondary1, Secondary2
        # idx:         1         3        4        5    6       8            9
        # Parental > Secondary > CTC
        # Primary > Secondary > CTC
        Parental = int(items[1])
        P2 = int(items[3])
        P3 = int(items[4])
        C1 = int(items[5])
        C2 = int(items[6])
        S1 = int(items[8])
        S2 = int(items[9])
        if Parental < S1 and Parental < S2 and P2 < S1 and P3 < S1 and P2 < S2 and P3 < S2 \
        and S1 < C1 and S2 < C2 and S1 < C2 and S2 < C1 and 2*(S1+S2) < (C1+C2) and 2*(P2+P3) < (S1+S2):
            ofd.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" % (items[0].upper(), Parental, P2, P3, C1, C2, S1, S2))
            cnt_down += 1
    ifd.close()
    ofd.close()

    print("%d Epithelial-Related genes, %d Mesenchymal-Related genes" % (cnt_up, cnt_down))

    return


def main(argv):
    ifn = ""

    try:
        opts, args = getopt.getopt(argv, "hi:j:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in "-i":
            ifn = arg
    if not ifn:
        usage()
        sys.exit(2)

    # Main Function
    filtering(ifn)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
