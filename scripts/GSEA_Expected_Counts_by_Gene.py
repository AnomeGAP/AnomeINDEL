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
# @file    DESeq2-to-GSEA.py
#
# @brief   Transform result format from DESeq2 to GSEA
#
# @author  Chung-Tsai Su(chungtsai.su@atgenomix.com)
#
# @date    2019/10/14
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os.path
from collections import defaultdict

TOP_N = 50

def usage():
    """ Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    """
    print("GSEA_Expected_Counts_by_Gene.py -i <DESeq2 normalized readcount file> -j <GSEA ranking list>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <DESeq2 normalized readcount file>")
    print("\t-j: <GSEA ranking list>")
    print("Usage:")
    print("\tpython ./GSEA_Expected_Counts_by_Gene.py -i ~/ballgown-P2O/normalized_read_counts.txt -j /Users/chungtsai_su/gsea_home/output/oct14/my_analysis.Gsea.1571046047332/ranked_gene_list_primary_versus_others_1571046047332.xls")
    print("\tpython ./GSEA_Expected_Counts_by_Gene.py -i ~/ballgown-P2O/normalized_read_counts.txt -j /Users/chungtsai_su/gsea_home/output/oct14/my_analysis.Gsea.1571045456843/ranked_gene_list_CTC_versus_others_1571045456843.xls")
    print("\tpython ./GSEA_Expected_Counts_by_Gene.py -i ~/ballgown-P2O/normalized_read_counts.txt -j /Users/chungtsai_su/gsea_home/output/oct14/my_analysis.Gsea.1571105821025/ranked_gene_list_primary_versus_CTC_1571105821025.xls")
    return


def transform(ifn, jfn):

    hExpectedCount = defaultdict(str)

    is_firstline = 1
    ifd = open(ifn, "r")
    for line in ifd:
        items = line.strip().split("\t")
        if is_firstline == 1:
            is_firstline = 0
        else:
            hExpectedCount[items[0].upper()] = line.strip()

    ifd.close()

    is_firstline = 1
    lRankedList = list()
    jfd = open(jfn, "r")
    for line in jfd:
        if is_firstline == 1:
            is_firstline = 0
        else:
            items = line.strip().split("\t")
            lRankedList.append(items[0])
    jfd.close()

    for i in range(0, TOP_N):
        if not str(lRankedList[i]).startswith("MSTRG."):
            print("%s" % str(hExpectedCount[lRankedList[i]]).upper())
    print("\n")
    for i in range(len(lRankedList)-1, len(lRankedList)-TOP_N-1, -1):
        if not str(lRankedList[i]).startswith("MSTRG."):
            print("%s" % str(hExpectedCount[lRankedList[i]]).upper())

    return


def main(argv):
    ifn = ""
    jfn = ""

    try:
        opts, args = getopt.getopt(argv, "hi:j:")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in "-i":
            ifn = arg
        elif opt in "-j":
            jfn = arg
    if not ifn or not jfn:
        usage()
        sys.exit(2)

    # Main Function
    transform(ifn, jfn)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
