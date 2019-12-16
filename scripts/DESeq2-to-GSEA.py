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


def usage():
    """ Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    """
    print("DESeq2-to-GSEA.py -i <DESeq2 normalized readcount file> -j <DESeq2 phenotype condition file> -o <GSEA input prefixs>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <DESeq2 readcount file>")
    print("\t-j: <DESeq2 phenotype condition file>")
    print("\t-o: <GSEA input prefixs>")
    print("Usage:")
    print("\tpython ./DESeq2-to-GSEA.py -i ~/ballgown-P2C/normalized_read_counts.txt -j ~/ballgown-P2C/condition.list -o ~/ballgown-P2C/P2C")
    print("\tpython ./DESeq2-to-GSEA.py -i ~/ballgown-P2S/normalized_read_counts.txt -j ~/ballgown-P2S/condition.list -o ~/ballgown-P2S/P2S")
    print("\tpython ./DESeq2-to-GSEA.py -i ~/ballgown-C2S/normalized_read_counts.txt -j ~/ballgown-C2S/condition.list -o ~/ballgown-C2S/C2S")

    return


def blocks(files, size=65536):
    while True:
        b = files.read(size)
        if not b: break
        yield b


def transform(ifn, jfn, oprefix):
    sample_count = 0
    phenotypefd = open(oprefix + ".cls", "w")
    jfd = open(jfn, "r")
    for line in jfd:
        items = line.strip().split(" ")
        sample_count = len(items)
        phenotypefd.write("%d 2 1\n" % sample_count)
        temp = set()
        for item in items:
            if item not in temp:
                if len(temp) == 0:
                    phenotypefd.write("#%s" % item)
                else:
                    phenotypefd.write(" %s" % item)
                temp.add(item)
        phenotypefd.write("\n")
        phenotypefd.write("%s" % line)
    jfd.close()
    phenotypefd.close()

    gene_count = -1
    # with open(ifn, "r") as f:
    #     for bl in blocks(f):
    #         gene_count += bl.count("\n")
    ifd = open(ifn, "r")
    for line in ifd:
        items = line.strip().split("\t")
        if items[0].startswith("MSTRG"):
            continue
        gene_count += 1
    ifd.close()

    is_firstline = 1
    gctfd = open(oprefix + ".gct", "w")
    ifd = open(ifn, "r")
    for line in ifd:
        items = line.strip().split("\t")
        if items[0].startswith("MSTRG"):
            continue
        elif is_firstline == 1:
            is_firstline = 0
            gctfd.write("#Atgenomix v1.0\n%d\t%d\nNAME\tDESCRIPTION" % (gene_count, sample_count))
            for i in range(0, len(items)):
                gctfd.write("\t%s" % items[i])
        else:
            gctfd.write("%s\tna" % items[0].upper())
            for i in range(1, len(items)):
                gctfd.write("\t%s" % items[i])
        gctfd.write("\n")

    ifd.close()
    gctfd.close()

    return


def main(argv):
    ifn = ""
    jfn = ""
    output_prefix = ""

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
        elif opt in "-j":
            jfn = arg
        elif opt in "-o":
            output_prefix = arg
    if not ifn or not jfn:
        usage()
        sys.exit(2)
    elif output_prefix == "":
        print("Please specify the GSEA output prefix you want to query")
        usage()
        sys.exit(3)

    # Main Function
    transform(ifn, jfn, output_prefix)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
