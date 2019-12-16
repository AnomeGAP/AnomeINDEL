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
# @file    CloudBlast.py
#
# @brief   Blast Query via Web
#
# @author  Chung-Tsai Su(chungtsai.su@atgenomix.com)
#
# @date    2019/07/05
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os.path
from Bio.Blast import NCBIWWW   # https://biopython-cn.readthedocs.io/zh_CN/latest/cn/chr07.html#sec-parsing-blast
from Bio.Blast import NCBIXML

# CONST
E_VALUE_THRESH = 0.000001

# Default Parameter


def usage():
    print("CloudBlast.py -i <Input FASTA file> -o <Output BLAST result>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: Input FASTA file")
    print("\t-o: Output BLAST result")
    print("Usage:")
    print("\tpython3 ./CloudBlast.py -i ../data/NA12878/U1000.alt.fa  "
          "-o ../data/NA12878/U1000.alt.fa.blast ")

    return


def query_fun(ifile, ofile):
    i = 0
    (num_both, num_bac, num_nui, num_unknown) = (0, 0, 0, 0)
    ofd = open(ofile, "w")
    ifd = open(ifile, "r")
    ofd.write("QUERY\tLength\tBIT-Score\tCONTIG\tLength\t#Matches\n")
    for item in ifd:
        item = item.strip()
        if item.startswith(">"):
            name = item
        else:
            query = "%s\n%s" % (name, item)
            hit_bac = 0
            hit_nui = 0
            result = NCBIWWW.qblast("blastn", "nt", query)
            blast_record = NCBIXML.read(result)
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    # print('****Alignment****')
                    # print(' sequence:', alignment.title)
                    # print('length:', alignment.length)
                    # print('e value:', hsp.expect)
                    # print('score: ', hsp.score)
                    # print('bits: ', hsp.bits)
                    # print('identities: ', hsp.identities)
                    ofd.write("%s\t%d\t%f\t%s\t%d\t%d\n" % (name[1:], len(item), hsp.bits, alignment.title,
                                                            alignment.length, hsp.identities))
                    if hsp.expect < E_VALUE_THRESH:
                        if str(alignment.title).find("Homo sapiens BAC") > -1:
                            hit_bac = 1
                            print("%s\t%d\t%f\t%s\t%d\t%d" % (name[1:], len(item), hsp.bits, alignment.title,
                                                              alignment.length, hsp.identities))
                        elif str(alignment.title).find("non-reference unique insertion sequence") > -1:
                            hit_nui = 1
                            print("%s\t%d\t%f\t%s\t%d\t%d" % (name[1:], len(item), hsp.bits, alignment.title,
                                                              alignment.length, hsp.identities))
                    break
            if hit_bac and hit_nui:
                num_both += 1
                print("%s\tBOTH" % (name[1:]))
            elif hit_bac:
                num_bac += 1
            elif hit_nui:
                num_nui += 1
            else:
                num_unknown += 1
                print("%s\tUNKNOWN" % (name[1:]))
        i += 1
        # if i > 3:
        #     break
    print("BAC+NUI\tBAC\tNUI\tUnknown")
    print("%d\t%d\t%d\t%d" % (num_both, num_bac, num_nui, num_unknown))
    ifd.close()
    ofd.close()
    return


def main(argv):
    ifile = ""
    ofile = ""

    try:
        opts, args = getopt.getopt(argv, "hi:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt == "-i":
            ifile = arg
        elif opt == "-o":
            ofile = arg

    if ofile == "":
        ofile = "%s.blast.log" % ifile

    # error handling for input parameters
    if ifile == "" :
        print("Error: '-i' is required")
        usage()
        sys.exit(2)
    elif not os.path.isfile(ifile):
        print("Error: input file(%s) is not existed" % ifile)
        usage()
        sys.exit(3)

    # Main Function
    query_fun(ifile, ofile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
