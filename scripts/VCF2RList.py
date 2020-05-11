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
# @file    VCF2RList.py
#
# @brief   transform variants from VCF file to R List
#
# @author  Chung-Tsai Su(chungtsai.su@atgenomix.com)
#
# @date    2020/05/11
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os
from collections import defaultdict

# CONSTANT
UNIT = 1000000
COLOR = ("#000000", "#000000", "#330000", "#990000", "#CC0000", "#FF0000")


def usage():
    print("VCF2RList.py -i <Input vcf file> -n <Sample Name>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: input file (vcf.gz)")
    print("\t-n: <Sample Name>")
    print("Usage:")
    print("\tpython ./VCF2RList.py -i ../data/WES/new-BALB-contam/Y02.PASS.vcf -n CDL1-TUMOR")
    print("\tpython ./VCF2RList.py -i ../data/WES/new-BALB-contam/Y02.PASS.vcf -n CDL1-TUMOR")
    return


def transform(ifile, sname):
    ifd = open(ifile, "rt")
    idx = 0
    str_chrom = "%s_chromsomes =c('1'" % sname
    str_coord = "%s_coordinates= c(0" % sname
    str_value = "%s_values = c(5" % sname
    str_color = "%s_colors = c(\"#AAAAAA\"" % sname
    max_count = 0
    is_first = 0

    h_vcf = defaultdict(str)
    for line in ifd:
        if not line.startswith("#") and not line.startswith("CHROM"):
            items = line.strip().split("\t")
            if items[0] == "X" or items[0] == "Y":
                h_vcf["%s_%10d" % (items[0]*2, int(items[1]))] = line.strip()
            else:
                h_vcf["%2s_%10d" % (items[0], int(items[1]))] = line.strip()
    variants = h_vcf.keys()
    variants.sort()

    pre_chrom = ""
    pre_coord = 0
    cur_chrom = ""
    cur_coord = 0
    count = 0

    for key in variants:
        if not line.startswith("#") and not line.startswith("CHROM"):
            items = h_vcf[key].split("\t")
            cur_chrom = items[0]
            cur_coord = int(int(items[1]) / UNIT)
            if pre_chrom == cur_chrom:
                if pre_coord == cur_coord:
                    count += 1
                    if count > max_count:
                        max_count += 1
                else:
                    if is_first and count:
                        str_chrom += "'%s'" % pre_chrom
                        str_coord += "%d" % (pre_coord * UNIT)
                        str_value += "%d" % count
                        str_color += "\"%s\"" % COLOR[count]
                        is_first = 0
                    elif count:
                        str_chrom += ", '%s'" % pre_chrom
                        str_coord += ", %d" % (pre_coord * UNIT)
                        str_value += ", %d" % count
                        str_color += ", \"%s\"" % COLOR[count]
                    count = 1
                    pre_coord = cur_coord
            else:
                if is_first and count:
                    str_chrom += "'%s'" % pre_chrom
                    str_coord += "%d" % (pre_coord * UNIT)
                    str_value += "%d" % count
                    str_color += "\"%s\"" % COLOR[count]
                    is_first = 0
                elif count:
                    str_chrom += ", '%s'" % pre_chrom
                    str_coord += ", %d" % (pre_coord * UNIT)
                    str_value += ", %d" % count
                    str_color += ", \"%s\"" % COLOR[count]
                count = 1
                pre_chrom = cur_chrom
                pre_coord = cur_coord
    if count:
        if is_first:
            str_chrom += "'%s')" % cur_chrom
            str_coord += "%d)" % (cur_coord * UNIT)
            str_value += "%d)" % count
            str_color += "\"%s\")" % COLOR[count]
            is_first = 0
        else:
            str_chrom += ", '%s')" % cur_chrom
            str_coord += ", %d)" % (cur_coord * UNIT)
            str_value += ", %d)" % count
            str_color += ",\"%s\")" % COLOR[count]

    ifd.close()
    print("%s" % str_chrom)
    print("%s" % str_coord)
    print("%s" % str_value)
    print("%s" % str_color)
    #print("max_value = %d" % max_count)

    return


def main(argv):
    ifile = ""
    sname = ""

    try:
        opts, args = getopt.getopt(argv, "hi:n:")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt == '-i':
            ifile = arg
        elif opt == '-n':
            sname = str(arg).replace("-", "_")

    if ifile == "":
        print("Error: '-i' is required")
        usage()
        sys.exit(2)
    elif not os.path.isfile(ifile):
        print("Error: input file(%s) is not existed" % ifile)
        usage()
        sys.exit(3)

    transform(ifile, sname)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
