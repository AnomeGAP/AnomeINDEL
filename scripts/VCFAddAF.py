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
# @file    VCFAddAF.py
#
# @brief   Add fake AF information into BALB_cJ.mgp.v5.all.sorted.vcf.gz
#
# @author  Chung-Tsai Su(chungtsai_su@atgenomix.com)
#
# @date    2020/01/02
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os
import gzip

# CONSTANT
AF = 0.5


def usage():
    print("VCFAddAF.py -i <Input vcf.gz file> -o <output vcf.gz file>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: input file (vcf.gz)")
    print("\t-o: output file (vcf.gz)")
    print("Usage:")
    print("\tpython ./VCFAddAF.py -i ../data/BALBcJ/BALB_cJ.mgp.v5.all.sorted.vcf.gz "
          " -o ../data/BALBcJ/BALB_cJ.mgp.v5.all.AF.sorted.vcf.gz")
    print("NOTE: need to repack the file by bgzip")
    print("root@cgbs:/seqslab/atsai/script# gzip -d /seqslab/atsai/data/BALBcJ/BALB_cJ.mgp.v5.all.AF.sorted.vcf.gz")
    print("root@cgbs:/seqslab/atsai/script# bgzip /seqslab/atsai/data/BALBcJ/BALB_cJ.mgp.v5.all.AF.sorted.vcf")
    print("root@cgbs:/seqslab/atsai/script# java -jar /usr/local/seqslab/third_party/gatk4.jar IndexFeatureFile -F /seqslab/atsai/data/BALBcJ/BALB_cJ.mgp.v5.all.AF.sorted.vcf.gz")
    return


def trimmer(ifile, ofile):
    ofd = gzip.open(ofile, "wb")
    ifd = gzip.open(ifile, "rt")
    idx = 0

    for line in ifd:
        if line.startswith("#"):
            ofd.write("%s" % line)
            if line.startswith("##INFO=<ID=DP4"):
                ofd.write("##INFO=<ID=AF,Number=G,Type=Float,Description=\"Faked by A-Tsai\">\n")
        else:
            items = line.strip().split("\t")
            num_alts = len(items[4].split(","))

            ofd.write("%s" % '\t'.join(items[:7]))
            af_str = "AF=%0.2f" % AF
            for i in range(1, num_alts):
                af_str += ",%0.2f" % AF
            ofd.write("\t%s;%s\n" % (af_str, '\t'.join(items[7:])))
            #print("%s\t%s\t%d\t%s" % (items[3], items[4], num_alts, af_str))
            #idx += 1
            #if idx > 1000:
            #    break

    ifd.close()
    ofd.close()
    print("%d variants are extracted" % idx)
    return


def main(argv):
    ofile = "output.vcf.gz"
    ifile = ""

    try:
        opts, args = getopt.getopt(argv, "hi:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt == '-i':
            ifile = arg
        elif opt == '-o':
            ofile = arg

    if ifile == "":
        print("Error: '-i' is required")
        usage()
        sys.exit(2)
    elif not os.path.isfile(ifile):
        print("Error: input file(%s) is not existed" % ifile)
        usage()
        sys.exit(3)

    trimmer(ifile, ofile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
