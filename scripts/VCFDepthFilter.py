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
MIN_DP = 100


def usage():
    print("VCFDepthFilter.py -i <Input vcf.gz file> -o <output vcf.gz file> -d <minimal depth>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: input file (vcf.gz)")
    print("\t-o: output file (vcf.gz)")
    print("\t-d: minal depth (default: %d)" % MIN_DP)
    print("Usage:")
    print("\tpython ./VCFDepthFilter.py -i ../data/WES/new-BALB-contam/Y02_Y01.BALB_cJ.contam.filtered2.vcf.gz"
          " -o ../data/WES/new-BALB-contam/Y02_Y01.BALB_cJ.contam.filtered2.d%d.vcf.gz -d %d" % (MIN_DP, MIN_DP))
    return


def depth(ifile, ofile, min_dp):
    ofd = gzip.open(ofile, "wb")
    ifd = gzip.open(ifile, "rt")
    idx = 0
    num_pass = 0

    for line in ifd:
        if line.startswith("#"):
            ofd.write("%s" % line)
        else:
            # 1       4351880 .       TA      T       .       artifact_in_normal;str_contraction      CONTQ=93;DP=141;ECNT=1;GERMQ=265;MBQ=30,30;MFRL=257,277;MMQ=60,60;MPOS=20;NALOD=-1.178e+00;NLOD=11.66;POPAF=6.00;RPA=11,10;RU=A;SAAF=0.051,0.051,0.065;SAPP=5.753e-03,0.021,0.973;STR;TLOD=5.48  GT:AD:AF:DP:F1R2:F2R1:OBAM:OBAMRC       0/0:56,2:0.049:58:28,1:28,1:false:false 0/1:72,5:0.075:77:25,3:47,2:false:false
            items = line.strip().split("\t")
            ad = int(items[9].split(":")[1].split(",")[1])
            dp1 = int(items[9].split(":")[3])
            dp2 = int(items[10].split(":")[3])

            if items[6] == "PASS":
                if dp1 >= min_dp and dp2 >= min_dp and ad == 0:
                    ofd.write("%s" % line)
                    #print("%s" % line)
                    num_pass += 1
                idx += 1

    ifd.close()
    ofd.close()
    print("%d/%d PASS variants are filtered" % (num_pass, idx))

    return


def main(argv):
    ofile = "output.vcf.gz"
    ifile = ""
    min_dp = MIN_DP
    try:
        opts, args = getopt.getopt(argv, "hi:o:d:")
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
        elif opt == '-d':
            min_dp = int(arg)

    if ifile == "":
        print("Error: '-i' is required")
        usage()
        sys.exit(2)
    elif not os.path.isfile(ifile):
        print("Error: input file(%s) is not existed" % ifile)
        usage()
        sys.exit(3)

    depth(ifile, ofile, min_dp)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
