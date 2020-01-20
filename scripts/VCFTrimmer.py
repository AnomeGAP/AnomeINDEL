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
# @file    FASTQTrimmer.py
#
# @brief   Trimming several mouse strains to PoN for MuTect2. Input file is from
#          ftp://ftp.sra.ebi.ac.uk/vol1/ERZ022/ERZ022024/mgp.v3.indels.rsIDdbSNPv137.vcf.gz
#
# @author  Chung-Tsai Su(chungtsai_su@atgenomix.com)
#
# @date    2019/12/20
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
NUM_TRIMMING = 25
MIN_LENGTH = 60
MAX_LENGTH = 180


def usage():
    print("VCFTrimmer.py -i <Input vcf.gz file> -o <output vcf.gz file>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: input file (vcf.gz)")
    print("\t-o: output file (vcf.gz)")
    print("Usage:")
    print("\tpython ./VCFTrimmer.py -i ../data/WES/BALB_cJ/mgp.v3.indels.rsIDdbSNPv137.vcf.gz "
          " -o ../data/WES/BALB_cJ/BALB_cJ.GRCm38.vcf.gz")

    return


def trimmer(ifile, ofile):
    ofd = gzip.open(ofile, "wb")
    ifd = gzip.open(ifile, "rt")
    idx = 0

    for line in ifd:
        if not line.startswith("##"):
            # print("%s" % line)
            items = line.strip().split("\t")
            if not line.startswith("#"):
                if items[14] != ".":
                    alleles = items[4].split(",")
                    # print("%s" % alleles)
                    types = items[14].split(":")[0].split("/")
                    # print("%s" % types)
                    items[2] = "."
                    if types[0] == types[1]:
                        ofd.write("chr%s\t%s\t.\t.\t.\n" % ('\t'.join(items[:4]), alleles[int(types[0]) - 1]))
                        # print("chr%s\t%s\t.\t.\t." % ('\t'.join(items[:4]), alleles[int(types[0]) - 1]))
                    elif types[0] != items[4] and types[1] != items[4]:
                        ofd.write("chr%s\t%s,%s\t.\t.\t.\n" % ('\t'.join(items[:4]), alleles[int(types[0]) - 1], alleles[int(types[1]) - 1]))
                        # print("chr%s\t%s,%s\t.\t.\t." % ('\t'.join(items[:4]), alleles[int(types[0]) - 1], alleles[int(types[1]) - 1]))
                    elif types[0] != items[4]:
                        ofd.write("chr%s\t%s\t.\t.\t.\n" % ('\t'.join(items[:4]), alleles[int(types[0]) - 1]))
                        # print("chr%s\t%s\t.\t.\t." % ('\t'.join(items[:4]), alleles[int(types[0]) - 1]))
                    else:
                        ofd.write("chr%s\t%s\t.\t.\t.\n" % ('\t'.join(items[:4]), alleles[int(types[1]) - 1]))
                        # print("chr%s\t%s\t.\t.\t." % ('\t'.join(items[:4]), alleles[int(types[1]) - 1]))

                    idx += 1
            else:
                # CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	129P2	129S1	129S5	AJ	AKRJ	BALBcJ	C3HHeJ	C57BL6NJ	CASTEiJ	CBAJ	DBA2J	FVBNJ	LPJ	NODShiLtJNZOHlLtJ	PWKPhJ	SPRETEiJ	WSBEiJ
                ofd.write("%s\n" % '\t'.join(items[:8]))
                # print("%s" % '\t'.join(items[:8]))

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
