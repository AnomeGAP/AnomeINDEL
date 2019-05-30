#!/bin/env python
#
# @note Copyright (C) 2018, Atgenomix Incorporated. All Rights Reserved.
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
# @file    BAMPileupScan.py
#
# @brief   Scan the coverage of each region in genome by pileup
#
# @author  Chung-Tsai Su(chungtsai_su@atgenomix.com)
#
# @date    2019/02/14
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os.path
import pysam  # http://pysam.readthedocs.org/en/latest/api.html


# CONST
chrID = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
         'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
         'chr21', 'chr22', 'chrX', 'chrY']
chrSize = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422,
           135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167,
           46709983, 50818468, 156040895, 57227415]

# Default Parameter


def usage():
    print("BAMPileupScan.py -i <Input BAM> -w <Window Size> -o <Non-mapped region>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: Input BAM file")
    print("\t-o: Output depth report")
    print("Usage:")
    print("\ttime python ./BAMPileupScan.py -i ~/NA12878-novaseq-chr6.bam -o ~/NA12878-novaseq-chr6.non-mapping-region.log")

    return


def pileupScan(ifn, ofn):

    samfile = pysam.AlignmentFile(ifn, "rb")
    total = 0
    total_mapped = 0

    # print("%d" % len(list(samfile.pileup("chr6", 5341131,5341353))))

    ofd = open(ofn, "w")
    for i in range(0, len(chrID)):
        # total += chrSize[i]
        # mapped = 0
        # for pileupcolumn in samfile.pileup(chrID[i], 0, chrSize[i]-1):
        #     # print("%s\t%s" % (pileupcolumn.pos, pileupcolumn.n))
        #     mapped += 1
        # print("%s\t%d\t%d\t%2.1f%%" % (chrID[i], mapped, chrSize[i], (chrSize[i] - mapped)*100/chrSize[i]))
        # total_mapped += mapped
        total += chrSize[i]
        mapped = 0
        prev = -1
        segment = 0
        for pileupcolumn in samfile.pileup(chrID[i], 0, chrSize[i]-1):
            # print("%s\t%s" % (pileupcolumn.pos, pileupcolumn.n))
            if prev+1 == pileupcolumn.pos:
                mapped += 1
            else:
                ofd.write("%s\t%s\t%s\t%d\n" % (chrID[i], prev+1, pileupcolumn.pos, pileupcolumn.pos-prev))
                segment += 1
            prev = pileupcolumn.pos
        print("%s\t%d\t%d\t%2.1f%%\t%d" % (chrID[i], mapped, chrSize[i], (chrSize[i] - mapped)*100/chrSize[i], segment))
        total_mapped += mapped

    ofd.close()
    print("Total\t%d\t%d\t%2.1f%%" % (total_mapped, total, (total - total_mapped)*100/total))
    samfile.close()

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

    # error handling for input parameters
    if ifile == "":
        print("Error: '-i' is required")
        usage()
        sys.exit(2)
    elif not os.path.isfile(ifile):
        print("Error: input file(%s) is not existed" % ifile)
        usage()
        sys.exit(3)

    if ofile == "":
        ofile = "%s.non-mapped-region.log" % ifile

    # Main Function
    pileupScan(ifile, ofile)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
