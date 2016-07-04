#!/bin/env python
#
# @note Copyright (C) 2015, Anome Incorporated. All Rights Reserved.
#       This program is an unpublished copyrighted work which is proprietary to
#       Anome Incorporated and contains confidential information that is not to
#       be reproduced or disclosed to any other person or entity without prior
#       written consent from Anome, Inc. in each and every instance.
#
# @warning Unauthorized reproduction of this program as well as unauthorized
#          preparation of derivative works based upon the program or distribution of
#          copies by sale, rental, lease or lending are violations of federal copyright
#          laws and state trade secret laws, punishable by civil and criminal penalties.
#
# @file    ExtractUnMappedReads.py
#
# @brief   Extracting unmapped reads from BAM file
#
# @author  Chung-Tsai Su(chungtsai_su@anome.com)
#
# @date    2015/03/26
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os.path
import math
from collections import defaultdict
import pysam    ##http://pysam.readthedocs.org/en/latest/api.html

#Default Parameter
MIN_SEGMENTLENGTH = 500
MAX_SEGMENTLENGTH = 50000
MIN_QUALITY = 10
MIN_COVERAGE = 0.9
MIN_SUPPORT = 5


def Usage():
    print("ExtractUnMappedReads.py -i <Input BAM> -o <Output SAM>") 
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: Input BAM file")
    print("\t-o: Output SAM file which contains all unmapped reads[Default=<inputfile>.unmapped.sam]")
    print("Usage:")
    print("\t./ExtractUnMappedReads.py -i ../example/example.bam -o ../example/example.unmapped.sam")

    return

def ExtractUnMappedReads(ifn, ofn):
    samfile = pysam.AlignmentFile(ifn, "rb")
    output = pysam.AlignmentFile(ofn, "wh", template=samfile)

    print("Number of mapped reads = %d" % (samfile.mapped))
    print("Number of unmapped reads = %d" % (samfile.unmapped))
    print("Unmapped ratio = %.2f%%" % ((float)(samfile.unmapped*100)/(samfile.mapped+samfile.unmapped) ))
    print("nocoordinate=%d" % (samfile.nocoordinate))
    print("")
    num_reads = 0
    num_unmapped = 0

    for read in samfile.fetch():
        if (read.is_unmapped == 1):
            output.write(read)
            num_unmapped += 1

        num_reads += 1

    print("Number of reads = %d, Number of mapped reads=%d, Number of unmapped reads=%d (%.2f%%)" % ( num_reads, num_reads - num_unmapped, num_unmapped, (float)(num_unmapped*100)/num_reads ))

    output.close()
    samfile.close()
    return

def main(argv):
    ifile = ""
    ofile = ""

    try:
        opts, args = getopt.getopt(argv,"hi:o:")
    except getopt.GetoptError:
        Usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            Usage()
            sys.exit()
        elif opt in ("-i"):
            ifile = arg
        elif opt in ("-o"):
            ofile = arg

    #error handling for input parameters
    if (ifile == ""):
        print("Error: '-i' is required")
        Usage()
        sys.exit(2)
    elif (os.path.isfile(ifile) == False): 
        print("Error: input file(%s) is not existed" % (ifile))
        Usage()
        sys.exit(3)
    if (ofile == ""):
        ofile = "%s.unmapped.sam" % (ifile)

    #Main Function
    ExtractUnMappedReads(ifile, ofile)

    return

if __name__ == '__main__':
    main(sys.argv[1:])

