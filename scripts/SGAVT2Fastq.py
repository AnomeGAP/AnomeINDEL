#!/bin/env python3
#
# @note Copyright (C) 2016, Anome Incorporated. All Rights Reserved.
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
# @file    VCFQualityFilter.py
#
# @brief   An example to filter Quality and Coverage on VCF file (Delly)
#
# @author  A-Tsai Su(chungtsai_su@anome.com.tw)
#
# @date    2016/01/23
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os.path


def Usage():
    ''' Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    '''
    print("SGAVT2Fastq.py -i <SGA VT file> -o <FASTQ output file>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <SGA VT file>")
    print("\t-o: <SGA VT file>")
    print("Usage:")
    print("\t./SGAVT2Fastq.py -i testing1.vt -o testing1.fq")

    return


def SGAVT2Fastq(ifn, ofn):

    ofd = open(ofn, "w+")
    ifd = open(ifn, "r")
    ##VT    123     ATCG    ss:01
    for line in ifd:
        line = line.strip()
        items = line.split('\t')
        ofd.write("@%s\n" % items[1])
        ofd.write("%s\n" % items[2])
        ofd.write("+%s\n" % items[1])
        ofd.write("%s\n" % ("F" * len(items[2])))

    ifd.close()
    ofd.close()

    return


def main(argv):
    infile = ""
    outfile = ""

    try:
        opts, args = getopt.getopt(argv,"hi:o:")
    except getopt.GetoptError:
        Usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            Usage()
            sys.exit()
        elif opt in "-i":
            infile = arg
            outfile = infile + ".fq"
        elif opt in "-o":
            outfile = arg

    if not infile:
        Usage()
        sys.exit(2)
    elif not os.path.isfile(infile):
        print("Error: input file(%s) is not existed" % infile)
        Usage()
        sys.exit(3)

    #Main Function
    SGAVT2Fastq(infile, outfile)

    return

if __name__ == '__main__':
    main(sys.argv[1:])
