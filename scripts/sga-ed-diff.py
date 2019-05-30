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
    print("SGAVT2Fastq.py -i <GraphSeq ED sorted file> -j <SGA ED ordered and sorted file>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <GraphSeq ED sorted file>")
    print("\t-j: <SGA ED ordered and sorted file>")
    print("Usage:")
    print("\t./sga-ed-diff.py -i all.ed.sorted -j all-sga.ed.ordered.sorted")

    return


def sga_ed_diff(ifn1, ifn2):

    ifd1 = open(ifn1, "r")
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
    in1file = ""
    in2file = ""

    try:
        opts, args = getopt.getopt(argv,"hi:j:")
    except getopt.GetoptError:
        Usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            Usage()
            sys.exit()
        elif opt in "-i":
            in1file = arg
        elif opt in "-j":
            in2file = arg

    if not infile:
        Usage()
        sys.exit(2)
    elif not os.path.isfile(in1file):
        print("Error: input file(%s) is not existed" % in1file)
        Usage()
        sys.exit(3)
    elif not os.path.isfile(in2file):
        print("Error: input file(%s) is not existed" % in2file)
        Usage()
        sys.exit(4)

    #Main Function
    sga_ed_diff(in1file, in2file)

    return

if __name__ == '__main__':
    main(sys.argv[1:])
