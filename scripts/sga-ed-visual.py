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
# @file    sga-ed-visual.py
#
# @brief   An example to transform SGA ED record to graph record for Cytoscope
#
# @author  A-Tsai Su(chungtsai_su@anome.com.tw)
#
# @date    2018/05/16
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
    print("sga-ed-visual.py -i <SGA ED file> -o <SGA ED visual file>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <SGA ED file>")
    print("\t-0: <SGA ED visual file>")
    print("Usage:")
    print("\t./sga-ed-visual.py -i m105.ed.txt -o m105.ed.tsv")

    return


def sga_ed_visual(ifn, ofn):

    ofd = open(ofn, "w")
    ifd = open(ifn, "r")
    ##ED    123 456 100 150 151
    for line in ifd:
        line = line.strip()
        items = line.split('\t')
        res = items[1].split(' ')
        ofd.write("%s\t%s\t%d\n" % (res[0], res[1], abs(int(res[2]) - int(res[3]))+ 1))

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
    sga_ed_visual(infile, outfile)

    return

if __name__ == '__main__':
    main(sys.argv[1:])
