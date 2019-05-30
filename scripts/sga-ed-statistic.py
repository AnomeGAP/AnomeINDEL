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
# @file    sga-ed-statistic.py
#
# @brief   An example to calculate the distribution of SGA ED
#
# @author  A-Tsai Su(chungtsai_su@anome.com.tw)
#
# @date    2018/05/18
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os.path
from collections import defaultdict

MAX_VALUE = 100
RC_VALUE = 10000000000


def Usage():
    """ Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    """
    print("sga-ed-statistic.py -i <SGA ED file>  -o <SGA ED distribution file>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <SGA ED file>")
    print("\t-0: <SGA ED distribution file>")
    print("Usage:")
    print("\t.python ./sga-ed-statistic.py -i m85.graphseq.ed.txt -o m85.graphseq.ed.dist.log")

    return


def sga_ed_statistic(ifn, ofn):
    total = 0
    hiread = defaultdict(int)
    horead = defaultdict(int)
    ifd = open(ifn, "r")
    # #ED    123 456 100 150 151 0 50 151 0 0
    for line in ifd:
        line = line.strip()
        items = line.split('\t')
        res = items[1].split(' ')
        if res[2] == '0':
            if res[5] == '0':
                horead[RC_VALUE+int(res[0])] += 1
                hiread[int(res[1])] += 1
                horead[RC_VALUE + int(res[1])] += 1
                hiread[int(res[0])] += 1
            else:
                horead[RC_VALUE+int(res[0])] += 1
                hiread[RC_VALUE+int(res[1])] += 1
                horead[int(res[1])] += 1
                hiread[int(res[0])] += 1
        else:
            if res[5] == '0':
                horead[int(res[0])] += 1
                hiread[int(res[1])] += 1
                horead[RC_VALUE+int(res[1])] += 1
                hiread[RC_VALUE+int(res[0])] += 1
            else:
                horead[int(res[0])] += 1
                hiread[RC_VALUE+int(res[1])] += 1
                horead[int(res[1])] += 1
                hiread[RC_VALUE + int(res[0])] += 1
        total += 2
    ifd.close()
    print("Total %d edges" % total)

    matrix = [[0 for x in range(MAX_VALUE+1)] for y in range(MAX_VALUE+1)]
    for i in hiread:
        # print("%d\t%d" % (i, hread[i]))
        x = hiread[i]
        if hiread[i] > MAX_VALUE:
            x = MAX_VALUE
        else:
            x = hiread[i]
        y = horead[i]
        if horead[i] > MAX_VALUE:
            y = MAX_VALUE
        else:
            y = horead[i]
        matrix[x][y] += 1

    for i in horead:
        if hiread[i] == 0:
            if horead[i] > MAX_VALUE:
                y = MAX_VALUE
            else:
                y = horead[i]
            matrix[0][y] += 1

    ofd = open(ofn, "w")
    for i in range(MAX_VALUE + 1):
        for j in range(MAX_VALUE + 1):
            ofd.write("%d\t" % matrix[i][j])
        ofd.write("\n")
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

    # Main Function
    sga_ed_statistic(infile, outfile)

    return

if __name__ == '__main__':
    main(sys.argv[1:])
