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
# @file    CliqueStatistic.py
#
# @brief   A program to statistic clique data
#
# @author  A-Tsai Su(chungtsai_su@anome.com.tw)
#
# @date    2016/04/18
#
# @version 1.0
#
# @remark
#

import sys
import getopt
import os.path
import re


def Usage():
    """ Description: Program Usage
        Argument:    NONE
        Return:	     NONE
    """
    print("TransfacStatistic.py -i <TF2SITE listing file> -j <SITE2TF listing file>")
    print("Argument:")
    print("\t-h: Usage")
    print("\t-i: <TF2SITE list file>")
    print("\t-j: <TF2SITE list file>")
    print("Usage:")
    print("\tpython ~/src/github/AnomeINDEL/scripts/TransfacStatistic.py -i ~/src/clique-tf/output/TF2SITE.listing.txt -j ~/src/clique-tf/output/SITE2TF.listing.txt")

    return

MAX_COUNT = 20


def TransfacStatistic(ifn1, ifn2):
    _num = 0
    _max_neighbors = 0
    _max_neighbor_name = ""
    _hcount = [0 for x in range(MAX_COUNT+1)]

    print(ifn1)
    ifd = open(ifn1, "r")
    # Read sequence and segment them with length of slength
    for line in ifd:
        if re.match("^#", line):
            # print "Skip comment: %s" % (line)
            continue

        line = line.strip()
        items = line.split('\t')

        lRight = items[1].split(';')
        num_right = len(lRight) - 1

        if num_right > _max_neighbors:
            _max_neighbors = num_right
            _max_neighbor_name = items[0]

        if num_right > MAX_COUNT:
            num_right = MAX_COUNT

        _hcount[num_right] += 1
        _num += 1

    ifd.close()

    # output
    print("Number of node = %d" % _num)
    print("Maximal neighbors=%d %s" % (_max_neighbors,_max_neighbor_name))
    print("%s" % "\n".join([str(x) for x in _hcount]))

    _num = 0
    _max_neighbors = 0
    _hcount = [0 for x in range(MAX_COUNT+1)]
    print(ifn2)
    ifd = open(ifn2, "r")
    # Read sequence and segment them with length of slength
    for line in ifd:
        if re.match("^#", line):
            # print "Skip comment: %s" % (line)
            continue

        line = line.strip()
        items = line.split('\t')

        lRight = items[1].split(';')
        num_right = len(lRight) - 1

        if num_right > _max_neighbors:
            _max_neighbors = num_right
            _max_neighbor_name = items[0]

        if num_right > MAX_COUNT:
            num_right = MAX_COUNT

        _hcount[num_right] += 1
        _num += 1

    ifd.close()

    # output
    print("Number of node = %d" % _num)
    print("Maximal neighbors=%d %s" % (_max_neighbors,_max_neighbor_name))
    print("%s" % "\n".join([str(x) for x in _hcount]))

    return


def main(argv):
    infile1 = ""

    try:
        opts, args = getopt.getopt(argv, "hi:j:")
    except getopt.GetoptError:
        Usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            Usage()
            sys.exit()
        elif opt in "-i":
            infile1 = arg
        elif opt in "-j":
            infile2 = arg
    if not infile1 or not infile2:
        Usage()
        sys.exit(2)
    elif not os.path.isfile(infile1):
        print("Error: input file(%s) is not existed" % infile1)
        Usage()
        sys.exit(3)
    elif not os.path.isfile(infile2):
        print("Error: input file(%s) is not existed" % infile2)
        Usage()
        sys.exit(4)

    # Main Function
    TransfacStatistic(infile1, infile2)

    return


if __name__ == '__main__':
    main(sys.argv[1:])
